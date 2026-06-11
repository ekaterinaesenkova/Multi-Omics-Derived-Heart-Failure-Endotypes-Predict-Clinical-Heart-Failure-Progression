################################################################################
# Multi-Omics-Derived Heart Failure Endotypes Predict Clinical Heart Failure Progression (MyoVasc Study)
#
# Description:
#   This script performs multi-omics integration (proteomics + lipidomics) using
#   Similarity Network Fusion (SNF) to identify patient subgroups in Stage C/D
#   heart failure. It evaluates clustering quality via survival analysis
#   (Harrell's C-index, log-rank tests), characterises clusters by clinical and
#   omics features, trains Random Forest classifiers, and validates cluster
#   assignments in an independent asymptomatic cohort.
#
# Study:       MyoVasc cohort
# Outcome:     Worsening heart failure (WHF), time-to-event
# Data:        Baseline visit proteomics (Olink) + lipidomics
#
# Authors:     Esenkova Ekaterina, Elisa Araldi, Elena Casiraghi
# Institution: Unimedizin Mainz
# Date:        11.06.2026
#
# Repository:  https://github.com/ekaterinaesenkova/Multi-Omics-Derived-Heart-Failure-Endotypes-Predict-Clinical-Heart-Failure-Progression
# Preprint/DOI: https://www.medrxiv.org/content/10.1101/2025.01.28.25321241v1
#
# NOTE: Raw data cannot be shared due to participant privacy. Contact the corresponding authors
# for data access requests.
#
# Usage:
#   1. Set PATHS in Section 0 to match your local environment.
#   2. Run sections sequentially (0 → 1 → 2 → ... → 9).
#   3. Section 5 onwards requires objects produced by earlier sections.
################################################################################


# ==============================================================================
# SECTION 0: Configuration — set all paths and global parameters here
# ==============================================================================

# ---- 0.1 User-adjustable paths -----------------------------------------------
# Replace these with your own paths before running.

PATH_WORKING  <- "path/to/your/working/directory"   # main working directory
PATH_DATA_SQL <- "path/to/sql/data/scripts"          # SQL data-loading scripts
PATH_R_FUNCS  <- "path/to/shared/R/functions"        # shared helper functions
PATH_SAVE     <- file.path(PATH_WORKING, "outputs")  # where figures/CSVs go

# ---- 0.2 Analysis parameters -------------------------------------------------
K_MIN            <- 3      # minimum number of clusters to test
K_MAX            <- 10     # maximum number of clusters to test
K_FINAL          <- 8      # final chosen cluster number (from C-index / NMI results)
CORR_CUTOFF      <- 0.9    # Pearson r cutoff for removing collinear features
NZV_FREQRATIO    <- 95/5   # near-zero-variance frequency ratio (passed to caret)
N_BOOTSTRAP      <- 100    # bootstrap replicates for stability analysis
BOOTSTRAP_FRAC   <- 0.80   # fraction of samples per bootstrap replicate
RF_NTREE         <- 1000   # number of trees in Random Forest models
RF_PROB_CUTOFF   <- 0.25   # minimum RF probability to assign a cluster label
MOFA_FACTORS     <- 10     # number of MOFA latent factors
MOFA_MAXITER     <- 1000   # maximum MOFA training iterations
RANDOM_SEED      <- 123    # global random seed for reproducibility

set.seed(RANDOM_SEED)

# Cluster colour palette (8 clusters)
CLUSTER_COLORS <- c(
  "1" = "#FF3355", "2" = "#C4961A", "3" = "#56B4E9", "4" = "#FF9900",
  "5" = "#0000FF", "6" = "#33FF66", "7" = "#293352", "8" = "#FFDB6D"
)


# ==============================================================================
# SECTION 1: Libraries
# ==============================================================================

# CRAN packages
required_packages <- c(
  "SNFtool", "tidyr", "ComplexHeatmap", "survival", "survminer",
  "flexclust", "ggplot2", "Rtsne", "rstatix", "stringr", "janitor",
  "ggsurvfit", "readxl", "caret", "clv", "tidycmprsk", "Hmisc",
  "aricode", "ggalluvial", "patchwork", "pheatmap", "randomForest",
  "mclust", "cmprsk", "kernlab", "Cairo", "circlize", "grid",
  "glmnet", "dplyr", "broom", "cowplot", "reshape2", "tidytext",
  "survival", "compareC", "labelled"
)

# Install missing packages automatically
new_pkgs <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

invisible(lapply(required_packages, library, character.only = TRUE))

# Bioconductor packages (install once with BiocManager::install("MOFA2"))
if (!requireNamespace("MOFA2", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("MOFA2")
}
library(MOFA2)


# ==============================================================================
# SECTION 2: Helper Functions
# ==============================================================================

#' Impute missing values in a numeric data frame
#' (wrapper around your custom imputa() function loaded from Imputa.R)
#'
#' @param df Numeric data frame (samples × features)
#' @return Data frame with imputed values
impute_data <- function(df) {
  imputa(df)
}


#' Remove near-zero-variance and highly correlated features
#'
#' @param data_imputed Numeric matrix/data.frame (samples × features), already imputed
#' @param corr_cutoff  Pearson r above which one of a pair of features is dropped (default 0.9)
#' @return List with:
#'   \item{data_norm}{Scaled matrix with problematic features removed}
#'   \item{train_mean}{Column means used for scaling (save for test-set projection)}
#'   \item{train_sd}{Column SDs used for scaling}
preprocess_omics <- function(data_imputed, corr_cutoff = CORR_CUTOFF) {
  
  # Remove near-zero-variance features
  nzv <- caret::nearZeroVar(data_imputed, names = TRUE)
  if (length(nzv) > 0) {
    message("  Removing ", length(nzv), " near-zero-variance features.")
    data_imputed <- data_imputed[, !colnames(data_imputed) %in% nzv]
  }
  
  # Z-score normalisation
  train_mean <- apply(data_imputed, 2, mean, na.rm = TRUE)
  train_sd   <- apply(data_imputed, 2, sd,   na.rm = TRUE)
  data_norm  <- scale(data_imputed, center = train_mean, scale = train_sd)
  
  # Remove highly correlated features
  cc          <- cor(data_norm, use = "pairwise.complete.obs", method = "pearson")
  select_corr <- caret::findCorrelation(cc, cutoff = corr_cutoff, exact = FALSE)
  if (length(select_corr) > 0) {
    message("  Removing ", length(select_corr), " highly correlated features.")
    data_norm <- data_norm[, -select_corr]
  }
  
  list(data_norm  = data_norm,
       train_mean = train_mean,
       train_sd   = train_sd)
}


#' Load and extract one omics panel from the MyoVasc SQL data object
#'
#' @param dall1       Data frame returned by read_SQL_MyoVasc_BL()
#' @param var_attr    Attribute name for variable codes (e.g. "lipid_vars_raw")
#' @param nam_attr    Attribute name for display names  (e.g. "lipid_nams_raw")
#' @param id_col      Name of the patient ID column
#' @param excl_cols   Optional character vector of column names to exclude
#' @return Data frame (samples × features) with rownames = patient IDs
extract_omics_panel <- function(dall1, var_attr, nam_attr,
                                id_col    = "v11_sid01",
                                excl_cols = NULL) {
  dattr     <- attributes(dall1)
  vars      <- dattr[[var_attr]]
  nams      <- dattr[[nam_attr]]
  panel     <- dall1[, colnames(dall1) %in% vars, drop = FALSE]
  panel     <- data.frame(dall1[[id_col]], panel)
  colnames(panel) <- c("id", nams)
  rownames(panel) <- panel$id
  panel     <- panel[, -1]
  
  # Replace NA in lipid-style panels with 0 (consistent with original code)
  panel[is.na(panel)] <- 0
  
  # Drop excluded columns (e.g. duplicate BNP measures in protein panel)
  if (!is.null(excl_cols)) {
    excl_found <- excl_cols[excl_cols %in% colnames(panel)]
    if (length(excl_found) > 0) panel <- panel[, !colnames(panel) %in% excl_found]
  }
  panel
}


#' Intersect samples across omics panels and restrict to a given patient subset
#'
#' @param panels_list List of data frames (each with rownames = patient IDs)
#' @param id_subset   Character vector of patient IDs to keep (e.g. Stage C/D)
#' @return List of data frames, each restricted to the common intersection
intersect_panels <- function(panels_list, id_subset = NULL) {
  # Keep only complete-case rows in each panel
  panels_cc <- lapply(panels_list, function(p) p[complete.cases(p), ])
  
  # Intersect patient IDs across all panels
  common_ids <- Reduce(intersect, lapply(panels_cc, rownames))
  
  # Optionally restrict to a clinical subgroup
  if (!is.null(id_subset)) {
    common_ids <- intersect(common_ids, id_subset)
  }
  
  lapply(panels_cc, function(p) p[rownames(p) %in% common_ids, ])
}


#' Build per-panel affinity matrices for SNF
#'
#' @param panels_list_norm List of pre-processed (scaled, filtered) matrices
#' @param K_affinity       Number of nearest neighbours for affinity matrix
#' @param sigma            Bandwidth parameter for affinity matrix
#' @return List of affinity matrices (one per omics panel)
build_affinity_matrices <- function(panels_list_norm,
                                    K_affinity = 20,
                                    sigma      = 0.5) {
  lapply(panels_list_norm, function(mat) {
    dist_mat <- (SNFtool::dist2(as.matrix(mat), as.matrix(mat)))^(1/2)
    W        <- SNFtool::affinityMatrix(dist_mat, K = K_affinity, sigma = sigma)
    # Ensure diagonal equals row maximum (standard SNF normalisation)
    W <- W - diag(apply(W, 1, max))
    W <- W + diag(apply(W, 1, max))
    W
  })
}


#' Compute spectral cluster memberships for a range of K
#'
#' @param W          Fused (or single-panel) affinity matrix
#' @param k_min      Minimum K
#' @param k_max      Maximum K
#' @return Data frame (samples × K_range), columns named "mem<K>"
compute_memberships <- function(W, k_min = K_MIN, k_max = K_MAX) {
  mem_list <- lapply(k_min:k_max, function(k) {
    SNFtool::spectralClustering(W, K = k)
  })
  mem_df <- as.data.frame(do.call(cbind, mem_list))
  colnames(mem_df) <- paste0("mem", k_min:k_max)
  mem_df
}


#' Rank cluster labels by median survival time (ascending = higher risk first)
#'
#' Re-labels clusters 1..K so that cluster 1 has the worst survival and
#' cluster K has the best, matching the convention used in the manuscript.
#'
#' @param membership Integer vector of cluster labels
#' @param time       Numeric vector of survival times
#' @param event      Numeric/integer event indicator (1 = event)
#' @return Integer vector with re-ordered labels
rank_memberships_by_survival <- function(membership, time, event) {
  model_fit <- survival::survfit(
    survival::Surv(time, event) ~ membership
  )
  median_surv <- summary(model_fit)$table[, "median"]
  old_labels  <- seq_along(median_surv)
  # rank so that shortest median → label 1
  new_labels  <- rank(median_surv, ties.method = "first")
  x           <- as.integer(membership)
  new_x       <- new_labels[match(x, old_labels)]
  new_x
}


#' Compute Harrell's C-index and log-rank p-value for a cluster membership vector
#'
#' @param surv_time   Numeric vector of event/censoring times
#' @param surv_event  Integer event indicator
#' @param cluster_vec Cluster membership vector (numeric or factor)
#' @return Named numeric vector: k, log10_p, c_index
compute_survival_metrics <- function(surv_time, surv_event, cluster_vec) {
  surv_obj <- survival::Surv(time = surv_time, event = surv_event)
  
  # Log-rank test
  lr       <- survival::survdiff(surv_obj ~ cluster_vec)
  df_lr    <- length(lr$n) - 1
  pval     <- stats::pchisq(as.numeric(lr$chisq), df = df_lr, lower.tail = FALSE)
  log10_p  <- -log10(max(pval, .Machine$double.xmin))
  
  # Harrell's C-index via Cox linear predictor
  cox_mod     <- survival::coxph(surv_obj ~ cluster_vec)
  risk_scores <- stats::predict(cox_mod, type = "lp")
  c_idx       <- as.numeric(Hmisc::rcorr.cens(-risk_scores, surv_obj)[1])
  
  c(log10_p = log10_p, c_index = c_idx)
}


#' Evaluate survival metrics across a range of K for one membership data frame
#'
#' @param memberships_df Data frame with columns "mem<K>" for K in k_min:k_max
#' @param surv_time      Numeric survival time vector (aligned with membership rows)
#' @param surv_event     Integer event indicator
#' @param k_min          Minimum K tested
#' @param k_max          Maximum K tested
#' @return Data frame with columns: k, log10_p, c_index
evaluate_clustering_survival <- function(memberships_df,
                                         surv_time, surv_event,
                                         k_min = K_MIN, k_max = K_MAX) {
  results <- do.call(rbind, lapply(k_min:k_max, function(k) {
    clust_col  <- paste0("mem", k)
    metrics    <- compute_survival_metrics(surv_time, surv_event,
                                           memberships_df[[clust_col]])
    data.frame(k = k, log10_p = metrics["log10_p"], c_index = metrics["c_index"],
               row.names = NULL)
  }))
  results
}


#' Run pairwise Wilcoxon tests for a single feature and symmetrise the p-matrix
#'
#' @param feature_vec Numeric vector of feature values
#' @param cluster_vec Cluster membership vector
#' @param p_thresh    Significance threshold
#' @return Integer vector: number of significantly different cluster pairs for each cluster
pairwise_wilcox_sig <- function(feature_vec, cluster_vec, p_thresh = 0.001) {
  make_symm <- function(m) {
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    m
  }
  wt  <- stats::pairwise.wilcox.test(feature_vec, cluster_vec,
                                     p.adjust.method = "BH")
  tab <- wt$p.value
  n   <- max(as.integer(cluster_vec), na.rm = TRUE)
  # Expand to full n×n matrix with NA diagonal
  tab_full <- as.matrix(cbind(rbind(rep(NA, n), tab), rep(NA, n)))
  tab_full <- make_symm(tab_full)
  tab_bin  <- ifelse(tab_full < p_thresh, 1, 0)
  colSums(tab_bin, na.rm = TRUE)
}


#' Run full SNF pipeline on a (sub)set of omics data
#'
#' Used both for the main analysis and for bootstrap stability testing.
#'
#' @param data_list    List of numeric matrices (samples × features), one per panel
#' @param K_cluster    Number of spectral clusters
#' @param K_affinity   Number of nearest neighbours for affinity
#' @param sigma        Bandwidth for affinity matrix
#' @return Integer vector of cluster labels
run_snf_pipeline <- function(data_list,
                             K_cluster  = K_FINAL,
                             K_affinity = 20,
                             sigma      = 0.5) {
  W_list <- lapply(data_list, function(mat) {
    data_imp  <- imputa(mat)
    data_norm <- scale(data_imp)
    dist_mat  <- (SNFtool::dist2(as.matrix(data_norm), as.matrix(data_norm)))^(1/2)
    SNFtool::affinityMatrix(dist_mat, K = K_affinity, sigma = sigma)
  })
  W_fused  <- SNFtool::SNF(W_list, K = K_affinity)
  clusters <- SNFtool::spectralClustering(W_fused, K = K_cluster)
  clusters
}


#' Bootstrap stability analysis: ARI and NMI over B subsamples
#'
#' @param data_list     List of pre-processed (scaled, filtered) matrices
#' @param reference_mem Integer reference cluster membership for the full cohort
#' @param B             Number of bootstrap replicates
#' @param frac          Fraction of samples per replicate
#' @param K_cluster     Number of clusters
#' @return Data frame with columns: b (replicate), ari, nmi
bootstrap_stability <- function(data_list, reference_mem,
                                B          = N_BOOTSTRAP,
                                frac       = BOOTSTRAP_FRAC,
                                K_cluster  = K_FINAL) {
  n <- nrow(data_list[[1]])
  results <- do.call(rbind, lapply(seq_len(B), function(b) {
    idx      <- sample(n, size = floor(frac * n), replace = FALSE)
    data_sub <- lapply(data_list, function(m) m[idx, ])
    cl_sub   <- run_snf_pipeline(data_sub, K_cluster = K_cluster)
    cl_ref   <- reference_mem[idx]
    data.frame(
      b   = b,
      ari = mclust::adjustedRandIndex(cl_ref, cl_sub),
      nmi = aricode::NMI(cl_ref, cl_sub)
    )
  }))
  results
}


#' Train a binary Random Forest for one cluster vs. rest
#'
#' @param train_matrix  Numeric matrix of training features (samples × features)
#' @param cluster_label Integer vector of cluster memberships (training cohort)
#' @param target_k      Cluster of interest (1 vs. rest)
#' @param ntree         Number of trees
#' @return Trained randomForest object
train_rf_one_vs_rest <- function(train_matrix, cluster_label,
                                 target_k, ntree = RF_NTREE) {
  df         <- as.data.frame(train_matrix)
  colnames(df) <- make.names(colnames(df))
  df$cluster <- as.factor(ifelse(cluster_label == target_k, 1, 0))
  randomForest::randomForest(cluster ~ ., data = df, ntree = ntree, importance = TRUE)
}


#' Assign samples to clusters using a list of one-vs-rest RF models
#'
#' @param rf_models    Named list of RF models (one per cluster 1..K)
#' @param new_data     Numeric matrix (samples × features) — test cohort
#' @param prob_cutoff  Minimum probability to assign; samples below get label K+1
#' @return Integer vector of assigned cluster labels (K+1 = "unassigned")
assign_clusters_rf <- function(rf_models, new_data, prob_cutoff = RF_PROB_CUTOFF) {
  new_df  <- as.data.frame(new_data)
  colnames(new_df) <- make.names(colnames(new_df))
  
  prob_mat <- do.call(cbind, lapply(rf_models, function(m) {
    preds <- predict(m, newdata = new_df, type = "prob")
    preds[, "1"]  # probability of belonging to target cluster
  }))
  colnames(prob_mat) <- names(rf_models)
  
  apply(prob_mat, 1, function(row) {
    if (max(row) > prob_cutoff) which.max(row) else length(rf_models) + 1
  })
}


#' Plot C-index vs. number of clusters for multiple methods
#'
#' @param c_index_df Data frame with columns: k, and one column per method
#' @param method_labels Named character vector mapping column names to display labels
#' @return ggplot object
plot_c_index <- function(c_index_df,
                         method_labels = c(Lipids   = "Lipids",
                                           Proteins = "Proteins",
                                           SNF      = "SNF (Lipids + Proteins)")) {
  df_long <- tidyr::pivot_longer(c_index_df, -k,
                                 names_to  = "Method",
                                 values_to = "C_index")
  df_long$Method <- factor(df_long$Method,
                           levels = names(method_labels),
                           labels = method_labels)
  means <- df_long |>
    dplyr::group_by(Method) |>
    dplyr::summarize(mean_c = mean(C_index), .groups = "drop")
  
  ggplot2::ggplot(df_long, ggplot2::aes(x = k, y = C_index, color = Method)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(data = means,
                        ggplot2::aes(yintercept = mean_c, color = Method),
                        linetype = "dashed", linewidth = 0.7, show.legend = FALSE) +
    ggplot2::scale_x_continuous(breaks = unique(df_long$k)) +
    ggplot2::labs(x = "Number of clusters", y = "C-index", color = "Data") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = element_blank(),
                   legend.position = "bottom")
}


#' Notched boxplot with significance stars vs. a reference cluster
#'
#' @param df          Data frame containing the feature and membership columns
#' @param y_var       Name of the feature column (string)
#' @param membership  Name of the membership column (string)
#' @param title       Plot title
#' @param y_label     Y-axis label
#' @param colors      Named character vector mapping cluster labels to colours
#' @param ref_group   Reference group for pairwise t-test (default "8")
#' @return ggplot object
plot_cluster_boxplot <- function(df, y_var, membership = "membership",
                                 title, y_label,
                                 colors    = CLUSTER_COLORS,
                                 ref_group = "8") {
  ggplot2::ggplot(df, ggplot2::aes(x = .data[[membership]],
                                   y = .data[[y_var]])) +
    ggplot2::geom_jitter(ggplot2::aes(color = .data[[membership]]),
                         width = 0.2, size = 2, shape = 1) +
    ggplot2::geom_boxplot(notch = TRUE, fill = NA) +
    ggpubr::theme_classic2() +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme(legend.position = "none",
                   text = ggplot2::element_text(size = 13)) +
    ggplot2::geom_hline(yintercept = mean(df[[y_var]], na.rm = TRUE),
                        linetype = 2) +
    rstatix::stat_compare_means(label = "p.signif",
                                method    = "t.test",
                                ref.group = ref_group) +
    ggplot2::labs(title = title, x = "Cluster", y = y_label)
}


#' Make a symmetric matrix from the lower triangle
#'
#' @param m Square matrix with values in the lower triangle
#' @return Symmetric matrix
make_symm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  m
}


# ==============================================================================
# SECTION 3: Data Loading — Discovery Cohort (Stage C/D, symptomatic HF)
# ==============================================================================
# NOTE: The SQL-loading scripts are institution-specific.
# For reproducibility, save the processed `d`, `Prot_lipids`, and
# `memberships_df` objects to an RData file (see end of Section 4) and
# load them directly if re-running downstream analyses.

message("Loading discovery cohort data ...")
setwd(PATH_DATA_SQL)
source("Read_SQL_MyoVasc_BL_v11.R")
dall1 <- read_SQL_MyoVasc_BL(proteins = TRUE, lipids = TRUE,
                             user = "ese0e", ATC_7_digits = "B01AC06")
setwd(PATH_WORKING)

# Extract omics panels
lipids   <- extract_omics_panel(dall1,
                                var_attr = "lipid_vars_raw",
                                nam_attr = "lipid_nams_raw")
proteins <- extract_omics_panel(dall1,
                                var_attr  = "vars_prot",
                                nam_attr  = "nams_prot",
                                excl_cols = c("NTproBNP (NT-proBNP)", "P16860 (BNP)"))

# Intersect to Stage C/D patients with complete data in both panels
cd_ids  <- dall1$v11_sid01[dall1$hf_stages == "Stage C/D"]
panels_cd <- intersect_panels(list(proteins, lipids), id_subset = cd_ids)
intersection_cd <- rownames(panels_cd[[1]])

# Clinical data frame restricted to the same patients
d <- dall1[dall1$v11_sid01 %in% intersection_cd, ]
message("  Discovery cohort: N = ", nrow(d), " patients")


# ==============================================================================
# SECTION 4: Pre-processing and SNF Clustering — Discovery Cohort
# ==============================================================================

message("Pre-processing omics panels ...")

# Pre-process each panel (imputation → NZV removal → scaling → decorrelation)
prep_results <- lapply(panels_cd, function(panel) {
  imp <- impute_data(panel)
  preprocess_omics(imp, corr_cutoff = CORR_CUTOFF)
})

# Processed (scaled, decorrelated) matrices — used for affinity and RF training
Prot_lipids  <- lapply(prep_results, `[[`, "data_norm")

# Save scaling parameters for projection to validation cohort (Section 7)
train_means  <- lapply(prep_results, `[[`, "train_mean")
train_sds    <- lapply(prep_results, `[[`, "train_sd")

message("Building affinity matrices and fusing networks (SNF) ...")

# Per-panel affinity matrices (used for NMI comparison and C-index evaluation)
W <- build_affinity_matrices(Prot_lipids)

# Fused network (SNF)
W_snf <- SNFtool::SNF(W, K = 20)
# Normalise fused matrix to [0, 1] with diagonal = row max
W_snf <- W_snf - diag(apply(W_snf, 1, max))
W_snf <- W_snf + diag(apply(W_snf, 1, max))
W_snf <- W_snf / max(W_snf)

message("Computing cluster memberships for K = ", K_MIN, " to ", K_MAX, " ...")

# Spectral clustering for individual panels and SNF
mem_proteins_raw  <- compute_memberships(W[[1]])
mem_lipids_raw    <- compute_memberships(W[[2]])
mem_snf_raw       <- compute_memberships(W_snf)

# Rank cluster labels by survival outcome (WHF), ascending risk
memberships_df <- as.data.frame(
  lapply(mem_snf_raw, function(mem) {
    rank_memberships_by_survival(mem, d$whf_all_time, d$whf_all_event)
  })
)
colnames(memberships_df) <- paste0("mem", K_MIN:K_MAX)
memberships_df$freq   <- 1
memberships_df$whf    <- as.factor(d$whf_all_event)
memberships_df$death  <- as.factor(d$tod_event)

# Set final membership used throughout the paper
d$membership <- as.factor(memberships_df[[paste0("mem", K_FINAL)]])

# Save intermediate objects for reuse
save(d, Prot_lipids, memberships_df, W, W_snf,
     train_means, train_sds,
     file = file.path(PATH_WORKING, "processed_discovery_cohort.RData"))
message("  Saved processed_discovery_cohort.RData")


# ==============================================================================
# SECTION 5: Cluster Quality Evaluation (C-index and Log-rank)
# ==============================================================================

message("Evaluating clustering quality (C-index, log-rank) ...")

# Evaluate each modality separately and the SNF fusion
results_prot   <- evaluate_clustering_survival(
  as.data.frame(lapply(mem_proteins_raw, function(m) {
    rank_memberships_by_survival(m, d$whf_all_time, d$whf_all_event)
  })),
  d$whf_all_time, d$whf_all_event
)
results_lipids <- evaluate_clustering_survival(
  as.data.frame(lapply(mem_lipids_raw, function(m) {
    rank_memberships_by_survival(m, d$whf_all_time, d$whf_all_event)
  })),
  d$whf_all_time, d$whf_all_event
)
results_snf    <- evaluate_clustering_survival(memberships_df,
                                               d$whf_all_time, d$whf_all_event)

# Combined C-index comparison plot (Figure 2A in manuscript)
c_index_res <- data.frame(
  k        = K_MIN:K_MAX,
  Lipids   = results_lipids$c_index,
  Proteins = results_prot$c_index,
  SNF      = results_snf$c_index
)
fig_cindex <- plot_c_index(c_index_res)
ggplot2::ggsave(file.path(PATH_SAVE, "fig_cindex_comparison.pdf"),
                fig_cindex, width = 6, height = 4)

# ---- NMI between SNF clusters and each single-panel clustering ---------------
nmi_res <- do.call(rbind, lapply(K_MIN:K_MAX, function(kk) {
  cl_snf  <- SNFtool::spectralClustering(W_snf, kk)
  cl_lip  <- SNFtool::spectralClustering(W[[2]], kk)
  cl_prot <- SNFtool::spectralClustering(W[[1]], kk)
  data.frame(k            = kk,
             NMI_lipid    = aricode::NMI(cl_snf, cl_lip),
             NMI_protein  = aricode::NMI(cl_snf, cl_prot))
}))
print(nmi_res)


# ==============================================================================
# SECTION 6: Survival Analysis — Discovery Cohort
# ==============================================================================

message("Running survival analyses ...")

# ---- 6.1 Kaplan-Meier curves -------------------------------------------------
membership_km <- memberships_df[[paste0("mem", K_FINAL)]]
km_fit <- ggsurvfit::survfit2(
  survival::Surv(whf_all_time, whf_all_event) ~ membership_km,
  data = d
)
fig_km <- km_fit |>
  ggsurvfit::ggsurvfit() +
  ggplot2::labs(x = "Time (years)", y = "Survival probability",
                title = paste0("Worsening of heart failure (K = ", K_FINAL, ")")) +
  ggplot2::scale_color_manual(values = CLUSTER_COLORS) +
  ggsurvfit::add_risktable() +
  ggplot2::ylim(0.5, 1)
ggplot2::ggsave(file.path(PATH_SAVE, "fig_km_whf.pdf"), fig_km, width = 8, height = 6)

# ---- 6.2 Cumulative incidence (competing risks: WHF vs. death) ---------------
# Construct competing-risk status: 1 = WHF within 4 years, 2 = death (competing)
d$membership <- as.factor(memberships_df[[paste0("mem", K_FINAL)]])

event_table <- d |>
  dplyr::mutate(
    cr_status = dplyr::case_when(
      whf_all_event == 1 & whf_all_time <= 4                               ~ 1L,
      tod_event == 1 & tod_time <= 4 &
        (whf_all_event == 0 | tod_time < whf_all_time)                     ~ 2L,
      TRUE                                                                  ~ 0L
    ),
    cr_time = pmin(
      ifelse(whf_all_event == 1, whf_all_time, Inf),
      ifelse(tod_event      == 1, tod_time,     Inf),
      4
    ),
    cluster = memberships_df[[paste0("mem", K_FINAL)]]
  )

cif_result <- cmprsk::cuminc(
  ftime   = event_table$cr_time,
  fstatus = event_table$cr_status,
  group   = event_table$cluster
)
pdf(file.path(PATH_SAVE, "fig_cumulative_incidence.pdf"))
plot(cif_result, lwd = 2)
dev.off()

# ---- 6.3 Cox regression: cluster + covariates --------------------------------
# Residualise NT-proBNP on eGFR to remove renal confounding
d$Nt_proBNP <- residuals(
  stats::lm(nt_pro_bnp_ln ~ egfr, data = d, na.action = na.exclude)
)
d$Cluster <- factor(d$membership,
                    levels = c(as.character(K_FINAL), as.character(seq_len(K_FINAL - 1))))

cox_cluster <- survival::coxph(
  survival::Surv(whf_all_time, whf_all_event) ~ Cluster + age + sex + Nt_proBNP,
  data = d
)
summary(cox_cluster)

# Compare MAGGIC-only vs. MAGGIC + cluster (C-index improvement)
fit_base <- survival::coxph(
  survival::Surv(whf_all_time, whf_all_event) ~ maggic,
  data = d, x = TRUE, y = TRUE
)
fit_ext <- survival::coxph(
  survival::Surv(whf_all_time, whf_all_event) ~ maggic + Cluster,
  data = d, x = TRUE, y = TRUE
)
anova(fit_base, fit_ext, test = "LRT")
AIC(fit_base, fit_ext)

# Harrell's C-index comparison (correlated samples)
lp_base  <- predict(fit_base, type = "lp")
lp_ext   <- predict(fit_ext,  type = "lp")
C_base   <- as.numeric(survival::survConcordance(
  survival::Surv(d$whf_all_time, d$whf_all_event) ~ lp_base)$concordance)
C_ext    <- as.numeric(survival::survConcordance(
  survival::Surv(d$whf_all_time, d$whf_all_event) ~ lp_ext)$concordance)
cc       <- compareC::compareC(d$whf_all_time, d$whf_all_event, lp_base, lp_ext)
message(sprintf("  C-base = %.3f | C-ext = %.3f | ΔC = %.3f | p = %.4f",
                C_base, C_ext, as.numeric(cc$deltaC), as.numeric(cc$pval)))


# ==============================================================================
# SECTION 7: Cluster Characterisation
# ==============================================================================

message("Characterising clusters ...")

# ---- 7.1 HF subtype enrichment (chi-squared + standardised residuals) --------
cluster_vec  <- memberships_df[[paste0("mem", K_FINAL)]]
hf_type_vec  <- droplevels(d$hf_a3)
tab_hf       <- table(cluster_vec, hf_type_vec)
tab_hf_filt  <- tab_hf[, colnames(tab_hf) != "Not classifiable (C/D)"]
chisq_hf     <- stats::chisq.test(tab_hf_filt)
std_res_hf   <- chisq_hf$stdres

# Heatmap of standardised residuals
df_res <- reshape2::melt(std_res_hf,
                         varnames = c("Cluster", "HF_subtype"),
                         value.name = "Std_residual")
fig_stdres <- ggplot2::ggplot(df_res,
                              ggplot2::aes(x = HF_subtype, y = Cluster, fill = Std_residual)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Std_residual)),
                     color = "black", size = 3) +
  ggplot2::scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                                midpoint = 0, limits = c(-4, 4),
                                name = "Std residual") +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 panel.grid = ggplot2::element_blank()) +
  ggplot2::labs(x = "HF subtype", y = "Cluster")
ggplot2::ggsave(file.path(PATH_SAVE, "fig_hf_subtype_residuals.pdf"),
                fig_stdres, width = 6, height = 5)

# ---- 7.2 Omics marker identification (Kruskal-Wallis + cluster z-scores) -----
# Combine both panels into one matrix for marker analysis
prot_lipids_combined <- cbind(Prot_lipids[[1]], Prot_lipids[[2]])
cluster_labels       <- memberships_df[[paste0("mem", K_FINAL)]]

kw_p   <- apply(prot_lipids_combined, 2,
                function(f) stats::kruskal.test(f ~ cluster_labels)$p.value)
kw_fdr <- stats::p.adjust(kw_p, "BH")

# Cluster-wise z-scores relative to global mean/SD
mu_all  <- colMeans(prot_lipids_combined, na.rm = TRUE)
sd_all  <- apply(prot_lipids_combined, 2, sd, na.rm = TRUE)
cl_levs <- levels(as.factor(cluster_labels))

means_per_cluster <- sapply(cl_levs, function(cl) {
  colMeans(prot_lipids_combined[cluster_labels == cl, , drop = FALSE], na.rm = TRUE)
})
Z_markers          <- sweep(means_per_cluster, 1, mu_all, FUN = "-")
Z_markers          <- sweep(Z_markers, 1, sd_all, FUN = "/")
Z_markers[!is.finite(Z_markers)] <- 0
colnames(Z_markers) <- cl_levs

best_cluster <- apply(Z_markers, 1, function(z) colnames(Z_markers)[which.max(abs(z))])
z_best       <- apply(Z_markers, 1, function(z) z[which.max(abs(z))])

# Annotate features as protein or lipid
annotation_type <- rep("Lipid", ncol(prot_lipids_combined))
annotation_type[seq_len(ncol(Prot_lipids[[1]]))] <- "Protein"

marker_tbl <- data.frame(
  feature      = colnames(prot_lipids_combined),
  type         = annotation_type,
  fdr          = kw_fdr,
  best_cluster = best_cluster,
  z_best       = z_best,
  score        = abs(z_best) * (-log10(pmax(kw_fdr, 1e-300))),
  stringsAsFactors = FALSE
)

# ---- 7.3 UpSet plot: cluster-specific features (significant in ≥ K-1 pairs) -
# Features differentially expressed vs. all other clusters
n_clusters <- length(cl_levs)
Proteins_clusters <- t(apply(Prot_lipids[[1]], 2, pairwise_wilcox_sig,
                             cluster_vec = cluster_labels))
Lipids_clusters   <- t(apply(Prot_lipids[[2]], 2, pairwise_wilcox_sig,
                             cluster_vec = cluster_labels))
comb_clusters     <- rbind(Proteins_clusters, Lipids_clusters)

# Extract cluster-exclusive features (significant vs. all others)
cl_specific <- lapply(seq_len(n_clusters), function(i) {
  rownames(comb_clusters)[comb_clusters[, i] == (n_clusters - 1)]
})
names(cl_specific) <- paste0("cl", seq_len(n_clusters))

lt <- cl_specific
m_upset <- ComplexHeatmap::make_comb_mat(lt)
ComplexHeatmap::UpSet(m_upset)

# ---- 7.4 Heatmap of cluster-exclusive features (Figure 3) -------------------
features_ordered    <- unlist(cl_specific)
features_heatmap_sep <- rep(seq_along(cl_specific),
                            lengths(cl_specific))

# Strip UniProt prefix from protein names for display
prot_lipids_combined_named <- prot_lipids_combined
only_prot_idx <- which(features_ordered %in% colnames(Prot_lipids[[1]]))
colnames(prot_lipids_combined_named)[only_prot_idx] <-
  gsub("^.{0,8}(.{1})$", "\\1",
       colnames(prot_lipids_combined_named)[only_prot_idx])

plm     <- prot_lipids_combined_named[, features_ordered]
plmt    <- t(plm)

annotation_features <- ifelse(features_ordered %in% colnames(Prot_lipids[[1]]),
                              "protein", "lipid")
col_fun      <- circlize::colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
omics_colors <- c("lipid" = "#1f78b4", "protein" = "#33a02c")

ra1 <- ComplexHeatmap::rowAnnotation(
  OMICs = annotation_features,
  col   = list(OMICs = omics_colors)
)
hm <- ComplexHeatmap::Heatmap(
  matrix           = plmt,
  name             = "Expression",
  heatmap_height   = grid::unit(40, "cm"),
  show_column_names = FALSE,
  show_row_names   = TRUE,
  row_split        = factor(features_heatmap_sep),
  cluster_row_slices = FALSE,
  cluster_rows     = TRUE,
  cluster_columns  = FALSE,
  column_split     = factor(memberships_df[[paste0("mem", K_FINAL)]]),
  use_raster       = TRUE,
  right_annotation = ra1,
  col              = col_fun
)
Cairo::CairoPNG(file.path(PATH_SAVE, "fig_heatmap_cluster_features.png"),
                width = 3500, height = 6000, res = 300)
ComplexHeatmap::draw(hm, merge_legend = TRUE,
                     padding = grid::unit(c(20, 10, 40, 10), "mm"))
dev.off()

# ---- 7.5 LASSO logistic regression: clinical predictors of cluster membership -
# Define clinical variable sets for LASSO
vars_clinical <- c("age", "sex", "fli_cat", "nic", "hyper", "diab",
                   "dyslip", "adipos", "ckd", "afib", "mi", "stroke",
                   "cad", "pad", "vte")
nams_clinical <- c("Age [y]", "Sex (Women)", "FLI", "Active smoking",
                   "Hypertension", "Diabetes", "Dyslipidemia", "Obesity",
                   "Chronic Kidney Disease", "AF", "MI", "Stroke",
                   "CAD", "PAD", "VTE")

d$fli_cat <- as.factor(ifelse(d$fli >= 60, "yes", "no"))
data_clinic <- d[, colnames(d) %in% vars_clinical]
data_clinic$cluster_membership <- memberships_df[[paste0("mem", K_FINAL)]]
data_clinic <- data_clinic[stats::complete.cases(data_clinic), ]

y_lasso     <- data_clinic$cluster_membership
X_lasso     <- data.matrix(data_clinic[, colnames(data_clinic) != "cluster_membership"])

set.seed(RANDOM_SEED)
cvfit_lasso <- glmnet::cv.glmnet(X_lasso, y_lasso, family = "multinomial",
                                 type.multinomial = "ungrouped",
                                 alpha = 1, nfolds = 10)
fit_lasso   <- glmnet::glmnet(X_lasso, y_lasso, family = "multinomial",
                              type.multinomial = "ungrouped",
                              alpha = 1, lambda = cvfit_lasso$lambda.min)

# Extract and reshape odds ratios
coefs_list  <- coef(fit_lasso)
coef_mat    <- do.call(cbind, lapply(coefs_list, as.matrix))
colnames(coef_mat) <- paste0("cluster", seq_along(coefs_list))
odds_mat    <- exp(coef_mat[-1, ])  # remove intercept row

odds_long <- as.data.frame(odds_mat) |>
  tibble::rownames_to_column("parameter") |>
  tidyr::pivot_longer(-parameter, names_to = "cluster", values_to = "value")

# Save odds ratios for manuscript table
utils::write.csv(odds_long,
                 file.path(PATH_SAVE, "lasso_odds_discovery.csv"),
                 row.names = FALSE)


# ==============================================================================
# SECTION 8: Bootstrap Stability Analysis
# ==============================================================================

message("Running bootstrap stability analysis (B = ", N_BOOTSTRAP, ") ...")
set.seed(RANDOM_SEED)

reference_mem <- as.integer(memberships_df[[paste0("mem", K_FINAL)]])

stability_results <- bootstrap_stability(
  data_list     = Prot_lipids,
  reference_mem = reference_mem,
  B             = N_BOOTSTRAP,
  frac          = BOOTSTRAP_FRAC,
  K_cluster     = K_FINAL
)
message(sprintf("  Median ARI = %.3f (IQR: %.3f–%.3f)",
                median(stability_results$ari),
                quantile(stability_results$ari, 0.25),
                quantile(stability_results$ari, 0.75)))
message(sprintf("  Median NMI = %.3f (IQR: %.3f–%.3f)",
                median(stability_results$nmi),
                quantile(stability_results$nmi, 0.25),
                quantile(stability_results$nmi, 0.75)))

utils::write.csv(stability_results,
                 file.path(PATH_SAVE, "bootstrap_stability.csv"),
                 row.names = FALSE)


# ==============================================================================
# SECTION 9: MOFA2 Comparison
# ==============================================================================

message("Running MOFA2 ...")
set.seed(RANDOM_SEED)

data_list_mofa <- list(
  proteins = t(Prot_lipids[[1]]),
  lipids   = t(Prot_lipids[[2]])
)
MOFAobject   <- MOFA2::create_mofa(data_list_mofa)
model_opts   <- MOFA2::get_default_model_options(MOFAobject)
model_opts$num_factors <- MOFA_FACTORS
train_opts   <- MOFA2::get_default_training_options(MOFAobject)
train_opts$maxiter <- MOFA_MAXITER
train_opts$seed    <- RANDOM_SEED

MOFAobject <- MOFA2::prepare_mofa(MOFAobject,
                                  data_options    = MOFA2::get_default_data_options(MOFAobject),
                                  model_options   = model_opts,
                                  training_options = train_opts)
MOFAmodel  <- MOFA2::run_mofa(MOFAobject,
                              outfile      = file.path(PATH_SAVE, "MOFA_model.hdf5"),
                              use_basilisk = TRUE)

# Extract latent factors (samples × factors)
factors <- MOFA2::get_factors(MOFAmodel, factors = 1:MOFA_FACTORS)[[1]]
Z_mofa  <- as.matrix(factors)

# Spectral clustering on MOFA factors
mofa_clusters <- as.numeric(kernlab::specc(Z_mofa, centers = K_FINAL))

# Compare MOFA vs. SNF cluster assignments (ARI)
snf_clusters <- as.integer(memberships_df[[paste0("mem", K_FINAL)]])
ari_mofa_snf <- mclust::adjustedRandIndex(snf_clusters, mofa_clusters)
message(sprintf("  ARI between SNF and MOFA clusters: %.3f", ari_mofa_snf))


# ==============================================================================
# SECTION 10: Validation Cohort (Stage A/B, asymptomatic)
# ==============================================================================

message("Loading validation cohort (Stage A/B) ...")
setwd(PATH_DATA_SQL)
source("Read_SQL_MyoVasc_BL_v11.R")
dall1_val <- read_SQL_MyoVasc_BL(proteins = TRUE, lipids = TRUE,
                                 user = "ese0e", ATC_7_digits = "B01AC06")
setwd(PATH_WORKING)

lipids_val   <- extract_omics_panel(dall1_val,
                                    var_attr = "lipid_vars_raw",
                                    nam_attr = "lipid_nams_raw")
proteins_val <- extract_omics_panel(dall1_val,
                                    var_attr  = "vars_prot",
                                    nam_attr  = "nams_prot",
                                    excl_cols = c("NTproBNP (NT-proBNP)", "P16860 (BNP)"))

ab_ids      <- dall1_val$v11_sid01[dall1_val$hf_stages %in% c("Stage A", "Stage B")]
panels_ab   <- intersect_panels(list(proteins_val, lipids_val), id_subset = ab_ids)
intersection_ab <- rownames(panels_ab[[1]])
ab <- dall1_val[dall1_val$v11_sid01 %in% intersection_ab, ]
message("  Validation cohort: N = ", nrow(ab), " patients")

# Project validation cohort onto discovery-cohort scaling parameters
Prot_lipids_ab <- lapply(seq_along(panels_ab), function(i) {
  imp       <- impute_data(panels_ab[[i]])
  # Use training means/SDs — never refit scaling on test data
  scaled    <- scale(imp,
                     center = train_means[[i]],
                     scale  = train_sds[[i]])
  scaled
})

# ---- Train one-vs-rest RF classifiers per cluster ----------------------------
message("Training Random Forest classifiers for cluster transfer ...")
# Feature names (cleaned) shared between discovery and validation
colnames_cleaned <- c(
  sub(".*\\(([^)]+)\\).*", "\\1", colnames(Prot_lipids[[1]])),
  sub(".*\\(([^)]+)\\).*", "\\1", colnames(Prot_lipids[[2]]))
)
# Full combined matrix (discovery)
train_matrix_full <- cbind(Prot_lipids[[1]], Prot_lipids[[2]])
colnames(train_matrix_full) <- colnames_cleaned

# Cluster-exclusive feature sets (from Section 7)
features_list_named <- setNames(cl_specific, seq_along(cl_specific))

set.seed(RANDOM_SEED)
rf_models <- lapply(seq_along(features_list_named), function(i) {
  feats <- features_list_named[[i]]
  mat_i <- train_matrix_full[, colnames(train_matrix_full) %in% feats, drop = FALSE]
  train_rf_one_vs_rest(mat_i,
                       cluster_label = as.integer(memberships_df[[paste0("mem", K_FINAL)]]),
                       target_k      = i)
})
names(rf_models) <- as.character(seq_along(rf_models))

# Predict cluster membership in the validation cohort
val_matrix_full <- cbind(Prot_lipids_ab[[1]], Prot_lipids_ab[[2]])
colnames(val_matrix_full) <- colnames_cleaned

cluster_vector_ab <- assign_clusters_rf(rf_models, val_matrix_full,
                                        prob_cutoff = RF_PROB_CUTOFF)

# Assign labels; samples with no high-probability assignment get label = K+1
ab$cluster        <- as.factor(cluster_vector_ab)
ab <- ab[ab$cluster != (K_FINAL + 1), ]   # remove "unassigned" samples
ab$cluster        <- droplevels(ab$cluster)
ab$binary_cluster <- ifelse(ab$cluster %in% c(1, 2, 3, 4), 1, 0)
message("  Validation cohort after cluster assignment: N = ", nrow(ab))

# ---- Survival analysis in validation cohort ----------------------------------
ab$whf_all_event <- as.numeric(as.character(ab$whf_all_event))
ab$whf_all_time  <- as.numeric(ab$whf_all_time)

cox_val <- survival::coxph(
  survival::Surv(whf_all_time, whf_all_event) ~ binary_cluster + age + sex + nt_pro_bnp,
  data = ab, x = TRUE, y = TRUE, model = TRUE
)
summary(cox_val)

# Forest plot of adjusted HRs (validation cohort)
tidy_cox_val <- broom::tidy(cox_val, exponentiate = TRUE, conf.int = TRUE)
label_map_val <- c(binary_cluster = "High-risk clusters",
                   age            = "Age (per year)",
                   sexWomen       = "Women (vs Men)",
                   nt_pro_bnp     = "NT-proBNP")
forest_df_val <- tidy_cox_val |>
  dplyr::mutate(var = dplyr::recode(term, !!!label_map_val),
                var = factor(var, levels = unname(label_map_val))) |>
  dplyr::arrange(var)

fig_forest_val <- ggplot2::ggplot(
  forest_df_val,
  ggplot2::aes(y = var, x = estimate, xmin = conf.low, xmax = conf.high)
) +
  ggplot2::geom_vline(xintercept = 1, linetype = 2) +
  ggplot2::geom_errorbarh(height = 0.18) +
  ggplot2::geom_point(size = 2.6) +
  ggplot2::scale_x_log10(name = "Hazard Ratio (log scale)") +
  ggplot2::ylab(NULL) +
  ggplot2::ggtitle("Adjusted Hazard Ratios — Validation Cohort") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
ggplot2::ggsave(file.path(PATH_SAVE, "fig_forest_validation.pdf"),
                fig_forest_val, width = 5, height = 4)

# Stage progression at follow-up
setwd(PATH_DATA_SQL)
source("Read_SQL_MyoVasc_FU2_v5.R")
dall_fu2 <- read_SQL_MyoVasc_FU2(proteins = FALSE, user = "ese0e")
setwd(PATH_WORKING)

ids_cluster   <- data.frame(id = ab$v11_sid01,
                            bl1 = ab$hf_stages,
                            binary_cluster = ab$binary_cluster)
fu2_stages    <- data.frame(id = dall_fu2$v11_sid01_bl,
                            fu = dall_fu2$hf_stages)
stage_change  <- merge(ids_cluster, fu2_stages, by = "id")
stage_change  <- stage_change[stage_change$fu != "Stage 0", ]

# Chi-squared test: binary cluster vs. follow-up stage
tab_stage <- table(stage_change$binary_cluster, stage_change$fu)
print(stats::chisq.test(tab_stage))

fig_stage <- ggplot2::ggplot(stage_change,
                             ggplot2::aes(x = factor(binary_cluster), fill = fu)) +
  ggplot2::geom_bar(position = "fill") +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::labs(x = "Binary cluster", y = "Proportion at follow-up",
                fill = "Follow-up HF stage") +
  ggplot2::theme_minimal()
ggplot2::ggsave(file.path(PATH_SAVE, "fig_stage_progression.pdf"),
                fig_stage, width = 5, height = 4)


# ==============================================================================
# SECTION 11: Clinical Boxplots for Selected Cluster Phenotypes
# ==============================================================================

message("Generating cluster-phenotype boxplots ...")

# Extract specific proteins from the combined matrix by column index
# (indices correspond to the combined protein+lipid matrix order)
prot_lipids_myo <- cbind(Prot_lipids[[1]], Prot_lipids[[2]])
d$membership    <- as.factor(memberships_df[[paste0("mem", K_FINAL)]])

# ---- Cluster 1: Kidney / liver / cardiac phenotype --------------------------
# (Add or adjust variable names to match your cohort column names)
fig_cl1 <- ggpubr::ggarrange(
  plot_cluster_boxplot(d, "crp_ln",     title = "CRP (log)",     y_label = "ln(CRP)"),
  plot_cluster_boxplot(d, "egfr",       title = "eGFR",          y_label = "eGFR (ml/min/1.73m²)"),
  plot_cluster_boxplot(d, "nt_pro_bnp_ln", title = "NT-proBNP",  y_label = "ln(NT-proBNP)"),
  ncol = 3, nrow = 1
)
ggplot2::ggsave(file.path(PATH_SAVE, "fig_cluster1_phenotype.pdf"),
                fig_cl1, width = 12, height = 5)

# ---- Cluster 4: Metabolic phenotype ----------------------------------------
d$tri_ln <- log10(d$v34_lbl55)
fig_cl4 <- ggpubr::ggarrange(
  plot_cluster_boxplot(d, "hba1c",    title = "HbA1c",        y_label = "HbA1c [%]",
                       colors = c("1" = "#a9a9a9","2" = "#a9a9a9","3" = "#a9a9a9",
                                  "4" = "#FF9900","5" = "#a9a9a9","6" = "#a9a9a9",
                                  "7" = "#a9a9a9","8" = "#a9a9a9")),
  plot_cluster_boxplot(d, "homa_ln",  title = "HOMA-IR",       y_label = "ln(HOMA-IR)",
                       colors = c("1" = "#a9a9a9","2" = "#a9a9a9","3" = "#a9a9a9",
                                  "4" = "#FF9900","5" = "#a9a9a9","6" = "#a9a9a9",
                                  "7" = "#a9a9a9","8" = "#a9a9a9")),
  plot_cluster_boxplot(d, "tri_ln",   title = "Triglycerides",  y_label = "ln(Triglycerides)",
                       colors = c("1" = "#a9a9a9","2" = "#a9a9a9","3" = "#a9a9a9",
                                  "4" = "#FF9900","5" = "#a9a9a9","6" = "#a9a9a9",
                                  "7" = "#a9a9a9","8" = "#a9a9a9")),
  ncol = 3, nrow = 1
)
ggplot2::ggsave(file.path(PATH_SAVE, "fig_cluster4_metabolic.pdf"),
                fig_cl4, width = 12, height = 5)


# ==============================================================================
# SECTION 12: Session Information (for reproducibility)
# ==============================================================================

message("\nSession information saved to session_info.txt")
sink(file.path(PATH_SAVE, "session_info.txt"))
utils::sessionInfo()
sink()