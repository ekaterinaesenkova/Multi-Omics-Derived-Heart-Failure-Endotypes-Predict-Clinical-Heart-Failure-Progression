
# Multi-Omics-Derived Heart Failure Endotypes Predict Clinical Heart Failure Progression

Proteomics + lipidomics integration via Similarity Network Fusion (SNF) to stratify
Stage C/D heart failure patients and validate cluster-based risk in an independent
asymptomatic cohort.

> **Data access:** Raw patient data cannot be shared due to privacy restrictions.
> Aggregated intermediate objects and summary statistics are provided in `/data`.
> Contact the corresponding author for data access requests.

---

## Cohort overview

| | |
|---|---|
| **Discovery cohort** | Stage C/D HF |
| **Validation cohort** | Stage A/B HF |
| **Omics layers** | Proteomics (Olink) + lipidomics |
| **Primary outcome** | Worsening heart failure (WHF), time-to-event |

---

## Analysis pipeline

1. **Data preprocessing** — Imputation, near-zero-variance removal, z-score scaling,
   and decorrelation (Pearson r > 0.9) applied independently to each omics panel.

2. **Network fusion (SNF)** — Per-panel affinity matrices fused via Similarity Network
   Fusion; spectral clustering tested over K = 3–10.

3. **Cluster quality evaluation** — Harrell's C-index, log-rank tests, and NMI used
   to select the optimal K and compare SNF against single-modality clusterings.

4. **Cluster characterisation** — Kruskal-Wallis + BH FDR on all features;
   cluster-wise z-scores; UpSet plots for exclusive markers; LASSO multinomial
   regression for clinical predictors.

5. **Stability & comparison** — Bootstrap subsampling (B = 100, 80% of samples)
   reporting ARI and NMI; MOFA2 run as a comparison multi-omics method.

6. **Validation cohort transfer** — Per-cluster Random Forest classifiers trained on
   the discovery cohort project cluster labels onto Stage A/B patients; survival and
   stage-progression outcomes tested in the transferred clusters.

---

## Repository structure
