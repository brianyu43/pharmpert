# Supplementary Material

## S1. Frozen cohort contracts

| Cohort | Lines | Cells | Role |
|---|---:|---:|---|
| core 24h strict | 94 | 16,588 | primary held-out benchmark |
| pooled-control sensitivity | 97 | 25,130 | control-source sensitivity |
| time-course all | 24 | 13,713 | coverage audit |
| time-course threshold-10 | 22 | 12,078 | temporal description, control/trametinib only |
| temporal external | 17 | 1,892 | transfer and latent-distribution extension |

Core eligibility requires at least 20 `cell_quality == normal` cells in both `DMSO_24hr_expt3` and `Trametinib_24hr_expt3`. Pooled sensitivity replaces the control with `DMSO_6hr + DMSO_24hr`. These cohorts are never silently mixed.

## S2. Split and leakage contract

- Seed: `20260827`
- Generalization unit: cell line
- Outer split: lineage-aware deterministic 5-fold, fold sizes 18–19
- Inner split: each outer-training set only, lineage-aware 4-fold
- Test treated response is never used for feature filtering, PCA, basis fitting, normalization fitting, or hyperparameter selection
- Sensitivity and mutation are interpretation-only
- Every prediction row records `treated_response_used_for_fit=False`
- W8 external transfer excludes the five core/time-course overlapping lines

The executable audit is in `docs/leakage_audit.md`; stage validators are in `results/logs/*_validation.json`.

## S3. Model definitions

| Model | Predictor | Target/operation | Selection |
|---|---|---|---|
| B0 | none | zero response | none |
| B1 | none | outer-training mean response | none |
| B2 | lineage | training lineage mean | B1 fallback |
| B3 | control PCA | nearest training line response | inner-CV PCA dimension |
| B4 | control PCA | full response multi-output ridge | nested-CV dimension and alpha |
| CCLR | control PCA | response-PC ridge and reconstruction | nested-CV dimension, rank, alpha |

B4 is the final primary model because it achieved the lowest prespecified held-out macro RMSE. W7 fixed variants are diagnostic and were not used to replace it after observing outer-test results.

## S4. W7 diagnostic matrix

Sixteen comparison rows include B1, B4, nested CCLR and 13 fixed CCLR variants. Fixed variants start from `control dimension=20`, `response rank=20`, `alpha=100` and change one factor:

- control dimension: 5, 10, 20, 30
- response rank: 2, 5, 10, 20, 30, 40, 50
- pathway input panel
- lineage one-hot appended to control PCs
- BRAF/KRAS binary features appended to control PCs

Selected results:

| Variant | RMSE | Gain vs B1 | Gain vs B4, 95% CI |
|---|---:|---:|---:|
| fixed full d20/r20 | 0.322350 | 0.001399 | -0.000315 [-0.000612, -0.000044] |
| pathway d20/r20 | 0.322286 | 0.001464 | -0.000251 [-0.000598, 0.000127] |
| lineage d20/r20 | 0.322329 | 0.001420 | -0.000295 [-0.000585, 0.000025] |
| BRAF/KRAS d20/r20 | 0.322162 | 0.001588 | -0.000127 [-0.000624, 0.000383] |

## S5. Temporal details

Immediate-early mean response was negative at every time point: 3h `-1.413`, 6h `-1.514`, 12h `-1.597`, 24h `-1.556`, 48h `-1.626`. E2F changed from `0.030` at 3h to `-0.544` at 48h; G2M changed from `0.004` to `-0.469`.

External B4 transfer gain relative to B1:

| Time | Gain | 95% CI | Interpretation |
|---:|---:|---:|---|
| 3h | -0.007310 | [-0.010007, -0.004820] | worse |
| 6h | -0.004712 | [-0.007077, -0.002399] | worse |
| 12h | -0.001242 | [-0.003043, 0.000515] | inconclusive |
| 24h | 0.002074 | [0.000311, 0.004036] | small improvement |
| 48h | 0.003839 | [0.002176, 0.005497] | small improvement |

## S6. Biological validation

Thirty-eight tests were corrected together with Benjamini–Hochberg; 22 were below FDR 0.05. All are interpretation-only. Key effect sizes:

- response PC1 vs sensitivity: Spearman `-0.668273`
- E2F vs sensitivity: `-0.617108`
- G2M vs sensitivity: `-0.601835`
- immediate early vs sensitivity: `-0.424831`
- B4 RMSE vs response RMS: `0.972806`
- B4 RMSE vs min-condition cells: `-0.774982`
- B4 RMSE vs baseline novelty: `0.033371`
- partial B4 RMSE vs response RMS controlling cells: `0.953481`
- partial B4 RMSE vs cells controlling response RMS: `-0.568730`

## S7. Robustness matrix

| Test | Result | Status |
|---|---|---|
| bootstrap seeds 20260827–29 | gain 0.001715; all CI > 0 | robust within frozen split |
| minimum cells 10/20/30 | gain 0.001470–0.001892; all CI > 0 | robust |
| top variable genes 1k/3k/5k/10k | gain 0.001246–0.001637; all CI > 0 | robust |
| remove sensitivity tails | gain 0.001314 [0.000843, 0.001780] | robust |
| equal 20 cells, five repeats | gain -0.001443 to -0.001677; 0/5 positive | sensitive |
| leave-one-lineage-out | gain 0.000514 [-0.000132, 0.001178] | inconclusive |
| split-half target | floor 0.279028; PCC 0.236373 | measurement limitation |

## S8. Single-cell distribution gate

The extension uses training-derived gene selection and PCA (`2,000` genes, `10` PCs) and 17 external lines. Each observed control cell receives the same line-level predicted shift. B4 improved Energy distance by `0.570597` and sliced-Wasserstein by `0.163040` relative to B1, with both paired CIs above zero.

This construction preserves within-line control covariance exactly. It is a test of line-level latent mean translation, not a generative single-cell model. Distribution metrics were not used to tune or select B4/CCLR.

## S9. Frozen deliverables

| Artifact | Rows/items | SHA-256 |
|---|---:|---|
| `results/final_predictions.parquet` | 564 | `5ffe460c0a3019fa6369d483fa6dfd2ec771494e472418ab401303c555ee6562` |
| `results/final_metrics.csv` | 64 | `bd1c448709a69b3858f0d91451615e8ed129f26450f574917a7df69297ab281f` |
| `results/final_cell_lines.csv` | 94 | `4fd26f2484a62027ac184759b1b71ecb9e23d91a0716aa6d3e73891af78027ed` |
| `results/figure_manifest.csv` | 10 | recorded in `release/final_config.json` |
| `results/table_manifest.csv` | 6 | recorded in `release/final_config.json` |

The prediction Parquet contains a fixed-size 32,738-vector in `predicted_delta_log1p_cpm`; model rows are B0–B4 and CCLR, 94 each.

## S10. Reproduction commands

```bash
uv sync --all-groups
uv run python -m yakseopdong smoke
uv run pytest -q
uv run python -m yakseopdong validate-cclr
uv run python -m yakseopdong validate-ablation
uv run python -m yakseopdong validate-temporal
uv run python -m yakseopdong validate-biology
uv run python -m yakseopdong validate-robustness
uv run python -m yakseopdong validate-distribution
uv run python -m yakseopdong freeze-release
uv run python -m yakseopdong validate-release
uv run python scripts/build_notebooks.py --execute
```

`data/raw/` and `data/processed/` are intentionally untracked. Source and processed manifests contain paths, sizes, shapes and SHA-256 hashes. A clean checkout therefore needs the declared source files before full analysis regeneration; validators can then confirm exact contracts.

## S11. Artifact index

- Main figures F1–F8: `results/figure_manifest.csv`
- Supplementary figures S1–S2: `results/figure_manifest.csv`
- Main tables T1–T5 and cell-line table S1: `results/table_manifest.csv`
- Notebook trail: `notebooks/00_data_audit.ipynb` through `notebooks/09_single_cell_extension.ipynb`
- Final configuration: `release/final_config.json`
- Final validation: `results/logs/release_validation.json`
