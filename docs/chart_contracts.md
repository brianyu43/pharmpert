# Chart contracts

## W8 temporal components and pathways

- Question: which shared transcriptional programs change across five measured times?
- Grain: cell-line pseudobulk response; 22 eligible lines per time.
- Encodings: x = discrete hours; y = mean PCA or pathway score; color = component/pathway;
  interval = 1.96 SE for descriptive PCA and cell-line bootstrap 95% CI for pathways.
- Comparison logic: same line set and same time-matched Tram−DMSO definition at every time.
- Caveats: PCA sign is arbitrary; five discrete observations are not a continuous trajectory;
  pathway means do not test enrichment.

## W8 heterogeneity and frozen-model transfer

- Question: does cross-line dispersion change, and when does a frozen 24h model transfer?
- Grain: left, 22-line descriptive cohort; right, 17 external DepMap IDs only.
- Encodings: x = discrete hours; y = dispersion or macro RMSE; color/shape = metric/model;
  interval = cell-line bootstrap 95% CI where applicable.
- Comparison logic: all transfer models were fit once on the strict 94-line 24h cohort;
  no time-course response was used for fitting or model selection.
- Caveats: source/batch and cohort differ from the 94-line training experiment, and small
  RMSE gains do not establish strong individualized prediction.

## W9 biological correlates

- Question: which observed 24h response programs covary with sensitivity and MAPK genotype?
- Grain: 94 cell lines; every point is one cell-line pseudobulk response.
- Encodings: scatter x = external author sensitivity, y = descriptive response PC1,
  color = BRAF/KRAS status; interval plot y = predefined pathway, x = Spearman rho.
- Comparison logic: all association variables are interpretation-only and BH-FDR corrected.
- Caveats: PCA sign is arbitrary, lineage/genotype are observational, and no association
  variable entered B4 or CCLR.

## W9 held-out error diagnostics

- Question: is absolute error associated with target scale or baseline-state novelty?
- Grain: one frozen B4 outer-fold prediction per cell line.
- Encodings: x = observed response RMS or nearest outer-training PC20 distance;
  y = held-out RMSE; annotations = three largest errors.
- Comparison logic: baseline PCA is fit independently inside each outer training fold.
- Caveats: response RMS and RMSE share scale; partial-rank analyses separately evaluate
  cell-count support, and all driver results remain observational.
