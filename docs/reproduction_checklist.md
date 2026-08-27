# 최종 재현 체크리스트

작성일: 2026-08-27

Release: 1.0.0

## 완료 게이트

- [x] 공식 experiment 3 Matrix Market의 32,738 genes, non-negative integer count, ordered gene 일치 검증
- [x] strict 94-line time-matched와 pooled-control 97-line 코호트 재현
- [x] 희소 pseudobulk, 24h response, processed manifest와 checksum 생성
- [x] DMSO 6h–24h source/time 차이 및 사전 marker 방향 검증
- [x] 94-line annotation 100% join과 lineage-aware outer/inner split 동결
- [x] B0–B4 470 held-out predictions, paired bootstrap, deterministic rerun
- [x] CCLR 94 held-out predictions, fold artifacts, reconstruction, 독립 validator
- [x] W7 16-model ablation, leakage audit, subspace/pathway 진단, 독립 validator
- [x] W8 24-line/13,713-cell time-course, 22-line 기술 분석, 17-line external transfer
- [x] W9 38개 biological/error 검정, BH-FDR 재계산, partial association
- [x] W10 threshold/gene-filter/bootstrap/subsampling/LOLO/noise-floor 강건성 분석
- [x] W11 사전 gate 통과 후 17-line/1,892-cell latent mean-shift 분포 분석
- [x] W12 564 final predictions, 64 metrics, 94 line table, 10 figures, 6 tables 동결
- [x] W13 최종 보고서, supplementary, limitations 작성
- [x] W14 전체 unit test·stage validator·notebook 재실행
- [x] W15 독립 release validator와 Git 상태/commit 기록

## 핵심 재현 결과

- Core: 94 lines, 16,588 normal cells, 32,738 genes
- DMSO 6h–24h median PCC: `0.9727`
- Control source/time RMSE: `0.4054`; trametinib response RMSE: `0.4536`; ratio `0.8879`
- B1 RMSE: `0.323750`
- B4 RMSE: `0.322035`; PCC-context: `0.107666`
- B4 gain vs B1: `0.001715` (95% CI `0.001168–0.002284`, about `0.53%`)
- CCLR RMSE: `0.322383`; gain vs B4: `-0.000348` (95% CI `-0.000584–-0.000105`)
- External B4 gain: 3h `-0.007310`, 6h `-0.004712`, 12h `-0.001242`, 24h `0.002074`, 48h `0.003839`
- Equal-20-cell B4 gain: negative in `5/5` repeats
- LOLO B4 gain: `0.000514` (95% CI `-0.000132–0.001178`)
- Split-half full-target floor: `0.279028`; PCC `0.236373`
- W11 B4 Energy gain: `0.570597`; sliced-Wasserstein gain: `0.163040`

## 명령

```bash
uv sync --all-groups
uv run ruff check src tests scripts
uv run pytest -q
uv run python -m yakseopdong smoke
uv run python -m yakseopdong validate-cclr
uv run python -m yakseopdong validate-ablation
uv run python -m yakseopdong validate-temporal
uv run python -m yakseopdong validate-biology
uv run python -m yakseopdong validate-robustness
uv run python -m yakseopdong validate-distribution
uv run python -m yakseopdong freeze-release
uv run python -m yakseopdong validate-release
uv run python scripts/build_notebooks.py --execute
git diff --check
```

## Frozen artifact hashes

- `final_predictions.parquet`: `5ffe460c0a3019fa6369d483fa6dfd2ec771494e472418ab401303c555ee6562`
- `final_metrics.csv`: `bd1c448709a69b3858f0d91451615e8ed129f26450f574917a7df69297ab281f`
- `final_cell_lines.csv`: `4fd26f2484a62027ac184759b1b71ecb9e23d91a0716aa6d3e73891af78027ed`

Figure/table manifest hashes are recorded in `release/final_config.json` because the study-design figure and corrected cohort counts are part of the final release freeze.

## 해석 제한

높은 control PCC만으로 pooling을 정당화할 수 없다. DMSO 6h와 24h는 별도 source이며 time-course도 core와 batch/source가 다르다. Pseudobulk와 W11은 동일 세포의 전후 trajectory가 아니다. B4 gain은 작고 equal-cell subsampling에 민감하므로 strong personalization claim은 금지한다. 자세한 경계는 `report/limitations.md`에 있다.
