# Leakage Audit

각 outer fold 실행 시 아래 항목을 자동 또는 수동 증거와 함께 채운다.

| 항목 | 상태 | 증거 |
|---|---|---|
| train/test cell-line ID disjoint | PASS | `test_splits.py`; `split_audit.json` |
| alias/pool을 통한 동일 line 중복 없음 | PASS | 94 unique cell line, 94 unique DepMap ID in `cell_line_annotations.csv` |
| test treated expression이 feature fit에 미사용 | PASS | `fit_predict_baselines` API는 outer-test response를 인자로 받지 않음 |
| gene filtering은 outer train에서만 fit | PASS | `fit_control_embedding`을 outer/inner training control마다 재호출 |
| control PCA는 outer train에서만 fit | PASS | outer/inner training control fit 후 test/validation transform |
| response basis는 outer train에서만 fit | N/A | B0–B4는 response PCA basis를 사용하지 않음; CCLR 단계에서 재감사 |
| normalization의 학습 통계가 test에 미의존 | PASS | 각 pseudobulk row별 고정 `log1p(CPM)`; cohort 통계 없음 |
| sensitivity는 주 모델 입력에 미사용 | PASS | `models.yaml=false`; `baseline_summary.json=false` |
| mutation은 주 모델 입력에 미사용 | PASS | `models.yaml=false`; `baseline_summary.json=false` |
| hyperparameter는 inner CV에서만 선택 | PASS | `baseline_inner_cv.csv`, fold별 `baseline_hyperparameters.csv` |
| test 결과를 보고 split/metric 미변경 | PASS | run log에 split/model config SHA-256 저장; deterministic prediction hash 확인 |

Baseline prediction Parquet은 동일 seed/config 재실행에서 SHA-256
`7d905af11b10f178954a2a9b5c0518a80efe2a2212054c4dbb23160356ad797c`로 일치했다.
