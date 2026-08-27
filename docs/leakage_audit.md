# Leakage Audit

각 outer fold 실행 시 아래 항목을 자동 또는 수동 증거와 함께 채운다.

| 항목 | 상태 | 증거 |
|---|---|---|
| train/test cell-line ID disjoint | PASS | `test_splits.py`; `split_audit.json` |
| alias/pool을 통한 동일 line 중복 없음 | PASS | 94 unique cell line, 94 unique DepMap ID in `cell_line_annotations.csv` |
| test treated expression이 feature fit에 미사용 | PASS | `fit_predict_baselines`, `fit_predict_cclr`, `predict_fixed_low_rank` API는 outer-test response를 인자로 받지 않음 |
| gene filtering은 outer train에서만 fit | PASS | `fit_control_embedding`을 outer/inner training control마다 재호출 |
| control PCA는 outer train에서만 fit | PASS | outer/inner training control fit 후 test/validation transform |
| response basis는 outer train에서만 fit | PASS | CCLR inner/outer training response마다 `fit_response_embedding` 재호출; test response는 예측 후 score 평가에만 사용 |
| W7 pathway panel이 response 결과에 미의존 | PASS | MSigDB Hallmark 2026.1.Hs 6개 세트+수동 8 markers를 실행 전 `genesets.yaml`에 동결; dataset 교집합 후 outer-training control 평균·분산만 사용 |
| W7 lineage/BRAF/KRAS encoding이 outer train에서만 fit | PASS | lineage category 및 BRAF/KRAS 평균·표준편차를 각 outer-training fold에서 fit; unseen lineage는 all-zero |
| normalization의 학습 통계가 test에 미의존 | PASS | 각 pseudobulk row별 고정 `log1p(CPM)`; cohort 통계 없음 |
| sensitivity는 주 모델 입력에 미사용 | PASS | `models.yaml=false`; baseline/CCLR summary=false |
| mutation은 주 모델 입력에 미사용 | PASS | baseline/CCLR은 `models.yaml=false`; W7의 명시적 BRAF/KRAS ablation에만 training-fitted 입력으로 사용 |
| hyperparameter는 inner CV 또는 사전 고정 | PASS | baseline/CCLR `inner_cv.csv`; W7은 결과 전 v1.2에서 d/rank/alpha와 전체 변형을 고정 |
| outer test로 W7 변형을 선택하지 않음 | PASS | 16개 고정 비교군 전체 보고; `selection_rule=report_all_fixed_variants_without_outer_test_selection` |
| test 결과를 보고 split/metric 미변경 | PASS | run log에 split/model config SHA-256 저장; deterministic prediction hash 확인 |

Baseline prediction Parquet은 동일 seed/config 재실행에서 SHA-256
`7d905af11b10f178954a2a9b5c0518a80efe2a2212054c4dbb23160356ad797c`로 일치했다.

CCLR prediction Parquet도 동일 seed/config 재실행에서 SHA-256
`68387db71ec768bde25e6934cd03193d131f930ca9c3c509f61820f6af2fd31b`로 일치했다.
fold artifact에서 직접 예측을 재구성한 최대 절대 오차는 `2.3e-7`이었다.

W7은 16개 변형 × 94 lines의 1,504 metric row를 모두 보존한다. 동일 코드·seed로
`ablation_metrics.csv` SHA-256이 재실행에서 일치하며, paired gain과 2,000회 bootstrap
CI는 `validate-ablation`에서 독립 재계산한다. W7에서는 outer-test 결과로 변형을
선택하지 않았고 sensitivity도 입력하지 않았다.
