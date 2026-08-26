# Leakage Audit

각 outer fold 실행 시 아래 항목을 자동 또는 수동 증거와 함께 채운다.

| 항목 | 상태 | 증거 |
|---|---|---|
| train/test cell-line ID disjoint | PENDING | split test |
| alias/pool을 통한 동일 line 중복 없음 | PENDING | metadata audit |
| test treated expression이 feature fit에 미사용 | PENDING | pipeline test |
| gene filtering은 outer train에서만 fit | PENDING | pipeline test |
| control PCA는 outer train에서만 fit | PENDING | pipeline test |
| response basis는 outer train에서만 fit | PENDING | pipeline test |
| normalization의 학습 통계가 test에 미의존 | PENDING | pipeline test |
| sensitivity는 주 모델 입력에 미사용 | PENDING | config audit |
| hyperparameter는 inner CV에서만 선택 | PENDING | run log |
| test 결과를 보고 split/metric 미변경 | PENDING | Git history/protocol hash |
