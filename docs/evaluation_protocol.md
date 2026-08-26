# Evaluation Protocol v0.1 — Pre-data freeze

이 문서는 모델 결과를 보기 전에 평가 원칙을 고정한다. 실제 split ID, gene universe, pathway gene set 버전은 Stage 1–4에서 원자료 감사 후 추가하고 Git commit으로 동결한다.

## 일반화 단위

- 일반화 단위는 cell이 아니라 `cell_line`이다.
- 같은 cell line의 cell, pseudobulk, pool alias가 outer train/test에 동시에 존재해서는 안 된다.
- outer test line의 treated expression은 예측 생성이 끝날 때까지 어떤 전처리 학습에도 사용하지 않는다.

## 분할

- outer: cell line 기준 deterministic 5-fold
- inner: outer-training cell line 안에서만 hyperparameter 선택
- lineage 빈도를 감사한 뒤 가능한 범위에서 fold 균형화
- 희소 lineage를 합치거나 별도 처리할 규칙은 성능 확인 전에 기록
- 모든 PCA, SVD, gene filtering, scaling, response basis는 outer-training fold에서 fit

시간축 데이터는 24시간 core 모델의 hyperparameter 선택에 사용하지 않는다.

## 표기

Outer-training fold의 평균 반응을 `mu_train`이라 한다.

```text
residual_observed = delta_observed - mu_train
residual_predicted = delta_predicted - mu_train
```

## 전체 반응 지표

- `RMSE-delta = sqrt(mean((predicted - observed)^2))`
- `NRMSE-delta = RMSE-delta / max(RMS(observed), 1e-8)`
- `PCC-delta`: 유전자별 Pearson correlation
- `Spearman-delta`: 유전자별 순위 상관
- `Pathway score MAE/correlation`
- `Signed top-k overlap`: training-fold gene universe에서 상향 50개와 하향 50개를 각각 비교

관측 또는 예측 벡터의 분산이 0이면 Pearson/Spearman은 `NA`다. 특히 B0의 correlation을 0으로 대체하지 않는다.

## 맥락 추가 정보 지표

- `PCC-context`: `residual_predicted`와 `residual_observed`의 Pearson correlation
- `RMSE gain vs B1 = RMSE(B1) - RMSE(model)`
- `NRMSE gain vs B1 = NRMSE(B1) - NRMSE(model)`

양의 gain은 해당 모델이 공통 평균 반응 B1보다 낫다는 뜻이다. 프로젝트의 맥락 예측 주장은 raw PCC만이 아니라 gain과 residual metric을 함께 근거로 한다.

## 집계와 불확실성

- 세포 수 가중치 없는 cell-line macro-average가 주 집계다.
- 모델 비교는 같은 held-out cell line의 paired difference로 계산한다.
- 기본 CI는 cell line을 2,000회 paired bootstrap한 percentile 95% CI다.
- 기본 bootstrap seed는 `20260827`이다.
- fold별 값과 세포주별 값을 모두 보존한다.

## 사전 판정 규칙

- B1 대비 개선 CI가 0을 포함하면 일관된 개선이라고 주장하지 않는다.
- 일부 lineage에서만 나타나는 개선은 전체 일반화 성능으로 표현하지 않는다.
- split, primary metric, top-k, bootstrap seed는 test 결과를 본 뒤 바꾸지 않는다.
- 규칙 변경이 불가피하면 이전 결과와 분리된 protocol version 및 이유를 남긴다.
