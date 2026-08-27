# Evaluation Protocol v1.3 — Frozen baseline through distribution extension

이 문서는 모델 결과를 보기 전에 정한 평가 원칙과 Stage 4에서 고정한 실제 split/feature 규칙을 기록한다. Baseline과 W6 CCLR의 주 결과는 v1.1, W7 ablation은 v1.2 규칙을 유지한다. W8–W11 확장은 주 모델을 소급 선택하지 않는 v1.3 해석·강건성 프로토콜이다.

## Frozen implementation

- seed: `20260827`
- outer: `split_assignments.csv`, lineage-aware 5-fold, fold size 18–19
- inner: `inner_split_assignments.csv`, 각 outer-train 안의 lineage-aware 4-fold, fold size 18–19
- lineage: Figshare v3 `all_CL_features.rds`의 저자 `Disease` 필드, 원래 21개 label 유지
- target: strict 94-line의 전체 32,738-gene `Trametinib 24h − DMSO 24h`
- B3/B4 control feature: 각 training partition에서 평균 `log1p(CPM) ≥ 0.1`, nonzero variance인 gene 중 분산 상위 5,000개
- B3/B4 PCA dimensions: `[5, 10, 20, 30]`; whitening 사용
- B4 ridge alpha: `[0.01, 0.1, 1, 10, 100]`
- B2: 같은 lineage training line이 2개 미만이면 B1 fallback
- sensitivity/mutation: predictor 사용 금지, 해석 전용

### W6 CCLR extension frozen before execution

- control feature filter/PCA: B4와 동일한 training-only 규칙
- control PCA dimensions: `[5, 10, 20, 30]`; whitening 사용
- response PCA ranks: `[2, 5, 10, 20]`; whitening 미사용
- CCLR ridge alpha: `[0.01, 0.1, 1, 10, 100]`
- selection metric: inner-validation cell-line macro RMSE-Δ
- exact-tie rule: 작은 response rank, 작은 control dimension, 큰 alpha 순
- response target: 전체 32,738 genes
- 비교: 같은 held-out line에서 B1 및 B4와 paired RMSE/NRMSE difference
- 결과를 본 뒤 W6 grid를 확장하지 않으며, 경계값 검사는 W7 ablation으로 분리

### W7 diagnostic ablation frozen before execution

- outer split, target, bootstrap 단위와 seed는 W5/W6와 동일하다.
- 모든 W7 변형은 `control dimension=20`, `response rank=20`, `ridge alpha=100`을 기준점으로 하며 한 요인만 바꾼다.
- control dimension은 `[5, 10, 20, 30]`, response rank는 `[2, 5, 10, 20, 30, 40, 50]`을 모두 보고한다.
- rank `30/40/50`은 W6의 모든 fold가 사전 grid 상한 `20`을 선택한 뒤 정한 **진단용 확장**이다. W6의 주 결과나 hyperparameter를 소급 변경하지 않는다.
- pathway 입력 panel은 결과 확인 전에 동결한 MSigDB Hallmark `2026.1.Hs`의 KRAS signaling up/down, E2F, G2M, apoptosis, EMT와 8개 immediate-early marker의 합집합이다.
- pathway 변형도 예측 target은 전체 32,738 genes이다. panel symbol을 데이터셋과 교집합한 뒤, outer-training control의 평균·분산 조건만 적용한다.
- W6 component enrichment는 각 fold/component/direction의 top-50 loading symbol에 대해 전체 dataset symbol을 배경으로 hypergeometric test하고, 1,200개 검정 전체에 Benjamini–Hochberg 보정을 적용한다.
- Fold response-subspace 안정성은 두 20차원 basis 사이 principal-angle cosine 제곱의 평균으로 기록한다. component 번호별 직접 정렬을 주장하지 않는다.
- lineage는 outer-training에서만 category를 fit한 one-hot으로 추가하며 test-only category는 all-zero이다.
- BRAF/KRAS는 저자 Figshare v3 annotation의 완전성을 확인한 두 binary field만 사용하고, outer-training 평균·표준편차로 변환한다. training에서 상수인 열은 제거한다.
- `B1`, `B4`, W6 `CCLR`, 모든 고정 변형을 함께 보고하되 outer-test 결과로 최적 변형을 선정하거나 주 모델을 교체하지 않는다.
- sensitivity는 predictor로 사용하지 않는다.

### W8 temporal extension

- clean hash-tag cell만 사용하고 condition/time별 최소 10 cells를 주 기술 기준으로 한다.
- 24h core model의 hyperparameter는 고정하며 time-course로 재선택하지 않는다.
- 외부 model transfer는 core와 겹치는 5개 line을 제외한 17개 line에서만 계산한다.
- 시간별 B1은 core outer-training mean이며 B4/CCLR은 고정된 24h fold artifact의 ensemble 예측이다.
- early 3–6h와 late 24–48h heterogeneity 차이는 line 단위 paired bootstrap으로 평가한다.
- time-course source 차이 때문에 순수 시간 인과효과를 주장하지 않는다.

### W9 biological interpretation

- sensitivity, lineage, BRAF/KRAS는 모두 interpretation-only다.
- 연속 변수는 cell-line Spearman, category는 permutation effect size를 사용한다.
- 38개 검정 전체를 Benjamini–Hochberg로 보정한다.
- prediction error association은 response magnitude와 cell support를 서로 보정한 partial Spearman을 함께 보고한다.

### W10 robustness

- bootstrap seed, inclusion threshold 10/20/30, variable-gene filter 1k/3k/5k/10k, response rank, extreme sensitivity 제거를 비교한다.
- 각 조건 20-cell subsampling은 5개 고정 seed로 반복한다.
- leave-one-lineage-out은 새로운 lineage에 대한 탐색적 외삽 검사이며 주 split을 대체하지 않는다.
- split-half pseudobulk 차이로 full-target measurement floor와 PCC를 근사한다.

### W11 gated distribution extension

- W10 완료 후에만 수행하며, 실행 전 gate는 외부 17-line에서 B4가 B1보다 Energy와 sliced-Wasserstein 중 적어도 하나를 paired CI로 개선하는지 확인한다.
- 2,000 training-derived genes, 10 control-derived PCs를 사용한다.
- 모든 control cell에 같은 line-level predicted shift를 적용하므로 covariance/shape와 paired-cell trajectory는 주장하지 않는다.
- 분포 metric은 B4/CCLR의 tuning이나 model selection에 사용하지 않는다.

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
- `RMSE gain vs B4 = RMSE(B4) - RMSE(CCLR)`
- `NRMSE gain vs B4 = NRMSE(B4) - NRMSE(CCLR)`

양의 gain은 기준 모델보다 낫다는 뜻이다. 프로젝트의 맥락 예측 주장은 raw PCC만이 아니라 gain과 residual metric을 함께 근거로 한다.

## 집계와 불확실성

- 세포 수 가중치 없는 cell-line macro-average가 주 집계다.
- 모델 비교는 같은 held-out cell line의 paired difference로 계산한다.
- 기본 CI는 cell line을 2,000회 paired bootstrap한 percentile 95% CI다.
- 기본 bootstrap seed는 `20260827`이다.
- fold별 값과 세포주별 값을 모두 보존한다.

## 사전 판정 규칙

- B1 대비 개선 CI가 0을 포함하면 일관된 개선이라고 주장하지 않는다.
- CCLR의 B4 대비 paired gain CI가 0 이하이면 CCLR 우위를 주장하지 않는다.
- 일부 lineage에서만 나타나는 개선은 전체 일반화 성능으로 표현하지 않는다.
- split, primary metric, top-k, bootstrap seed는 test 결과를 본 뒤 바꾸지 않는다.
- 규칙 변경이 불가피하면 이전 결과와 분리된 protocol version 및 이유를 남긴다.
