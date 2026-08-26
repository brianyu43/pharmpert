# Data Audit Findings — Probe 0

작성일: 2026-08-26

## 확인된 사실

- pertpy 1.3.0 loader 다운로드 및 로드 성공
- 파일 크기: 1,459,410,830 bytes
- SHA-256: `94a7240047b9ce6822dc2aa1e7a66c1fd8b00be842e5a95fc72cbeb6f2834d2c`
- AnnData shape: 182,875 × 32,738
- `adata.X`: CSC sparse `float32`, 711,380,841 nonzero entries
- `.layers`: 없음
- `.raw`: 없음
- nonzero 100,000개 표본은 모두 정수형 값
- clean time-course hashtag cell은 13,713개로 논문 수치와 일치

## 핵심 발견

단순히 다음 조건을 적용하면 99-line broad experiment가 아니라 여러 24시간 실험이 섞인다.

```text
time == "24"
perturbation in {"control", "Trametinib"}
cell_quality == "normal"
```

이 필터에는 169개 cell line이 들어온다. 그러므로 프로젝트 청사진의 약 99개 core cohort에 바로 사용할 수 없다.

또한 time-course cell의 `time`은 개별 시간이 아니라 문자열 `"3, 6, 12, 24, 48"`로 저장되어 있다. 실제 DMSO/trametinib/untreated 및 시간은 `hash_tag`로 복원해야 한다.

## 첫 probe에서 남은 판정 게이트

Stage 1의 다음 작업은 저자 코드와 Figshare manifest에서 다음 원 실험을 명시적으로 찾아 AnnData cell과 연결하는 것이다.

- `Trametinib_24hr_expt3`
- `DMSO_24hr_expt3`

이 매핑이 검증되기 전에는 24시간 cohort, 97-line threshold, dose를 동결하지 않기로 했다. 아래 교차검증으로 이 게이트를 해소했다.

## Figshare/저자 코드 교차검증 결과

- Supplementary Table 2: experiment 3은 99개 cell line 입력
- Supplementary Table 3: trametinib 농도는 0.1 µM
- 저자 코드의 `trametinib_24hr_expt3` 비교:
  - `control_1 = DMSO_24hr_expt3`
  - `control_2 = DMSO_6hr_expt3`
  - `treat_1 = Trametinib_24hr_expt3`
- 저자 코드와 동일하게 `cell_quality == "normal"`을 적용한 결과:
  - DMSO 24h only 대 trametinib 24h, 양 조건 각각 20개 이상: 94 lines
  - DMSO 6h+24h pooled 대 trametinib 24h, 양 조건 각각 20개 이상: 97 lines

따라서 주 예측은 시간 일치 94-line cohort로 정의하고, 97-line pooled-control 결과는 논문 재현 및 민감도 분석으로 분리한다. 이 선택은 모델 성능을 보기 전에 이루어졌다.

## Cell-line feature 파일 감사

Figshare v3의 `all_CL_features.rds`(3,097,885 bytes, MD5
`07ca14ced15e00468d08ffd7b145a871`)를 공식 API file ID `23322536`에서 받았다.
저자 코드와 일치하게 `Trametinib_24hr_expt3` 객체의 `sens=1-AUC_avg`,
BRAF/KRAS/HRAS/NRAS hotspot field와 `metadata` 객체의 Disease/Subtype을 사용한다.

- strict 94 DepMap ID 전부 정확히 한 번 연결
- Disease 21개, 결측 없음
- 결합 sensitivity 94/94, 네 mutation field 94/94 완전
- sensitivity/mutation은 baseline predictor가 아니라 탐색·해석 전용
