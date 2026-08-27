# 약섭동 프로젝트 설계서

> **프로젝트명:** 약섭동 — 암세포 맥락에 따른 trametinib 전사 반응 예측  
> **버전:** Scope Lock v1.1  
> **시작일:** 2026-08-27  
> **최종 동결일:** 2026-12-15  
> **최종 산출물:** 재현 가능한 코드 저장소 + 분석 결과 + 핵심 그림 + 기술 보고서

---

## 0. 한 문장으로 고정한 프로젝트

**동시에 측정한 vehicle-control(DMSO) 암세포주의 전사 상태를 이용해, 보지 못한 세포주에서 MEK1/2 억제제 trametinib이 24시간 후 일으킬 유전자 발현 변화를 예측하고, 그 반응이 3–48시간 동안 공통 표적 반응에서 세포주별 증식 억제 반응으로 어떻게 분화하는지 분석한다.**

수식으로는 다음 문제다.

\[
\text{vehicle-control state } c_l
\quad\longrightarrow\quad
\widehat{\Delta}_{l,24h}
\]

\[
\Delta_{l,t}
=
\operatorname{Expr}(l,\text{trametinib},t)
-
\operatorname{Expr}(l,\text{control},t)
\]

여기서 \(l\)은 세포주, \(t\)는 처리 시간이다. 핵심은 **학습에 사용하지 않은 세포주**의 반응을 맞히는 것이다.

---

## 1. 최종 범위 결정

### 1.1 반드시 수행할 핵심 범위

| 축 | 최종 결정 |
|---|---|
| 세포 | 암세포주 패널 |
| 약물 | **Trametinib 한 종류만** |
| 섭동 | MEK1/2 억제 |
| 핵심 데이터 | 약 99개 세포주의 DMSO 대 trametinib, 24시간 처리 데이터 |
| 핵심 예측 대상 | 세포주별 24시간 전사 반응 \(\Delta_{l,24h}\) |
| 일반화 조건 | 처리 세포를 전혀 보지 않은 **held-out cell line** |
| 시간 확장 | 24개 세포주, 3·6·12·24·48시간 |
| 표현 단위 | 1차: cell-line pseudobulk, 2차: single-cell distribution |
| 주 평가 | 반응 상관, 오차, 경로 점수, 반응 유전자 일치도 |
| 최종 마감 | 2026-12-15 |

원 논문은 24시간 trametinib/DMSO 실험에서 99개 세포주를 다뤘고, 그중 97개는 양 조건에서 각각 최소 20개 세포가 회수되었다. 별도의 시간축 실험은 24개 세포주를 3–48시간의 다섯 시점에서 측정했으며 총 13,713개 세포를 포함한다.[1]

### 1.2 일부러 제외하는 범위

아래는 이 프로젝트의 첫 보고서에서는 하지 않는다.

- 여러 약물의 동시 비교
- dose-response 모델링
- sci-Plex 전체 188개 화합물 분석
- 환자 반응이나 임상 효능의 직접 예측
- 거대 foundation model 또는 virtual cell 구축
- 동일 세포의 실제 전후 궤적을 복원했다는 주장
- 새로운 wet-lab 실험
- 모든 최신 perturbation 모델을 재현하는 대규모 벤치마크

### 1.3 확장 순서

핵심 결과가 완성된 뒤에만 다음 순서로 확장한다.

1. single-cell 분포 수준 평가
2. 세포주 lineage·BRAF/KRAS 상태를 추가한 ablation
3. 다른 MEK 억제제 또는 다른 약물군
4. dose 축
5. 다중 약물 perturbation world model

**12월 보고서 전에는 3번 이후로 넘어가지 않는다.**

---

## 2. 연구 질문과 가설

### 2.1 주 연구 질문

> **보지 못한 암세포주의 vehicle-control(DMSO) 전사 상태는 그 세포주의 trametinib 24시간 반응을 얼마나 예측할 수 있는가?**

### 2.2 부 연구 질문

1. 모든 세포주에 공통적인 MEK 억제 반응은 무엇인가?
2. 세포주마다 달라지는 반응은 어떤 유전자 프로그램으로 구성되는가?
3. 3–48시간 동안 초기 표적 억제와 후기 cell-cycle/viability 반응은 어떻게 분리되는가?
4. 예측이 잘되는 세포주와 실패하는 세포주는 무엇이 다른가?
5. 관측된 반응 성분은 trametinib sensitivity, lineage, BRAF/KRAS 맥락과 연결되는가?

### 2.3 사전 가설

- **H1 — 공통 반응:** 전체 세포주의 평균 반응만 사용해도 EGR1, DUSP6 등 초기 MAPK 표적 억제는 상당 부분 포착된다.
- **H2 — 맥락의 추가 정보:** baseline 전사 상태를 사용한 모델은 `no-change`와 `global mean response`보다 held-out 세포주의 반응을 더 잘 예측한다.
- **H3 — 시간적 분화:** 초기 3–6시간 반응은 세포주 간 공통성이 크고, 후기 12–48시간 반응은 cell-cycle arrest와 viability 차이 때문에 더 이질적이다.
- **H4 — 생물학적 정합성:** 주요 반응 축은 MAPK signaling, E2F targets, G2M checkpoint 및 약물 감수성과 연관된다.

### 2.4 어떤 결과도 연구 결과가 된다

- 주 모델이 단순 기준선을 확실히 이기면: **baseline context가 약물 반응 이질성을 설명한다.**
- 단순 평균 반응이 주 모델과 비슷하거나 더 좋으면: **이 데이터 규모에서는 복잡한 맥락 모델보다 공유된 약물 반응이 지배적이다.**
- 특정 lineage에서만 개선되면: **context 정보의 유용성이 계층별로 다르다.**

결과를 긍정적으로 만들기 위해 평가 기준이나 데이터 분할을 사후 변경하지 않는다.

---

## 3. 데이터와 근거 자료

### 3.1 원자료

주 데이터는 McFarland et al.의 MIX-seq 연구다.[1]

- **24시간 광범위 패널:** 약 99개 암세포주, DMSO 대 trametinib
- **시간축 패널:** 24개 세포주, 3·6·12·24·48시간
- **관측값:** single-cell RNA-seq UMI count matrix와 cell-line/condition/time metadata
- **부가 정보:** 약물 감수성 및 세포주 특성 일부

원자료는 Figshare에 공개되어 있고,[2] 저자 분석 코드는 `broadinstitute/mix_seq_ms`에 공개되어 있다.[3] 전처리된 AnnData는 scPerturb/pertpy를 통해 불러올 수 있다.[4,5]

### 3.2 가장 먼저 확인해야 할 항목

데이터를 받자마자 아래를 데이터 사전으로 기록한다.

- `adata.X`가 raw count인지 정규화 값인지
- raw count가 `layers` 또는 `raw`에 별도로 존재하는지
- cell-line 열 이름
- drug/perturbation 열 이름
- treatment time 열 이름과 단위
- control 표기 방식: DMSO, vehicle, untreated 등
- trametinib dose와 단위
- batch/pool/replicate 정보
- tissue/lineage 정보
- sensitivity와 mutation 정보의 출처
- gene ID가 symbol인지 Ensembl ID인지

**Trametinib 농도는 기억이나 2차 문헌으로 고정하지 않는다. 원 metadata와 Supplementary Data 3에서 확인한 값을 `data_dictionary.md`에 기록한다.**

### 3.3 빠른 데이터 로드 예시

```python
import pertpy as pt

adata = pt.dt.mcfarland_2020()
print(adata)
print(adata.obs.columns.tolist())
print(list(adata.layers.keys()))
print(adata.raw)
```

pertpy 로더가 실패하거나 필요한 metadata가 누락되면 Figshare 원자료와 scPerturb의 `McFarlandTsherniak2020.h5ad`를 직접 사용한다.

### 3.4 데이터 단위에 대한 중요한 제한

single-cell 데이터라고 해서 각 세포가 독립적인 생물학적 반복은 아니다. 이 프로젝트의 일반화 단위는 **세포 하나가 아니라 세포주 하나**다.

따라서:

- 학습·검증·시험 분할은 반드시 cell line 단위로 한다.
- 세포를 무작위로 train/test에 나누지 않는다.
- pseudobulk는 안정적인 반응 요약으로 사용한다.
- 세포 수가 많다는 이유만으로 과도하게 작은 p-value를 만들지 않는다.
- 불확실성은 cell-line bootstrap과 cell subsampling으로 평가한다.
- 독립 biological replicate가 없는 비교에서는 유전자별 검정을 확정적 인과 증거처럼 쓰지 않는다.[8]

---

## 4. 분석 대상의 수학적 정의

세포주 \(l\), 조건 \(a\), 시간 \(t\), 유전자 \(g\)에 대해 세포별 count를 \(x_{i,g}\)라 하자.

### 4.1 Pseudobulk 표현

조건별 UMI를 합산한 뒤 library-size normalization과 log 변환을 적용한다.

\[
B_{l,a,t,g}
=
\log\left(
1+s\frac{\sum_{i\in(l,a,t)}x_{i,g}}
{\sum_{g'}\sum_{i\in(l,a,t)}x_{i,g'}}
\right)
\]

여기서 \(s\)는 예를 들어 \(10^6\)인 scaling factor다. 실제 구현에서는 CPM/logCPM 또는 검증된 pseudobulk 정규화 방식을 하나로 고정한다.

### 4.2 반응 벡터

\[
\Delta_{l,t,g}
=
B_{l,\mathrm{trametinib},t,g}
-
B_{l,\mathrm{control},t,g}
\]

24시간 핵심 목표는 다음이다.

\[
\widehat{\Delta}_{l,24h}
=
f(c_l)
\]

여기서 \(c_l\)은 해당 세포주의 동시에 측정한 vehicle-control(DMSO) 상태에서 얻은 특징이다. 이는 동일 세포의 처리 전후 종단 측정이 아니므로, 보고서에서는 실제 paired trajectory 또는 엄밀한 의미의 처리 전 상태로 표현하지 않는다.

### 4.3 예측 가능 정보의 경계

시험 세포주에 대해 사용할 수 있는 정보:

- 그 세포주의 control expression
- 미리 정의된 tissue/lineage 정보: ablation에서만
- 공개 mutation 정보: ablation에서만

시험 세포주에 대해 사용할 수 없는 정보:

- trametinib-treated expression
- trametinib sensitivity를 주 모델의 입력으로 사용하는 것
- 시험 세포주의 반응으로 계산한 PCA 또는 feature selection 결과

---

## 5. 분석 파이프라인 전체 구조

```text
원자료 확보
   ↓
metadata·count layer 감사
   ↓
QC 및 세포주×조건×시간 표 작성
   ↓
cell-line pseudobulk 생성
   ↓
24h 반응 행렬 Δ 생성
   ↓
cell-line 단위 nested CV 고정
   ↓
단순 기준선 B0–B4
   ↓
주 모델: context-conditioned low-rank response
   ↓
ablation·불확실성·오류 분석
   ↓
24개 세포주 시간축 분석
   ↓
생물학적 해석
   ↓
그림·표 동결
   ↓
최종 기술 보고서 및 재현성 점검
```

---

## 6. 단계별 실행 계획과 종료 기준

### Stage 0. 프로젝트 고정과 저장소 생성

#### 목적

연구 질문과 완료 기준을 먼저 고정해 분석 중 범위가 계속 커지는 것을 막는다.

#### 작업

- 이 문서를 저장소의 `docs/project_blueprint.md`로 복사
- Git 저장소 생성
- Python 환경 생성 및 lockfile 작성
- 디렉터리 구조 생성
- 랜덤 시드, 결과 파일명 규칙, 실험 ID 규칙 정의
- `scope.md`에 포함/제외 범위 기록

#### 산출물

- `README.md`
- `pyproject.toml` 또는 `environment.yml`
- `uv.lock` 또는 동등한 lockfile
- `docs/scope.md`
- 빈 분석 디렉터리 구조

#### 종료 기준

다른 사람이 저장소를 clone한 뒤 환경을 만들고 빈 smoke test를 실행할 수 있다.

---

### Stage 1. 데이터 확보와 데이터 감사

#### 목적

모델을 만들기 전에 실제 데이터 구조와 분석 가능한 표본 수를 확정한다.

#### 작업

1. pertpy/scPerturb AnnData 다운로드
2. 원 Figshare 파일 목록과 비교
3. count layer 확인
4. metadata 열과 값의 frequency table 생성
5. trametinib, DMSO, untreated 조건만 추출
6. `cell_line × condition × time × batch` 세포 수 표 생성
7. dose와 단위를 원자료에서 검증
8. 결측 metadata와 중복 cell barcode 검사
9. 파일 checksum과 다운로드 URL 기록

#### 산출물

- `data_manifest.csv`
- `data_dictionary.md`
- `cell_count_matrix.csv`
- `00_data_audit.ipynb`
- 감사 로그

#### 종료 기준

아래 질문에 모두 답할 수 있어야 한다.

- 24시간 핵심 분석에 포함 가능한 세포주는 몇 개인가?
- 각 세포주에서 control/treated 세포 수는 몇 개인가?
- 시간축의 모든 시점에서 충분한 세포가 있는 세포주는 몇 개인가?
- raw count는 정확히 어디에 있는가?
- trametinib dose는 무엇인가?
- batch/pool을 어떻게 처리할 것인가?

#### 실패 시 대안

- pertpy metadata가 부족하면 Figshare 원자료 사용
- raw count를 찾지 못하면 scPerturb h5ad의 count layer 확인
- 특정 조건의 세포 수가 지나치게 적으면 해당 세포주를 사전 기준에 따라 제외하고 이유를 기록

---

### Stage 2. QC와 원 논문의 최소 재현

#### 목적

새 모델보다 먼저 데이터가 원 연구의 핵심 현상을 재현하는지 확인한다.

#### 필수 QC

- cell-line·condition·time별 세포 수
- library size와 detected genes
- mitochondrial fraction: raw gene annotation이 있을 때
- control 상태의 PCA/UMAP
- cell-line identity가 주 변동을 설명하는지
- DMSO control의 시간별 변화
- trametinib 처리 후 반응 크기 \(\|\Delta_l\|\)

#### 최소 재현 목표

- 평균적으로 EGR1, DUSP6 등 MAPK immediate-response gene의 억제 확인
- trametinib sensitivity와 반응 축의 연관 방향 확인
- 후기 시점에서 E2F/G2M/cell-cycle 프로그램 변화 확인
- 시간축 데이터의 3–48시간 구조 확인

원 논문은 공유된 MAPK 반응과 세포주별 viability/cell-cycle 반응을 보고했으며, 24시간 반응 PCA가 감수성 및 BRAF/KRAS 맥락과 연결됨을 보였다.[1]

#### 산출물

- `01_qc.ipynb`
- `qc_summary.csv`
- QC 그림 초안
- `reproduction_checklist.md`

#### 종료 기준

원 논문의 방향과 명백하게 모순되는 결과가 없고, 있다면 데이터 버전·정규화·조건 필터 차이로 설명할 수 있다.

---

### Stage 3. 분석 코호트와 반응 행렬 고정

#### 목적

모델링 전에 inclusion rule과 target matrix를 동결한다.

#### 24시간 핵심 코호트

- 주 분석은 시간 일치 `DMSO_24hr_expt3` 대 `Trametinib_24hr_expt3`이며, `cell_quality == normal`에서 양 조건 각각 **20 cells 이상**인 94개 세포주를 사용
- 원 논문/저자 코드 재현 분석은 `DMSO_6hr_expt3 + DMSO_24hr_expt3` pooled control을 사용하며 같은 threshold에서 97개 세포주를 포함
- pooled-control 결과는 주 분석을 대체하지 않고 DMSO time effect와 함께 민감도 분석으로 보고
- 특정 세포주 제외 이유를 표로 기록
- 주 분석과 민감도 분석의 threshold를 분리

#### 시간축 코호트

- 기본적으로 cell-line×time×condition별 10 cells 이상
- 20 cells 이상 조건에서도 결론이 유지되는지 민감도 분석
- DMSO time effect가 미미하다는 것을 확인한 뒤에만 control pooling

#### 유전자 필터

- training fold에서만 저발현 유전자 제거
- control 기반 feature PCA도 training fold에서만 fit
- 결과 해석용 고정 pathway panel은 분석 전에 정의

#### 고정 pathway panel

- MAPK/KRAS signaling
- E2F targets
- G2M checkpoint
- apoptosis
- epithelial–mesenchymal transition
- immediate early response

#### 산출물

- `processed/pseudobulk_24h.parquet`
- `processed/response_24h.parquet`
- `processed/pseudobulk_timecourse.parquet`
- `config/cohort.yaml`
- `config/genesets.yaml`

#### 종료 기준

같은 입력 파일과 config에서 동일한 반응 행렬이 재생성된다.

---

### Stage 4. 데이터 분할과 평가 프로토콜 동결

#### 목적

시험 데이터에 맞춰 모델이나 특징을 조정하는 누수를 막는다.

#### 분할 원칙

- **cell 단위 random split 금지**
- outer split: cell line 기준 5-fold GroupKFold
- inner split: 학습 cell line 안에서 hyperparameter 선택
- lineage가 한 fold에 과도하게 몰리지 않도록 가능한 범위에서 균형화
- 모든 preprocessing, PCA, response basis, feature selection은 outer training fold 안에서 학습

#### 시간축 분석

24개 세포주이므로 leave-one-cell-line-out 또는 6-fold grouped CV를 사용한다. 단일 16/4/4 분할보다 fold 평균과 불확실성을 보고한다.

#### 사전 등록할 평가 지표

#### 주 지표

1. **PCC-Δ:** 세포주별 predicted/observed gene response Pearson correlation
2. **Spearman-Δ:** 반응 유전자 순위 상관
3. **NRMSE-Δ:** 반응 크기로 정규화한 RMSE
4. **Pathway score MAE/correlation**
5. **Signed top-k overlap:** 가장 강한 상·하향 유전자의 방향 일치

#### 맥락 추가 정보 지표

공통 반응만 맞힌 것과 세포주별 이질성을 예측한 것을 분리하기 위해 outer-training fold의 평균 반응 \(\mu_{train}\)을 기준으로 다음 잔차를 계산한다.

\[
r_l=\Delta_l-\mu_{train},
\qquad
\widehat r_l=\widehat\Delta_l-\mu_{train}
\]

- **PCC-context:** \(\widehat r_l\)과 \(r_l\)의 유전자별 Pearson correlation
- **RMSE gain vs B1:** \(\operatorname{RMSE}(B1)-\operatorname{RMSE}(model)\); 양수일수록 주 모델이 유리
- B0처럼 상수 예측인 모델의 Pearson/Spearman은 `NA`로 두고 0으로 대체하지 않음
- NRMSE 분모, top-k의 k, 상수 벡터 처리 규칙은 `docs/evaluation_protocol.md`에서 모델 결과를 보기 전에 고정

#### 보조 지표

- response norm calibration
- Common-DEG 또는 thresholded response gene overlap
- single-cell 확장 시 Energy distance와 Wasserstein distance

Perturbation benchmark들은 MSE, response correlation, Wasserstein, KL, Energy distance, Common-DEG 같은 서로 다른 관점의 지표를 함께 사용할 것을 강조한다.[6,7]

#### 집계 방식

- cell 수로 가중한 micro-average가 아니라 **cell-line macro-average**를 주 결과로 사용
- fold별 값과 전체 cell-line bootstrap 95% CI 보고
- 기준선과의 차이는 paired cell-line bootstrap으로 비교

#### 산출물

- `config/splits.yaml`
- 고정 split ID 목록
- `metrics.py`
- `evaluation_protocol.md`

#### 종료 기준

모델 결과를 보기 전에 split과 지표가 Git commit으로 고정되어 있다.

---

### Stage 5. 기준선 구축

복잡한 모델보다 먼저 아래 기준선을 모두 구현한다.

| ID | 모델 | 정의 | 목적 |
|---|---|---|---|
| B0 | No change | \(\widehat\Delta_l=0\) | 약물이 없다고 가정하는 최저선 |
| B1 | Global mean | training lines의 평균 \(\bar\Delta\) | 공유 약물 반응의 강도 |
| B2 | Lineage mean | 같은 lineage의 평균 반응 | 단순한 맥락 효과 |
| B3 | Nearest neighbor | baseline transcriptome이 가장 가까운 학습 세포주의 반응 | 비모수 context baseline |
| B4 | Direct ridge | baseline PC에서 gene-level response를 직접 다중회귀 | 저복잡도 지도학습 |

B2는 fold 안에 동일 lineage의 충분한 학습 세포주가 있을 때만 계산하고, 없으면 B1로 fallback한다.

최근 benchmark에서 복잡한 perturbation 모델이 단순 선형 기준선을 항상 능가하지 못한다는 결과가 있으므로, 이 단계는 형식적인 비교가 아니라 프로젝트의 핵심이다.[6]

#### 산출물

- `03_baselines.ipynb`
- `baseline_metrics.csv`
- 세포주별 예측 파일
- 실행 시간과 parameter 수 표

#### 종료 기준

모든 outer fold에서 B0–B4가 오류 없이 실행되고, 각 cell line의 예측이 저장된다.

---

### Stage 6. 주 모델 — Context-Conditioned Low-Rank Response

#### 선택 이유

세포주는 약 97개 수준이므로 고차원 deep generative model을 주 모델로 삼기에는 context 수가 적다. 대신 반응을 몇 개의 공유 프로그램으로 압축하고, baseline 상태가 각 프로그램의 세기를 예측하게 한다.

#### 모델

학습 세포주의 반응 행렬에 대해

\[
\Delta_l \approx \mu + W s_l
\]

- \(\mu\): 학습 세포주의 평균 trametinib 반응
- \(W\): 반응 프로그램 basis, response PCA 또는 truncated SVD
- \(s_l\): 세포주별 프로그램 점수

control expression에서 context feature를 만든다.

\[
z_l=\operatorname{PCA}_{train}(c_l)
\]

ridge regression으로 반응 점수를 예측한다.

\[
\widehat{s}_l=A z_l+b
\]

최종 반응은

\[
\widehat{\Delta}_l=\mu+W\widehat{s}_l
\]

이다.

#### 학습할 hyperparameter

- control PCA 차원: 예 5, 10, 20, 30
- response basis 차원: 예 2, 5, 10, 20
- ridge penalty \(\lambda\)
- 저발현 유전자 threshold

모든 선택은 inner CV에서만 한다.

#### 왜 이 모델이 보고서에 적합한가

- global mean baseline을 자연스럽게 포함한다.
- 각 response component를 pathway와 연결할 수 있다.
- 데이터 규모에 비해 parameter 수를 통제할 수 있다.
- 실패 이유를 component별로 해석할 수 있다.
- 향후 시간축·다중 약물 모델로 확장하기 쉽다.

#### 산출물

- `src/yakseopdong/models.py`
- `04_main_model.ipynb`
- fold별 model artifact
- component loading/score 표
- 주 모델 성능 표

#### 종료 기준

- 모든 outer fold의 held-out cell line 예측 완성
- 주 지표와 95% CI 생성
- 동일 seed/config에서 결과 재현
- B0–B4와 공정한 비교 완료

#### 실행 결과 (2026-08-27)

- strict 94-line, outer 5-fold / inner 4-fold에서 CCLR held-out 예측 94개를 완성했다.
- CCLR RMSE는 `0.322383`, B1 대비 gain은 `0.001367`이다.
- B4 대비 paired RMSE gain은 `-0.000348` (95% CI `-0.000584–-0.000105`)로,
  저차원 response basis가 direct ridge를 능가한다는 가설은 지지되지 않았다.
- 동일 seed/config prediction SHA-256 재현, fold artifact 5개, component loading/score 표,
  `04_main_model.ipynb`를 생성했다.
- 모든 fold가 response rank 20과 alpha 100을 선택한 경계값 현상은 결과를 본 뒤 W6
  grid를 바꾸지 않고 Stage 7 ablation 대상으로 이관한다.

---

### Stage 7. Ablation과 누수 감사

#### 필수 ablation

1. 평균 반응 \(\mu\)만 사용
2. context-dependent component만 추가
3. control PCA 차원 변화
4. response rank 변화
5. pathway gene panel 대 전체 필터 유전자
6. tissue/lineage feature 추가
7. BRAF/KRAS feature 추가: metadata 신뢰성이 확인될 때만

#### 누수 체크리스트

- test line의 treated expression이 PCA fit에 들어갔는가?
- response-variable gene 선택에 test response가 사용됐는가?
- test line의 sensitivity가 predictor에 들어갔는가?
- 동일 cell line이 다른 pool 이름으로 train/test에 중복됐는가?
- normalization parameter가 전체 데이터에서 학습됐는가?
- hyperparameter를 outer test 성능을 보고 선택했는가?

#### 산출물

- `leakage_audit.md`
- `ablation_metrics.csv`
- 모델 복잡도 대 성능 그림
- `05_ablation.ipynb`
- pathway enrichment와 fold response-subspace 안정성 표

#### 종료 기준

누수 항목이 모두 `PASS`이거나, 예외가 있다면 결과에서 명시적으로 분리된다.

#### 실행 결과 (2026-08-27)

- v1.2 protocol에서 B1/B4/W6 CCLR과 13개 fixed low-rank 변형, 총 16개 비교군을
  같은 94-line outer 5-fold에 평가했다.
- 모든 fixed 변형은 B1보다 일관되게 나았지만 B4를 확실히 이긴 변형은 없었다.
- full d20/r20의 B4 대비 gain은 `-0.000315` (95% CI `-0.000612–-0.000044`)였다.
- pathway, lineage, BRAF/KRAS 변형의 B4 대비 gain CI는 모두 0을 포함했다.
- rank를 30/40/50으로 늘려도 B4보다 나았다는 근거가 없어 W6 rank 상한 뒤에 큰
  미탐색 이득이 있다는 가설은 지지되지 않았다.
- MSigDB Hallmark 2026.1.Hs 기반 panel은 1,039 symbols 중 1,001 symbols가 데이터에
  매핑됐고, fold별 training filter 뒤 917–924 columns가 사용됐다.
- W6 response subspace의 fold-pair mean squared cosine은 평균 `0.676`이었다.
- 1,504개 line-level metric, paired bootstrap CI, parameter count, pathway enrichment,
  10개 fold-pair 안정성 값을 저장했다. 누수 항목은 모두 PASS다.

---

### Stage 8. 시간축 분석

#### 목적

24시간 예측 결과를 단일 시점 성능표로 끝내지 않고, 반응이 어떻게 형성되는지 설명한다.

#### 분석 질문

- 3시간에 이미 나타나는 공통 프로그램은 무엇인가?
- 6–12시간 사이에 새로 등장하는 반응은 무엇인가?
- 24–48시간에 sensitivity와 연결되는 프로그램은 무엇인가?
- 초기 반응 크기가 후기 cell-cycle arrest를 예고하는가?

#### 권장 방법

1. 각 시점의 \(\Delta_{l,t}\) 계산
2. 시간별 global mean response와 cell-line deviation 분리
3. response PCA 또는 공통 low-rank basis 추정
4. component score를 `log(1 + time)` 위에서 spline 또는 단순 선으로 시각화
5. early 3–6h와 late 24–48h의 세포주 간 분산 비교
6. 24시간 주 모델이 시간축에서 어느 시점부터 유효해지는지 평가

시간축에서는 복잡한 recurrent model을 만들지 않는다. 목적은 **동역학의 해석과 24시간 결과의 맥락화**다.

#### 산출물

- `06_temporal.ipynb`
- 시간별 반응 행렬
- component trajectory 그림
- early/late heterogeneity 표

#### 종료 기준

초기 공통 반응과 후기 context-specific 반응을 최소 두 개의 정량 지표와 두 개의 그림으로 구분할 수 있다.

---

### Stage 9. 생물학적 해석

#### 분석 대상

- EGR1, ETV4/5, DUSP4/5/6, SPRY2/4
- E2F targets
- G2M checkpoint
- apoptosis
- MAPK/KRAS signaling
- MCM 계열 등 cell-cycle 관련 유전자

#### 연결 분석

- response component score 대 trametinib sensitivity
- component score 대 lineage
- component score 대 BRAF/KRAS 상태
- 예측 오차 대 baseline 상태
- 예측 오차 대 관측 반응 크기
- 예측 오차 대 세포 수: 측정 불확실성과 모델 실패 구분

#### 통계

- cell-line 단위 Spearman correlation
- bootstrap 95% CI
- category 비교 시 permutation test 또는 effect size 우선
- 다중 검정은 탐색 분석임을 명시하고 FDR 사용
- 독립 반복이 부족한 single-cell gene-level p-value를 과도하게 해석하지 않음

#### 산출물

- `biological_validation.csv`
- pathway enrichment 표
- best/median/worst 예측 사례
- 오류 유형 분류

#### 종료 기준

모델 성능이 왜 나왔는지를 최소 세 가지 생물학적·측정적 요인으로 설명할 수 있다.

---

### Stage 10. 강건성 분석

#### 필수 분석

- cell-line bootstrap
- cell subsampling: 각 조건의 세포 수를 동일하게 줄인 뒤 반복
- inclusion threshold 10/20/30 cells 비교
- gene filter 변화
- response rank 변화
- Pearson 대 Spearman 결과 비교
- 하나의 lineage를 통째로 제외하는 leave-one-lineage-out 탐색
- 극단적으로 민감하거나 저항성인 세포주 제거 후 결론 확인

#### 주요 구분

관측 반응이 noisy해서 예측이 어려운 것과, context 모델이 구조를 못 배우는 것을 분리한다.

이를 위해 같은 세포주·조건의 세포를 두 집단으로 반복 분할해 **관측 가능한 내부 재현성 상한선**을 근사한다.

#### 산출물

- `robustness_metrics.csv`
- noise ceiling 그림
- 주요 결론별 robustness table

#### 종료 기준

주 결론이 특정 threshold, seed, 소수 세포주 하나에 의존하지 않는다. 의존한다면 보고서의 제한점으로 승격한다.

---

### Stage 11. 선택적 single-cell 분포 확장

이 단계는 **핵심 모델과 강건성 분석이 모두 끝난 경우에만** 수행한다.

#### 가능한 최소 확장

- control cell들을 baseline PCA space에 표현
- 예측된 pseudobulk delta를 cell-level expression에 이동량으로 적용
- observed treated distribution과 비교
- PCA space에서 Energy distance와 Wasserstein distance 계산

#### 하지 않을 것

- scGen, CPA, diffusion model을 처음부터 대규모로 구현
- distribution metric 하나만 좋아지도록 모델 튜닝
- paired cell trajectory를 복원했다고 주장

#### 종료 기준

2주 이내에 명확한 추가 그림을 만들지 못하면 본문이 아니라 향후 연구로 남긴다.

---

### Stage 12. 결과 동결과 보고서 작성

#### 결과 동결 전에 필요한 파일

- `final_metrics.csv`
- `final_predictions.parquet`
- `final_cell_lines.csv`
- `figure_manifest.csv`
- `table_manifest.csv`
- `limitations.md`
- 최종 config와 commit hash

#### 보고서의 중심 주장 후보

> **Trametinib의 공유된 MAPK 표적 반응은 평균 반응만으로도 상당 부분 예측되지만, 후기 cell-cycle/viability 반응은 baseline cellular context에 의해 달라지며 저차원 맥락 모델이 그 이질성의 일부를 회수한다.**

실제 결과가 이 주장과 다르면 결과에 맞추어 문장을 바꾼다.

---

## 7. 실험 매트릭스

| 실험 | 질문 | 데이터 | 모델/비교 | 필수 여부 |
|---|---|---|---|---|
| E0 | 데이터가 원 논문의 핵심 현상을 재현하는가? | 24h + time-course | PCA, pathway score | 필수 |
| E1 | 공유 반응만으로 얼마나 예측되는가? | 24h broad panel | B0, B1 | 필수 |
| E2 | 단순 context가 도움이 되는가? | 24h broad panel | B2–B4 | 필수 |
| E3 | 저차원 주 모델이 기준선을 이기는가? | 24h broad panel | CCLR 대 B0–B4 | 필수 |
| E4 | 어떤 정보가 개선을 만드는가? | 24h broad panel | ablation | 필수 |
| E5 | 반응이 시간에 따라 어떻게 분화하는가? | 24-line time-course | temporal PCA/pathway | 필수 |
| E6 | 결과가 sampling/QC에 강건한가? | 양 데이터 | bootstrap/subsample | 필수 |
| E7 | 단일세포 분포도 개선되는가? | time-course 또는 24h | Energy/Wasserstein | 선택 |

CCLR은 이 문서의 `Context-Conditioned Low-Rank Response` 모델을 뜻한다.

---

## 8. 성공 기준

### 8.1 최소 성공

- 데이터와 metadata를 완전히 감사했다.
- 24시간 broad panel의 반응 landscape를 재현했다.
- held-out cell-line 기준 B0–B4와 주 모델을 비교했다.
- 시간축의 early/late 반응을 분석했다.
- 모든 결과가 한 명령 또는 명시된 명령 순서로 재현된다.
- 보고서가 null result도 정직하게 설명한다.

### 8.2 강한 성공

아래 중 둘 이상을 만족한다.

- 주 모델이 B1 global mean보다 PCC-Δ와 pathway metric에서 일관되게 우수
- 개선량의 paired bootstrap 95% CI가 0을 넘음
- response component가 sensitivity/lineage/MAPK genotype과 정합적인 관계
- early shared response와 late heterogeneous response의 분산 차이가 강건함
- best/median/worst 사례에서 예측 성공·실패를 생물학적으로 설명 가능

### 8.3 실패로 간주하지 않는 결과

주 모델이 기준선을 이기지 못해도 아래가 충족되면 완성된 연구다.

- 누수 없는 엄격한 context split
- 강한 기준선과 공정한 비교
- null result의 원인 분석
- 데이터 규모·노이즈 상한선·일반화 한계의 정량화

---

## 9. 최종 그림과 표

### Figure 1. 연구 설계

- 99-line 24h core와 24-line time course의 관계
- control state → response prediction 구조
- held-out cell-line split 강조

**전달할 주장:** 세포가 아니라 세포주가 일반화 단위다.

### Figure 2. 데이터 품질과 코호트

- cell-line×condition cell counts
- QC distributions
- 최종 포함/제외 흐름

**전달할 주장:** 성능 차이가 표본 수 불균형만으로 설명되지 않는다.

### Figure 3. 관측된 24시간 반응 구조

- global mean response
- cell-line response PCA
- MAPK 및 cell-cycle pathway score

**전달할 주장:** 공유 반응과 context-specific 반응이 동시에 존재한다.

### Figure 4. Held-out cell-line benchmark

- B0–B4와 주 모델의 metric distribution
- paired difference와 95% CI

**전달할 주장:** baseline context가 단순 평균 이상으로 주는 정보의 크기.

### Figure 5. 예측 사례

- best, median, worst cell line의 observed 대 predicted response
- top up/down genes와 pathway score

**전달할 주장:** 평균 성능 뒤에 어떤 성공과 실패가 있는지.

### Figure 6. 모델 해석

- response component loading
- baseline context feature와 response score 관계
- pathway enrichment

**전달할 주장:** 모델이 어떤 반응 프로그램을 예측했는지.

### Figure 7. 시간축

- 3–48h pathway/component trajectories
- early 대 late cross-line variance

**전달할 주장:** 초기 표적 반응에서 후기 운명 반응으로 이질성이 커진다.

### Figure 8. 생물학적 검증과 오류 분석

- response score 대 sensitivity
- lineage/BRAF/KRAS별 분포
- error 대 cell count 및 response magnitude

**전달할 주장:** 예측 가능성과 실패가 생물학적 맥락 및 측정 품질과 연결된다.

### 핵심 표

- Table 1: 분석 코호트와 세포 수
- Table 2: 모델 정의와 parameter 수
- Table 3: 주 benchmark 성능과 95% CI
- Table 4: ablation 결과
- Table 5: robustness 결과
- Supplementary Table: 세포주별 성능

---

## 10. 최종 보고서 목차

### 보고서 제목

**Predicting Context-Dependent Transcriptional Responses to MEK Inhibition across Cancer Cell Lines**

한국어 부제:

**암세포주의 초기 전사 상태로부터 trametinib 반응을 예측하는 저차원 통계 모델**

### 1. Technical Summary

- 질문
- 데이터
- 방법
- 가장 중요한 결과
- 한계와 결론

### 2. Introduction

- 동일 약물에도 세포 맥락에 따라 반응이 달라지는 문제
- single-cell perturbation prediction의 목표
- 기존 모델에서 unseen context generalization이 어려운 이유
- 본 연구가 하나의 약물과 많은 세포주를 선택한 이유
- 연구 질문과 가설

### 3. Data and Cohort Definition

- MIX-seq 자료
- 24h broad panel
- 3–48h time-course
- QC, inclusion rule, metadata
- pseudobulk를 주 분석으로 선택한 이유

### 4. Methods

- 전처리
- response 정의
- cell-line grouped nested CV
- B0–B4
- CCLR 주 모델
- 지표
- bootstrap과 robustness
- pathway 및 sensitivity 연계

### 5. The Observed Trametinib Response Landscape

- 공유 MAPK 반응
- 세포주별 이질성
- sensitivity·lineage와의 관계

### 6. Predicting Responses in Unseen Cell Lines

- benchmark
- 주 모델 개선량
- 세포주별 성능
- calibration

### 7. What the Model Learns

- response components
- baseline context feature
- pathway 해석
- ablation

### 8. Temporal Emergence of Context Dependence

- 3–48h trajectory
- early/late heterogeneity
- 24h 모델의 시간적 의미

### 9. Robustness and Failure Analysis

- subsampling
- threshold
- noise ceiling
- best/median/worst 사례
- null 또는 불안정한 결과

### 10. Discussion

- 얻은 답
- virtual cell/perturbation modeling에 대한 함의
- 왜 단순 기준선이 강하거나 약했는지
- 세포주 데이터의 한계
- 환자·조직·다중 약물로 확장할 때 필요한 것

### 11. Conclusion

한 문단으로 답한다.

> untreated state는 trametinib 반응의 어느 부분을 예측하며, 어느 부분은 여전히 예측하지 못하는가?

### 12. Reproducibility

- 데이터 버전과 URL
- 환경
- 실행 명령
- seed
- commit hash
- 결과 파일 manifest

### 13. References

---

## 11. 저장소 구조

```text
yak-seopdong/
├── README.md
├── pyproject.toml
├── uv.lock
├── Makefile
├── config/
│   ├── cohort.yaml
│   ├── splits.yaml
│   ├── models.yaml
│   └── genesets.yaml
├── data/
│   ├── raw/              # Git 제외
│   ├── interim/
│   └── processed/
├── docs/
│   ├── project_blueprint.md
│   ├── scope.md
│   ├── data_dictionary.md
│   ├── evaluation_protocol.md
│   └── leakage_audit.md
├── notebooks/
│   ├── 00_data_audit.ipynb
│   ├── 01_qc.ipynb
│   ├── 02_response_landscape.ipynb
│   ├── 03_baselines.ipynb
│   ├── 04_main_model.ipynb
│   ├── 05_ablation.ipynb
│   ├── 06_temporal.ipynb
│   └── 06_report_figures.ipynb
├── src/yakseopdong/
│   ├── io.py
│   ├── preprocess.py
│   ├── pseudobulk.py
│   ├── splits.py
│   ├── models.py
│   ├── metrics.py
│   ├── pathways.py
│   └── plots.py
├── tests/
│   ├── test_splits.py
│   ├── test_pseudobulk.py
│   └── test_metrics.py
├── results/
│   ├── figures/
│   ├── tables/
│   ├── predictions/
│   ├── models/
│   └── logs/
└── report/
    ├── final_report.md
    ├── references.bib
    └── supplementary.md
```

### 실행 인터페이스 예시

```bash
make data
make preprocess
make benchmark
make temporal
make figures
make report
make reproduce-all
```

처음에는 notebook으로 탐색해도 되지만, 최종 수치를 만드는 코드는 `src/`와 명령형 pipeline으로 옮긴다.

---

## 12. 2026년 주차별 일정

| 기간 | 단계 | 반드시 끝낼 산출물 |
|---|---|---|
| **8/27–8/30** | W0 범위 고정 | repo, 환경, 디렉터리, 데이터 다운로드 smoke test |
| **8/31–9/6** | W1 데이터 감사 | manifest, dictionary, cell count matrix, dose 확인 |
| **9/7–9/13** | W2 QC·최소 재현 | QC 그림, 원 논문 핵심 반응 재현 |
| **9/14–9/20** | W3 pseudobulk·반응 행렬 | 24h/time-course processed matrices |
| **9/21–9/27** | W4 split·metric 동결 | grouped nested CV, evaluation protocol |
| **9/28–10/4** | W5 기준선 | B0–B4 전체 결과 |
| **10/5–10/11** | W6 주 모델 | CCLR outer-CV 결과 |
| **10/12–10/18** | W7 ablation·누수 감사 | ablation table, leakage audit |
| **10/19–10/25** | W8 시간축 | early/late trajectory 결과 |
| **10/26–11/1** | W9 생물학적 해석 | sensitivity/pathway/error analysis |
| **11/2–11/8** | W10 robustness | bootstrap, subsampling, noise ceiling |
| **11/9–11/15** | W11 선택 확장 | single-cell 분포 결과 또는 명시적 중단 |
| **11/16–11/22** | W12 그림·표 동결 | Figure 1–8, Table 1–5 초안 |
| **11/23–11/29** | W13 보고서 v1 | 전 섹션이 채워진 첫 완성본 |
| **11/30–12/6** | W14 검증·개정 | clean run, 수치 대조, 보고서 v2 |
| **12/7–12/13** | W15 최종화 | README, limitations, release candidate |
| **12/14–12/15** | 최종 동결 | final report, code tag, artifact bundle |

---

## 13. 범위 절단 규칙

### 10월 11일 게이트

B0–B4와 주 모델의 outer-CV가 끝나지 않았다면:

- single-cell generative modeling 전부 삭제
- 외부 mutation/lineage feature 확장 보류
- 24시간 core benchmark 완성에 집중

### 11월 1일 게이트

시간축 분석이 끝나지 않았다면:

- 본문은 24시간 core로 완성
- 시간축은 단순 descriptive figure 또는 supplement로 축소

### 11월 15일 게이트

single-cell 분포 확장이 완전하지 않다면:

- 즉시 중단
- 향후 연구 한 문단으로 이동
- 본문 결과에는 포함하지 않음

### 절대 규칙

- Figure/Table freeze 이후 모델 추가 금지
- 시험 성능을 보고 split 변경 금지
- 새로운 약물 추가 금지
- dose 축 추가 금지
- 결과가 약하다는 이유로 primary metric 변경 금지

---

## 14. 위험 관리표

| 위험 | 조기 신호 | 대응 |
|---|---|---|
| raw count/metadata 불명확 | `adata.X` 의미를 확정 못함 | Figshare 원자료와 저자 코드로 교차검증 |
| cell 수 불균형 | 일부 조건의 세포 수가 매우 적음 | inclusion threshold, subsampling, macro-average |
| 진짜 biological replicate 부족 | gene-level p-value가 과도하게 작음 | pseudobulk effect size, cell-line bootstrap, 제한 명시 |
| leakage | random cell split에서만 성능이 높음 | GroupKFold, fold 내부 전처리 테스트 |
| 모델이 global mean을 못 이김 | B1과 차이가 거의 없음 | null result로 정리, 예측 가능한 component만 분석 |
| 시간축 n=24가 작음 | 복잡한 모델이 불안정 | descriptive low-rank/pathway trajectory로 제한 |
| sensitivity metadata 혼선 | 출처·방향 정의 불명확 | validation 분석에서 제외하거나 출처별 분리 |
| 계산 자원 문제 | sparse matrix가 dense로 변환됨 | sparse 연산, pseudobulk-first, chunk 저장 |
| 보고서가 늦어짐 | 11/22에도 그림이 바뀜 | model freeze, 미완성 확장 삭제 |
| 생물학적 과해석 | 세포주 결과를 환자에게 일반화 | in-vitro limitation을 summary와 discussion에 명시 |

---

## 15. 처음 72시간 작업 체크리스트

### Day 1 — 저장소와 데이터

- [ ] `yak-seopdong` Git 저장소 생성
- [ ] Python 환경과 JupyterLab 설정
- [ ] pertpy로 데이터 로드 smoke test
- [ ] AnnData 구조를 텍스트 파일로 저장
- [ ] 원자료 Figshare와 저자 GitHub 주소 기록
- [ ] raw data를 Git에서 제외
- [ ] 이 문서를 `docs/project_blueprint.md`로 추가

### Day 2 — Metadata 감사

- [ ] `obs` 모든 열과 unique value 요약
- [ ] trametinib/control/time 필터 후보 확인
- [ ] cell-line×condition×time cell count matrix 생성
- [ ] count layer 위치 확인
- [ ] gene ID 종류 확인
- [ ] dose를 Supplementary Data 3 또는 metadata에서 확인
- [ ] batch/pool 정보 확인

### Day 3 — 첫 그림과 첫 결정

- [ ] 24시간 코호트 후보 목록 생성
- [ ] 시간축 24개 세포주 목록 생성
- [ ] cell count heatmap
- [ ] pseudobulk 한 번 생성
- [ ] EGR1/DUSP6 평균 response 확인
- [ ] `data_dictionary.md` 작성
- [ ] inclusion threshold 초안 결정
- [ ] 첫 Git tag `data-audit-v0` 생성

### 72시간 종료 산출물

```text
README.md
docs/project_blueprint.md
docs/data_dictionary.md
data_manifest.csv
cell_count_matrix.csv
notebooks/00_data_audit.ipynb
results/figures/cell_count_heatmap.png
```

---

## 16. 각 주 작업 세션의 기본 리듬

### 월요일

- 주 질문 하나 고정
- 성공 기준과 결과 파일명 선언

### 화–수요일

- 분석 코드 작성
- 중간 표와 그림 생성
- 오류/가정 기록

### 목요일

- robustness 또는 누수 점검
- notebook 결과를 함수/스크립트로 이동

### 금요일

- 주간 결과 1쪽 작성
- 다음 주 blocker와 범위 결정
- Git tag 또는 milestone commit

### 주간 완료 기준

“코드를 많이 썼다”가 아니라 다음 셋을 남긴다.

1. 재실행 가능한 결과 파일
2. 그 결과가 답하는 한 문장
3. 다음 단계로 넘어가도 되는 종료 판정

---

## 17. 최종 재현성 점검

최종 제출 전 깨끗한 디렉터리에서 확인한다.

- [ ] 환경을 처음부터 설치 가능
- [ ] 데이터 다운로드 또는 수동 배치 절차가 명시됨
- [ ] `make reproduce-all` 또는 문서화된 명령 순서가 성공
- [ ] 모든 표의 숫자가 `results/tables` 파일과 일치
- [ ] 모든 그림이 코드로 재생성됨
- [ ] split과 seed가 저장됨
- [ ] test cell line 목록이 공개됨
- [ ] fold 내부 전처리 여부 테스트 통과
- [ ] 데이터 라이선스와 출처 표기
- [ ] Git commit hash와 release tag 기록
- [ ] 개인정보·API key·로컬 절대경로 없음

---

## 18. Definition of Done

아래 조건이 모두 충족되어야 프로젝트가 끝난 것이다.

- [ ] 연구 질문이 한 문장으로 명확하다.
- [ ] Trametinib 한 약물, 24시간 core, 3–48시간 extension 범위를 유지했다.
- [ ] 데이터 사전과 inclusion/exclusion 표가 있다.
- [ ] held-out cell-line 평가를 사용했다.
- [ ] 최소 다섯 개의 기준선과 주 모델을 비교했다.
- [ ] 주 지표, CI, ablation, robustness가 있다.
- [ ] early/late temporal response를 정량화했다.
- [ ] 모델의 생물학적 해석과 실패 분석이 있다.
- [ ] 결과가 약해도 숨기지 않고 설명했다.
- [ ] Figure 1–8과 Table 1–5가 최종 파일로 존재한다.
- [ ] 보고서의 모든 주장에 대응하는 그림·표·코드가 있다.
- [ ] 새 환경에서 재현 실행을 통과했다.
- [ ] `final_report.md`와 Git release가 2026-12-15 이전에 동결되었다.

---

## 19. 이 프로젝트의 최종 포지셔닝

이 프로젝트는 “거대한 virtual cell을 만들었다”는 주장이 아니다. 더 정확한 기여는 다음과 같다.

1. **하나의 명확한 섭동 아래 cellular context generalization을 엄격히 정의한다.**
2. **세포 무작위 분할이 아니라 보지 못한 세포주에 대한 일반화를 평가한다.**
3. **복잡한 모델보다 강한 통계적 기준선과 해석 가능한 저차원 모델을 비교한다.**
4. **24시간 예측과 3–48시간 반응 형성을 하나의 이야기로 연결한다.**
5. **무엇이 예측되고 무엇이 아직 예측되지 않는지를 정직하게 분리한다.**

린메트로가 알고리즘의 정당성을 형식적으로 검증하는 프로젝트라면, 약섭동은 noisy biological system에서 **공유 구조와 맥락 의존성을 통계적으로 분해하는 프로젝트**다.

---

## 20. 참고문헌 및 데이터 링크

1. McFarland, J. M. et al. *Multiplexed single-cell transcriptional response profiling to define cancer vulnerabilities and therapeutic mechanism of action.* Nature Communications 11, 4296 (2020). DOI: [10.1038/s41467-020-17440-w](https://doi.org/10.1038/s41467-020-17440-w)
2. McFarland et al. MIX-seq data, Figshare. [https://figshare.com/articles/dataset/MIX-seq_data/10298696](https://figshare.com/articles/dataset/MIX-seq_data/10298696)
3. Broad Institute. `mix_seq_ms` analysis code. [https://github.com/broadinstitute/mix_seq_ms](https://github.com/broadinstitute/mix_seq_ms)
4. Peidli, S. et al. *scPerturb: harmonized single-cell perturbation data.* Nature Methods 21, 531–540 (2024). DOI: [10.1038/s41592-023-02144-y](https://doi.org/10.1038/s41592-023-02144-y)
5. Heumos, L. et al. *Pertpy: an end-to-end framework for perturbation analysis.* Nature Methods (2026). DOI: [10.1038/s41592-025-02909-7](https://doi.org/10.1038/s41592-025-02909-7)
6. Wei, Y. et al. *Systematic benchmarking of computational methods for single-cell perturbation response prediction.* Nature Methods (2025). DOI: [10.1038/s41592-025-02980-0](https://doi.org/10.1038/s41592-025-02980-0)
7. Wu, J. et al. *PerturBench: Benchmarking machine learning models for cellular perturbation analysis.* arXiv:2408.10609 / NeurIPS Datasets and Benchmarks. [https://arxiv.org/abs/2408.10609](https://arxiv.org/abs/2408.10609)
8. Squair, J. W. et al. *Confronting false discoveries in single-cell differential expression.* Nature Communications 12, 5692 (2021). DOI: [10.1038/s41467-021-25960-2](https://doi.org/10.1038/s41467-021-25960-2)
