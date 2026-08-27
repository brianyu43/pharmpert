# 약섭동 (Yak-seopdong)

암세포주의 vehicle-control(DMSO) 전사 상태로부터 보지 못한 세포주의 24시간 trametinib 전사 반응을 예측하는 재현 가능한 분석 프로젝트다.

## 현재 상태

- Scope Lock: v1.1
- 단계: Stage 1–6 완료; 24h 반응 행렬, B0–B4, CCLR 주 모델 완료
- 핵심 일반화 단위: cell line
- 핵심 비교: B1 global mean response 대비 context 정보의 추가 이득
- GPU: 핵심 분석에는 불필요

공식 experiment 3 원자료에서 strict 94-line 24시간 pseudobulk와 반응 행렬을 생성했다. EGR1/DUSP6 등 사전 지정 MAPK marker의 억제 방향을 재현했고, lineage-aware 5-fold에서 B0–B4와 CCLR을 평가했다. B4 direct ridge는 B1 global mean보다 RMSE를 `0.001715` 개선했다. CCLR도 B1보다 `0.001367` 개선했지만 B4보다는 `0.000348` 나빴다. 맥락의 추가 정보는 검출되지만 강한 개인화 예측이나 임상적 결과는 주장하지 않는다.

## 빠른 시작

필수 도구는 Python 3.13과 [uv](https://docs.astral.sh/uv/)다.

```bash
make install
make smoke
make test
```

데이터 다운로드를 포함한 probe는 명시적으로 실행한다.

```bash
make data-probe
make core-audit
make pseudobulk
make metadata
make landscape
make splits
make baselines
make cclr
make validate-cclr
make notebooks
```

## 문서

- `docs/project_blueprint.md`: 전체 연구 청사진
- `docs/scope.md`: 포함·제외 범위와 연구 질문
- `docs/evaluation_protocol.md`: 분할, 지표, 누수 방지 규칙
- `docs/data_dictionary.md`: Stage 1에서 채울 데이터 계약
- `docs/leakage_audit.md`: fold별 누수 감사표
- `docs/reproduction_checklist.md`: 현재 최소 재현 결과와 남은 항목
- `processed_manifest.csv`: 생성된 Parquet의 shape, 값 열, SHA-256
- `cell_line_annotations.csv`: 저자 제공 disease/sensitivity/mutation의 94-line join
- `split_assignments.csv`, `inner_split_assignments.csv`: 고정 outer/inner cell-line split

## 검증 게이트

Stage 0 완료는 다음 명령이 모두 성공할 때만 선언한다.

```bash
uv sync --all-groups
uv run python -m yakseopdong smoke
uv run pytest
```

`make pseudobulk`는 공식 expt3 zip의 integer raw count를 희소 상태에서 세포주별 합산하고, `data/processed/`에 vector-valued Parquet을 생성한다. 전체 single-cell 행렬을 dense로 변환하지 않는다.

현재 core cohort는 공식 Figshare experiment 3 원자료로 분리한다. 주 분석은 시간 일치 `DMSO_24hr_expt3` 대 `Trametinib_24hr_expt3`의 94개 eligible line이며, 저자 분석 재현용 pooled-control(`DMSO_6hr + DMSO_24hr`) 97개 line은 별도 민감도 분석으로 유지한다.

현재 QC에서 DMSO 6h–24h 중앙 PCC는 0.9727이지만, control 시간/source RMSE는 trametinib 반응 RMSE의 중앙 0.8879배다. 따라서 높은 상관만으로 pooling을 주 분석에 사용하지 않고 24시간 DMSO를 유지한다.

## 현재 모델과 W6 결론

B4는 각 fold의 training control에서만 분산 상위 5,000 genes를 고르고 whitened PCA를 fit한 뒤, PCA score에서 32,738-gene response로 multi-output ridge를 학습한다. PCA 차원과 ridge alpha는 outer-train 내부 4-fold CV로 선택한다. 외부 sensitivity와 mutation은 해석에만 쓰며 predictor에는 넣지 않는다.

```text
DMSO 24h expression → train-only gene filter/PCA → ridge → predicted trametinib Δ
```

5-fold macro 결과에서 B1 RMSE는 `0.323750`, B4 RMSE는 `0.322035`, B4 PCC-context는 `0.107666`이다. B4의 B1 대비 RMSE gain 95% paired-bootstrap CI는 `0.001168–0.002284`다. 통계적으로 방향은 일관되지만 절대·상대 개선이 작으므로 “강한 개인화 예측”이 아니라 다음 저차원 모델을 시험할 근거로 해석한다.

CCLR은 training response를 PCA 프로그램으로 압축하고 control PCA에서 각 프로그램 점수를 예측한다.

```text
DMSO 24h → train-only control PCA → ridge → response-PC scores → reconstructed trametinib Δ
```

CCLR RMSE는 `0.322383`, PCC-context는 `0.096521`이다. B1 대비 RMSE gain은 `0.001367` (95% CI `0.000823–0.001933`)이지만, B4 대비 paired gain은 `-0.000348` (95% CI `-0.000584–-0.000105`)이다. 따라서 W6에서는 저차원 response basis가 B4보다 낫다는 가설을 지지하지 않는다. 모든 fold가 response rank 20과 alpha 100을 선택했다는 경계값 현상은 W7 ablation에서 확인하되 W6 결과를 보고 grid를 바꾸지는 않는다.
