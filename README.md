# 약섭동 (Yak-seopdong)

암세포주의 vehicle-control(DMSO) 전사 상태로부터 보지 못한 세포주의 24시간 trametinib 전사 반응을 예측하는 재현 가능한 분석 프로젝트다.

## 현재 상태

- Scope Lock: v1.1
- 단계: Stage 1 완료; Stage 2 최소 재현 및 Stage 3 24h 반응 행렬 완료
- 핵심 일반화 단위: cell line
- 핵심 비교: B1 global mean response 대비 context 정보의 추가 이득
- GPU: 핵심 분석에는 불필요

공식 experiment 3 원자료에서 strict 94-line 24시간 pseudobulk와 반응 행렬을 생성했다. EGR1/DUSP6 등 사전 지정 MAPK marker의 억제 방향을 재현했지만, 아직 held-out cell-line 모델 성능이나 임상적 결과는 주장하지 않는다.

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
