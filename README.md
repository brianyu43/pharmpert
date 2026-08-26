# 약섭동 (Yak-seopdong)

암세포주의 vehicle-control(DMSO) 전사 상태로부터 보지 못한 세포주의 24시간 trametinib 전사 반응을 예측하는 재현 가능한 분석 프로젝트다.

## 현재 상태

- Scope Lock: v1.1
- 단계: Stage 1 — 원자료 metadata 감사
- 핵심 일반화 단위: cell line
- 핵심 비교: B1 global mean response 대비 context 정보의 추가 이득
- GPU: 핵심 분석에는 불필요

이 저장소는 아직 분석 결과를 주장하지 않는다. 데이터 접근, metadata, dose, 코호트 수, split 가능성을 실제 원자료로 검증한 뒤 Stage 1로 넘어간다.

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
```

## 문서

- `docs/project_blueprint.md`: 전체 연구 청사진
- `docs/scope.md`: 포함·제외 범위와 연구 질문
- `docs/evaluation_protocol.md`: 분할, 지표, 누수 방지 규칙
- `docs/data_dictionary.md`: Stage 1에서 채울 데이터 계약
- `docs/leakage_audit.md`: fold별 누수 감사표

## 검증 게이트

Stage 0 완료는 다음 명령이 모두 성공할 때만 선언한다.

```bash
uv sync --all-groups
uv run python -m yakseopdong smoke
uv run pytest
```

Stage 1 진입 전에는 `data-probe`가 실제 AnnData 구조를 기록하고 raw count, condition, time, dose, pool/channel 필드를 확인해야 한다.

현재 core cohort는 공식 Figshare experiment 3 원자료로 분리한다. 주 분석은 시간 일치 `DMSO_24hr_expt3` 대 `Trametinib_24hr_expt3`의 94개 eligible line이며, 저자 분석 재현용 pooled-control(`DMSO_6hr + DMSO_24hr`) 97개 line은 별도 민감도 분석으로 유지한다.
