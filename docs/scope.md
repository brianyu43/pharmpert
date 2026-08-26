# Scope Lock v1.1

## 연구 질문

보지 못한 암세포주의 vehicle-control(DMSO) 전사 상태는 그 세포주의 24시간 trametinib 전사 반응을 얼마나 예측하는가?

입력은 동일 세포의 처리 전 측정값이 아니라 동시에 측정한 별도 control 세포들의 요약이다. 따라서 결과를 paired trajectory, 실제 동일 세포의 전후 변화, 환자 반응 또는 임상 효능으로 해석하지 않는다.

## 핵심 estimand

세포주 `l`과 유전자 `g`에 대해 다음을 예측한다.

```text
delta[l, g] = logCPM(trametinib, 24h)[l, g] - logCPM(control, 24h)[l, g]
```

주 분석의 control은 시간 일치 `DMSO_24hr_expt3`다. 저자 코드가 사용한 `DMSO_6hr + DMSO_24hr` pooled control은 원 논문 재현 및 민감도 분석으로 분리한다. 정규화 방식과 gene universe는 Stage 1 데이터 감사 후 모델 결과를 보기 전에 고정한다.

## 포함 범위

- Trametinib 한 약물
- 24시간 broad panel의 cell-line pseudobulk 예측
- 보지 못한 cell line 단위 outer cross-validation
- B0–B4와 context-conditioned low-rank response(CCLR) 비교
- 3–48시간 time-course의 해석적 확장
- cell-line bootstrap, cell subsampling, 누수 감사
- 결과가 null이어도 포함하는 기술 보고서

## 제외 범위

- 다중 약물·dose-response 모델
- 환자 또는 임상 효능 예측
- 동일 세포의 실제 전후 trajectory 주장
- foundation model, diffusion model, 대규모 deep generative benchmark
- wet-lab 검증
- core benchmark 완료 전 mutation/lineage 입력 확장

## 성공의 최소 단위

1. 원자료와 metadata를 감사한다.
2. 누수 없는 held-out cell-line split을 고정한다.
3. 강한 기준선 B0–B4를 재현 가능하게 실행한다.
4. CCLR의 B1 대비 추가 이득 또는 null result를 정량화한다.
5. 모든 숫자를 명령형 pipeline에서 재생성한다.

설계 문서나 notebook이 존재하는 것만으로는 완료가 아니다.
