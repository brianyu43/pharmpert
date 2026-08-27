# Limitations and Claim Boundaries

이 문서는 최종 결과에서 허용되는 주장과 금지되는 과대해석을 고정한다.

1. **효과 크기가 작다.** B4의 B1 대비 RMSE 개선은 `0.001715` (95% CI `0.001168–0.002284`)로 방향은 일관되지만 상대 개선은 약 `0.53%`다. 이는 baseline context 정보의 검출을 지지하지만 강한 개인화 예측을 뜻하지 않는다.
2. **저차원 주 모델의 우위는 없었다.** CCLR은 B1보다 `0.001367` 개선했지만 B4보다 `0.000348` 나빴다. response basis를 명시적으로 압축하는 것이 direct ridge보다 낫다는 가설은 지지되지 않았다.
3. **세포 표본 수에 민감하다.** 조건당 20개 cell로 맞춘 5회 subsampling에서 B4 gain은 5/5 모두 음수였다. 원래 결과의 일부는 세포 수와 pseudobulk 측정 정밀도에 의존한다.
4. **측정 잡음이 절대 오차의 큰 부분이다.** split-half에서 추정한 full-target noise floor는 약 `0.2790` (95% CI `0.2716–0.2865`)이고 split-half PCC는 `0.2364`였다. B4 RMSE `0.3220`을 순수 모델 실패로만 해석할 수 없다.
5. **새로운 lineage 일반화는 확정되지 않았다.** leave-one-lineage-out gain은 `0.000514` (95% CI `-0.000132–0.001178`)로 0을 포함한다.
6. **시간축 외부 평가는 작고 source가 다르다.** 기술적 시간축은 24 lines 중 threshold 10을 만족한 22 lines, 모델 전이는 core와 겹치지 않는 17 lines에 제한된다. core experiment와 time-course source/batch가 달라 순수한 시간 인과효과로 분리할 수 없다.
7. **후기에 이질성이 더 커진다는 가설은 지지되지 않았다.** early 대비 late heterogeneity difference는 `-0.009666`이고 95% CI가 `-0.0342–0.0121`로 0을 포함했다.
8. **W11은 위치 이동 검사다.** 모든 control cell에 한 세포주별 latent shift를 동일하게 더한다. 따라서 covariance나 분포 모양은 그대로이며 cell-specific state, fate, paired trajectory를 예측하지 않는다.
9. **생물학 연관은 관찰적이다.** sensitivity, lineage, BRAF/KRAS와의 연관은 모델 predictor가 아니라 사후 해석 변수다. 독립 biological replicate가 제한되어 임상 효능이나 인과 메커니즘을 주장할 수 없다.
10. **RMSE의 단위를 잊으면 안 된다.** 전체 32,738-gene `Δ log1p(CPM)`의 gene-wise 평균 제곱 오차다. 특정 marker, pathway, viability, 임상 반응의 직접 오차가 아니다.

따라서 최종 허용 문장은 다음과 같다.

> Leakage-safe held-out-cell-line 평가에서 baseline 전사 상태는 global mean을 넘어서는 작고 해석 가능한 신호를 제공했지만, 그 절대 개선은 sampling-sensitive하며 강한 개인화 예측을 입증하기에 충분하지 않았다.
