# 최소 재현 체크리스트

작성일: 2026-08-26

## 완료

- [x] 공식 experiment 3 원 Matrix Market이 32,738 genes × cells의 non-negative integer count임을 전체 검증
- [x] 세 archive의 ordered `genes.tsv` 일치 확인
- [x] `cell_quality == normal` 기준 94-line time-matched 주 코호트와 97-line pooled 민감도 코호트 재현
- [x] strict 94-line `pseudobulk_24h.parquet`과 `response_24h.parquet` 생성
- [x] DMSO 6h–24h control 시간/source 차이 정량화
- [x] EGR1, ETV4/5, DUSP4/5/6, SPRY2/4의 평균 억제 방향 확인
- [x] 2,000회 cell-line bootstrap으로 marker 평균의 95% 구간 계산
- [x] cell-count, control-time, marker-response 그림 생성 및 시각 검수

## 관찰된 최소 재현 결과

- DMSO 6h–24h 중앙 Pearson correlation: `0.9727` (97 lines, 16,843 descriptive QC genes)
- strict 94-line 중앙 control 시간/source RMSE: `0.4054`
- strict 94-line 중앙 trametinib response RMSE: `0.4536`
- 중앙 control/treatment RMSE 비율: `0.8879`
- EGR1 평균 response: `-2.6698`, 음수 비율 `100%`
- DUSP6 평균 response: `-2.3944`, 음수 비율 `100%`

높은 상관만으로 DMSO pooling을 정당화할 수 없다. 시간/source 차이의 RMSE가 약물 반응 RMSE에 비해 작지 않으므로 주 분석은 DMSO 24h time-matched control을 유지한다. pooled control은 원 논문 재현 및 민감도 분석에만 사용한다.

## 남은 항목

- [ ] trametinib sensitivity와 반응 PCA 축의 연관 방향
- [ ] 후기 E2F/G2M/cell-cycle 프로그램 재현
- [ ] 3–48시간 time-course pseudobulk 및 시점별 최소 세포 기준 감사
- [ ] control 상태 PCA와 cell-line identity 설명력 정량화

## 해석 제한

DMSO 6h와 24h는 별도 archive에서 측정되었다. 따라서 이 비교는 시간, sequencing depth, pool/batch 차이가 섞인 관찰적 QC이며 순수한 시간 인과효과가 아니다. pseudobulk는 독립 세포 집단의 요약이므로 동일 세포의 전후 trajectory로 해석하지 않는다.
