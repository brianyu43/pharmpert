# 최소 재현 체크리스트

작성일: 2026-08-27

## 완료

- [x] 공식 experiment 3 원 Matrix Market이 32,738 genes × cells의 non-negative integer count임을 전체 검증
- [x] 세 archive의 ordered `genes.tsv` 일치 확인
- [x] `cell_quality == normal` 기준 94-line time-matched 주 코호트와 97-line pooled 민감도 코호트 재현
- [x] strict 94-line `pseudobulk_24h.parquet`과 `response_24h.parquet` 생성
- [x] DMSO 6h–24h control 시간/source 차이 정량화
- [x] EGR1, ETV4/5, DUSP4/5/6, SPRY2/4의 평균 억제 방향 확인
- [x] 2,000회 cell-line bootstrap으로 marker 평균의 95% 구간 계산
- [x] cell-count, control-time, marker-response 그림 생성 및 시각 검수
- [x] 저자 `all_CL_features.rds`를 94 lines에 100% 연결하고 source checksum 검증
- [x] response PCA PC1과 외부 trametinib sensitivity 연관 정량화
- [x] lineage-aware outer 5-fold와 nested inner 4-fold 동결
- [x] B0–B4를 모든 outer fold에서 실행하고 470개 예측 저장
- [x] paired cell-line bootstrap으로 B1 대비 gain 95% CI 계산
- [x] 동일 seed/config 재실행에서 baseline prediction SHA-256 일치 확인
- [x] CCLR response basis와 ridge를 모든 nested outer fold에서 training-only로 학습
- [x] CCLR 94개 held-out 예측, component loading/score, fold artifact 생성
- [x] CCLR과 B1/B4의 paired cell-line bootstrap 비교
- [x] 동일 seed/config 재실행에서 CCLR prediction SHA-256 일치 확인
- [x] 저장한 fold artifact만으로 CCLR 예측 재구성 검증

## 관찰된 최소 재현 결과

- DMSO 6h–24h 중앙 Pearson correlation: `0.9727` (97 lines, 16,843 descriptive QC genes)
- strict 94-line 중앙 control 시간/source RMSE: `0.4054`
- strict 94-line 중앙 trametinib response RMSE: `0.4536`
- 중앙 control/treatment RMSE 비율: `0.8879`
- EGR1 평균 response: `-2.6698`, 음수 비율 `100%`
- DUSP6 평균 response: `-2.3944`, 음수 비율 `100%`
- response PC1 대 sensitivity: Pearson `-0.5989`, Spearman `-0.6683` (PCA 부호는 임의)
- B1 global mean RMSE: `0.323750`
- B4 direct ridge RMSE: `0.322035`; PCC-context: `0.107666`
- B4 RMSE gain vs B1: `0.001715` (95% CI `0.001168–0.002284`, 상대 약 `0.53%`)
- CCLR RMSE: `0.322383`; PCC-context: `0.096521`
- CCLR RMSE gain vs B1: `0.001367` (95% CI `0.000823–0.001933`)
- CCLR RMSE gain vs B4: `-0.000348` (95% CI `-0.000584–-0.000105`)
- CCLR은 94개 중 33개 line에서 B4보다 낮은 RMSE였고 macro-average에서는 B4보다 나빴다.
- 모든 outer fold가 response rank `20`, ridge alpha `100`을 선택했고 control dimension은 `5–30`이었다.

높은 상관만으로 DMSO pooling을 정당화할 수 없다. 시간/source 차이의 RMSE가 약물 반응 RMSE에 비해 작지 않으므로 주 분석은 DMSO 24h time-matched control을 유지한다. pooled control은 원 논문 재현 및 민감도 분석에만 사용한다.

## 남은 항목

- [ ] 후기 E2F/G2M/cell-cycle 프로그램 재현
- [ ] 3–48시간 time-course pseudobulk 및 시점별 최소 세포 기준 감사
- [ ] W7 response-rank/control-dimension/pathway-panel ablation
- [ ] W7 component pathway enrichment와 fold 안정성 해석

## 해석 제한

DMSO 6h와 24h는 별도 archive에서 측정되었다. 따라서 이 비교는 시간, sequencing depth, pool/batch 차이가 섞인 관찰적 QC이며 순수한 시간 인과효과가 아니다. pseudobulk는 독립 세포 집단의 요약이므로 동일 세포의 전후 trajectory로 해석하지 않는다.
