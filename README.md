# 약섭동 (Yak-seopdong)

암세포주의 동시에 측정된 vehicle-control(DMSO) 전사 상태로부터, 보지 못한 세포주의 24시간 trametinib 전사 반응을 예측하고 그 한계를 검증한 재현 가능한 분석 프로젝트다.

## 최종 상태

- Stage 0–12 분석·동결 완료, W13–W15 보고·재현·릴리스 검증 수행
- 주 일반화 단위: `cell_line`; lineage-aware outer 5-fold / nested inner 4-fold
- 주 코호트: 94 cell lines, 32,738 genes, 16,588 normal cells
- 주 모델: B4 direct multi-output ridge; 주 비교군: B1 outer-training mean response
- 최종 예측: B0–B4와 CCLR, 총 564 held-out line-model rows
- GPU: 필요 없음; 모든 결과는 CPU에서 생성

주 결과는 명확하지만 크기는 작다. B1 RMSE는 `0.323750`, B4 RMSE는 `0.322035`이며 paired improvement는 `0.001715` (95% CI `0.001168–0.002284`, 약 `0.53%`)이다. CCLR은 B1보다 낫지만 B4보다 `0.000348` 나빴다. inclusion threshold와 gene filter에는 방향이 유지됐지만, 조건당 20개 cell로 균등 subsampling하면 5/5 반복에서 B4 gain이 음수가 됐다. 따라서 결론은 **baseline context의 추가 정보는 검출되지만, 강한 개인화 예측 성능은 입증되지 않았다**이다.

## 빠른 시작

필수 도구는 Python 3.13과 [uv](https://docs.astral.sh/uv/)다.

```bash
make install
make smoke
make test
```

원자료 확보 후 전체 파이프라인은 다음 순서로 실행한다.

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
make ablation
make validate-ablation
make temporal
make validate-temporal
make biology
make validate-biology
make robustness
make validate-robustness
make distribution
make validate-distribution
make freeze-release
make validate-release
make notebooks
```

`make pseudobulk`는 공식 experiment 3의 integer raw count를 희소 상태에서 세포주별 합산한다. 전체 single-cell 행렬을 dense로 변환하지 않는다. 원자료와 중간 processed data는 Git에 넣지 않으며 checksum과 manifest로 식별한다.

## 모델과 평가

```text
DMSO 24h expression
  → outer-training-only variable-gene filter/PCA
  → nested-CV multi-output ridge
  → predicted trametinib 24h − DMSO 24h
```

B4는 각 training partition에서만 분산 상위 5,000 genes를 선택하고 whitened PCA를 fit한 뒤, control PC score에서 전체 32,738-gene response로 ridge regression을 학습한다. 외부 sensitivity와 mutation은 해석에만 쓰며 predictor에는 넣지 않는다. CCLR은 response를 training-only PCA로 한 번 더 압축하지만, 이 추가 low-rank 제약은 B4를 개선하지 못했다.

시간축 17-line 외부 전이에서 B4는 3h/6h에는 B1보다 나쁘고, 12h에는 차이가 불확실하며, 24h와 48h에만 작게 나았다. W11 single-cell 분포 분석은 17 lines, 1,892 cells에서 통과했지만 모든 control cell에 동일한 latent shift를 더하는 검사이므로 분포 모양이나 개별 세포 trajectory를 예측한 결과가 아니다.

## 주요 산출물

- `report/final_report.md`: 최종 기술 보고서
- `report/supplementary.md`: 보충 방법·결과·재현 정보
- `report/limitations.md`: 주장 경계를 정한 제한점 목록
- `results/final_predictions.parquet`: B0–B4/CCLR held-out 예측
- `results/final_metrics.csv`: 핵심·시간축·강건성·분포 지표
- `results/final_cell_lines.csv`: 94-line annotation과 오류 지표
- `results/figure_manifest.csv`, `results/table_manifest.csv`: 번호·경로·checksum manifest
- `release/final_config.json`: 동결된 release contract
- `results/logs/release_validation.json`: 독립 릴리스 검증 결과

세부 데이터 계약과 누수 방지 규칙은 `docs/data_dictionary.md`, `docs/evaluation_protocol.md`, `docs/leakage_audit.md`를 따른다. 전체 완료 게이트는 `docs/reproduction_checklist.md`에 기록한다.

## 인용

원자료와 생물학적 맥락은 McFarland et al., *Nature Communications* 11, 4296 (2020), DOI `10.1038/s41467-020-17440-w`를 따른다. 이 저장소의 결과는 독립 반복이 제한된 관찰적 예측 분석이며 임상 효능이나 인과효과를 주장하지 않는다.
