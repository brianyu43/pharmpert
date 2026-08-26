# Data Dictionary

상태: Stage 1–5 원자료, pseudobulk, annotation, split, baseline 계약 확인 완료. 아래 항목은 추정값이 아니라 실행으로 확인한 값만 기록한다.

## Source manifest

| 항목 | 확인값 |
|---|---|
| Primary URL | `https://exampledata.scverse.org/pertpy/mcfarland_2020.h5ad` |
| Download timestamp | `2026-08-26T11:34:45.250989+00:00` |
| File name | `data/mcfarland_2020.h5ad` |
| File size | `1,459,410,830 bytes` |
| SHA-256 | `94a7240047b9ce6822dc2aa1e7a66c1fd8b00be842e5a95fc72cbeb6f2834d2c` |
| pertpy version | `1.3.0` |
| AnnData shape | `182,875 cells × 32,738 genes` |
| Figshare dataset license | CC BY 4.0 |

## Matrix contract

| 항목 | 확인값 |
|---|---|
| `adata.X` 의미 | pertpy 통합 파일에서는 raw UMI count 후보; 비영값 100,000개 표본이 모두 정수형 |
| raw count 위치 | 별도 layer 없음; `adata.raw` 없음; 후보는 `adata.X`뿐 |
| sparse format/dtype | CSC sparse, `float32`, 711,380,841 nonzero entries |
| gene identifier | `var_names`는 gene symbol, `var["ensembl_id"]`에 Ensembl ID |
| official expt3 raw count | 각 Figshare zip의 `matrix.mtx`; Matrix Market `integer general`, non-negative 전체 검증 |
| raw matrix orientation | 32,738 genes × cells |
| duplicated gene IDs | expt3 `genes.tsv`의 Ensembl ID 32,738개는 결측·중복 없음 |

## Processed matrix contract

| 파일 | 행 | 값 열 | 의미 |
|---|---:|---|---|
| `pseudobulk_24h.parquet` | 188 | `log1p_cpm[32738]` | strict 94 lines × control/trametinib |
| `response_24h.parquet` | 94 | `delta_log1p_cpm[32738]` | trametinib 24h − DMSO 24h |
| `pseudobulk_control_time.parquet` | 194 | `log1p_cpm[32738]` | 97 lines × DMSO 6h/24h |
| `pseudobulk_pooled_sensitivity.parquet` | 194 | `log1p_cpm[32738]` | 97 lines × pooled DMSO/trametinib |
| `response_pooled_sensitivity.parquet` | 97 | `delta_log1p_cpm[32738]` | trametinib 24h − pooled DMSO |
| `gene_metadata.parquet` | 32,738 | `gene_id`, `gene_symbol` | 모든 vector 값의 고정 순서 |

각 expression/response 값은 Arrow fixed-size list로 저장한다. pseudobulk raw UMI 합을 library size로 나눈 뒤 `log1p(CPM)`을 적용하며 scaling factor는 1,000,000이다. 처리 파일은 Git에서 제외되고 `make pseudobulk`로 재생성한다.

## Annotation and model-output contract

| 파일 | 행 | 의미 |
|---|---:|---|
| `cell_line_annotations.csv` | 94 | DepMap ID, 저자 Disease/Subtype, Trametinib AUC/sensitivity, BRAF/RAS 변이 |
| `split_assignments.csv` | 94 | 각 cell line의 유일한 outer test fold |
| `inner_split_assignments.csv` | 376 | 5개 outer train 각각의 4-fold inner assignment |
| `baseline_predictions.parquet` | 470 | 94 lines × B0–B4, `predicted_delta_log1p_cpm[32738]` |
| `baseline_metrics_by_line.csv` | 470 | cell-line 단위 RMSE/NRMSE/PCC/Spearman/top-50/context/gain |

`all_CL_features.rds`는 Figshare v3의 `Trametinib_24hr_expt3`와 `metadata` 객체를 사용한다. strict 94 lines에서 lineage, sensitivity, 네 hotspot mutation field의 필수값은 모두 완전하다. PRISM AUC 2개와 GDSC AUC 39개는 개별 source에서 결측이지만 저자의 결합 `AUC_avg`와 `sens=1-AUC_avg`는 94개 모두 존재한다.

## Observation fields

| 개념 | 실제 열 | 단위/값 | 결측/주의사항 |
|---|---|---|---|
| cell line | `cell_line` | 209 unique | 결측 없음; core 실험 식별 후 재집계 필요 |
| perturbation | `perturbation` | `Trametinib`, `control` 등 18개 | 결측 없음 |
| control | `perturbation`, time-course에서는 `hash_tag` | `control`; `DMSO_*hr`; `Untreated_48hr` | time-course의 `perturbation`만으로 조건 구분 불가 |
| time | `time`, time-course에서는 `hash_tag` | `24`, `6`, `72, 96`, `3, 6, 12, 24, 48` | `time`은 일부 행에서 실험 전체 시점 목록이므로 cell-level 시간이 아님 |
| dose | `dose_value`, `dose_unit` | trametinib `0.1 µM`, control `0.0 µM` | Supplementary Table 3과 AnnData metadata가 일치 |
| pool/channel | `channel`, `hash_assignment`, `hash_tag` | 문자열 `nan`, 숫자/문자 channel, hash condition | 99-line expt3 식별 규칙이 아직 미확정 |
| quality/singlet | `cell_quality`, `hash_tag` | `normal` 등 5개; time-course `multiplet`/`unknown` | time-course clean tag filtering 필요 |
| lineage | 저자 RDS `metadata.Disease` | 21 unique | strict 94 lines 모두 연결; split balancing과 B2에만 사용 |

## Cohort counts

- 24시간 strict time-matched core: normal cell 기준 양 조건 각각 20개 이상인 94개 line
- 저자 pooled-control 재현 cohort: DMSO 6h+24h와 trametinib 24h에서 각각 20개 이상인 97개 line
- time-course clean hashtag cell 수: 13,713개로 원 논문 수치와 일치
- time-course의 모든 필수 시점을 만족하는 cell line: 미확인(후속 시간축 단계)
- strict 24h 제외 규칙: DMSO 24h 또는 trametinib 24h 정상 세포 20개 미만
- strict 24h 제외 목록: `JHOM1_OVARY`(DMSO 24h 11 cells),
  `SNU410_PANCREAS`(18), `SNU61_LARGE_INTESTINE`(18)

## 첫 감사에서 확인한 주의사항

1. `time == "24"`와 `perturbation in {control, Trametinib}`만으로는 99-line core가 분리되지 않는다.
2. time-course row의 `time` 값은 시점 목록 전체이며 실제 조건은 `hash_tag`에서 읽어야 한다.
3. 저자 코드는 experiment 3 trametinib 분석에서 `DMSO_6hr_expt3`와 `DMSO_24hr_expt3`를 control로 함께 사용한다. 이 pooling이 97-line 기준을 재현한다.
4. DMSO 6h와 24h는 별도 archive이므로 두 조건의 차이는 순수 시간효과가 아니라 시간과 source/batch 차이가 섞인 관찰적 QC다.
