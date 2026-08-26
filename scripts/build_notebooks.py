"""Build the reader-facing Stage 1/2 notebooks with nbformat."""

from __future__ import annotations

from pathlib import Path

import nbformat as nbf


def _notebook(cells: list[nbf.NotebookNode]) -> nbf.NotebookNode:
    notebook = nbf.v4.new_notebook(cells=cells)
    notebook.metadata.kernelspec = {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    }
    notebook.metadata.language_info = {"name": "python", "version": "3.13"}
    return notebook


ROOT_CELL = """from pathlib import Path
import json
import pandas as pd
from IPython.display import Image, display

root = next(
    path for path in [Path.cwd(), *Path.cwd().parents]
    if (path / \"pyproject.toml\").exists() and (path / \"cell_count_matrix.csv\").exists()
)
"""


def build_data_audit(root: Path) -> None:
    cells = [
        nbf.v4.new_markdown_cell(
            """# Experiment 3 data audit

## tl;dr

공식 Figshare experiment 3 원자료에서 시간 일치 24시간 주 코호트 94개와
pooled-control 민감도 코호트 97개를 확인했다. 이 노트북은 원자료 감사 로그와
추적 가능한 표를 다시 읽어 핵심 불변조건을 검증한다."""
        ),
        nbf.v4.new_code_cell(ROOT_CELL),
        nbf.v4.new_code_cell(
            """core = json.loads((root / \"results/logs/core_audit_summary.json\").read_text())
counts = pd.read_csv(root / \"cell_count_matrix.csv\")

audit_summary = pd.DataFrame({
    \"metric\": [
        \"genes\", \"normal cells: DMSO 6h\", \"normal cells: DMSO 24h\",
        \"normal cells: trametinib 24h\", \"strict primary lines\", \"pooled sensitivity lines\"
    ],
    \"value\": [
        32_738,
        core[\"archives\"][\"dmso_6h\"][\"cells_normal\"],
        core[\"archives\"][\"dmso_24h\"][\"cells_normal\"],
        core[\"archives\"][\"trametinib_24h\"][\"cells_normal\"],
        core[\"primary_strict_eligible_lines\"],
        core[\"reproduction_pooled_eligible_lines\"],
    ],
})
display(audit_summary)
assert core[\"primary_strict_eligible_lines\"] == 94
assert core[\"reproduction_pooled_eligible_lines\"] == 97
assert counts[\"cell_line\"].is_unique
"""
        ),
        nbf.v4.new_markdown_cell(
            """## Context & Methods

### Key Assumptions

- 분석 단위는 세포가 아니라 cell line이다.
- `cell_quality == "normal"`만 포함한다.
- 주 분석은 `DMSO_24hr_expt3` 대 `Trametinib_24hr_expt3`다.
- 각 조건에서 정상 세포 20개 이상인 세포주만 주 코호트에 포함한다.
- DMSO 6h+24h pooling은 주 분석이 아니라 민감도 분석이다."""
        ),
        nbf.v4.new_markdown_cell("## Data"),
        nbf.v4.new_code_cell(
            """manifest = pd.read_csv(root / \"data_manifest.csv\")
display(manifest[[\"source_role\", \"bytes\", \"sha256\", \"status\"]])
display(counts.head(10))
"""
        ),
        nbf.v4.new_markdown_cell("## Results"),
        nbf.v4.new_code_cell(
            """condition_totals = counts[[
    \"dmso_6h_normal\", \"dmso_24h_normal\", \"trametinib_24h_normal\"
]].sum().rename(\"normal_cells\").to_frame()
display(condition_totals)
display(Image(filename=str(root / \"results/figures/cell_count_heatmap.png\")))
"""
        ),
        nbf.v4.new_markdown_cell(
            """## Takeaways

- 단순 통합 AnnData 필터 대신 experiment 3 원자료를 직접 사용해야 한다.
- 24시간 time-matched 주 코호트는 94개 세포주다.
- 저자 방식의 pooled-control 재현 코호트는 97개지만 별도 민감도 분석으로 유지한다.
- 원 Matrix Market 행렬은 32,738 genes × cells의 non-negative integer raw count다."""
        ),
    ]
    nbf.write(_notebook(cells), root / "notebooks" / "00_data_audit.ipynb")


def build_qc(root: Path) -> None:
    cells = [
        nbf.v4.new_markdown_cell(
            """# Pseudobulk QC and minimal biological reproduction

## tl;dr

94개 주 코호트의 `log1p(CPM)` 반응 행렬을 생성했다. DMSO 6h–24h 상관은
높지만 시간/원자료 차이의 크기가 약물 반응과 비슷하므로 pooled control은
민감도 분석에만 사용한다. 사전 지정한 8개 immediate-early MAPK marker는
평균적으로 모두 억제 방향이다."""
        ),
        nbf.v4.new_code_cell(ROOT_CELL),
        nbf.v4.new_code_cell(
            """summary = json.loads((root / \"results/logs/qc_summary.json\").read_text())
control_qc = pd.read_csv(root / \"results/tables/control_time_qc.csv\")
markers = pd.read_csv(root / \"results/tables/marker_response_summary.csv\")

key_metrics = pd.DataFrame({
    \"metric\": [
        \"QC genes\", \"DMSO 6h–24h median PCC\", \"DMSO time/source median RMSE\",
        \"Trametinib response median RMSE\", \"control/treatment RMSE ratio\"
    ],
    \"value\": [
        summary[\"qc_gene_rule\"][\"gene_count\"],
        summary[\"control_time_all_97\"][\"median_pcc\"],
        summary[\"strict_94_comparison\"][\"median_control_time_rmse\"],
        summary[\"strict_94_comparison\"][\"median_trametinib_response_rmse\"],
        summary[\"strict_94_comparison\"][\"median_control_to_treatment_rmse_ratio\"],
    ],
})
display(key_metrics)
"""
        ),
        nbf.v4.new_markdown_cell(
            """## Context & Methods

### Key Assumptions

- pseudobulk는 cell line × condition별 raw UMI 합이다.
- 정규화는 `log1p(1e6 × count / library_size)`로 고정한다.
- DMSO 시간/원자료 QC 유전자는 두 control에서 총 count ≥100이고 각 시점
  10개 이상 세포주에서 발현된 16,843개다.
- 이 QC gene mask는 기술적 비교용이며 모델 feature selection에는 사용하지 않는다.
- marker 평균의 구간은 seed 20260827, 2,000회 cell-line bootstrap이다."""
        ),
        nbf.v4.new_markdown_cell("## Data"),
        nbf.v4.new_code_cell(
            """from yakseopdong.pseudobulk import read_vector_parquet

primary_meta, primary_matrix = read_vector_parquet(
    root / \"data/processed/pseudobulk_24h.parquet\", \"log1p_cpm\"
)
response_meta, response_matrix = read_vector_parquet(
    root / \"data/processed/response_24h.parquet\", \"delta_log1p_cpm\"
)
print(\"pseudobulk_24h\", primary_meta.shape, primary_matrix.shape)
print(\"response_24h\", response_meta.shape, response_matrix.shape)
display(pd.read_csv(root / \"processed_manifest.csv\"))
assert primary_matrix.shape == (188, 32_738)
assert response_matrix.shape == (94, 32_738)
"""
        ),
        nbf.v4.new_markdown_cell("## Results"),
        nbf.v4.new_code_cell(
            """display(markers)
assert (markers[\"mean_delta_log1p_cpm\"] < 0).all()
assert summary[\"pooling_used_for_primary\"] is False
display(Image(filename=str(root / \"results/figures/control_time_qc.png\")))
display(Image(filename=str(root / \"results/figures/marker_response.png\")))
"""
        ),
        nbf.v4.new_markdown_cell(
            """## Takeaways

- DMSO 6h–24h 중앙 PCC는 약 0.973으로 전체 발현 순서는 안정적이다.
- 그러나 strict 94-line에서 control 시간/원자료 RMSE는 trametinib 반응
  RMSE의 중앙 약 0.888이다. 이는 pooling 효과를 무시하기 어렵다는 뜻이다.
- EGR1과 DUSP6는 94개 세포주 모두 억제 방향이며, 사전 지정한 8개 marker의 평균 반응도 모두 음수다.
- 따라서 주 target은 DMSO 24h time-matched 반응으로 동결하고
  pooled-control은 민감도 분석으로만 유지한다.
- 이 비교는 시간과 archive/batch 차이가 섞인 관찰적 QC이므로 순수한 시간
  인과효과로 해석하지 않는다."""
        ),
    ]
    nbf.write(_notebook(cells), root / "notebooks" / "01_qc.ipynb")


def build_response_landscape(root: Path) -> None:
    cells = [
        nbf.v4.new_markdown_cell(
            """# Control and trametinib response landscape

## tl;dr

94개 세포주의 DMSO 24h control과 trametinib response를 각각 PCA로 요약했다.
Response PC1은 저자 제공 외부 민감도와 Pearson `-0.599`, Spearman `-0.668`로
연관된다. PCA 축의 부호는 임의이므로 부호보다 연관 강도를 본다. 이 분석은 전체
코호트 탐색용이며 모델 선택에는 사용하지 않는다."""
        ),
        nbf.v4.new_code_cell(ROOT_CELL),
        nbf.v4.new_code_cell(
            """summary = json.loads((root / "results/logs/landscape_summary.json").read_text())
landscape = pd.read_csv(root / "results/tables/response_landscape.csv")
directions = pd.read_csv(root / "results/tables/direction_comparison.csv")
display(pd.DataFrame({
    "metric": [
        "cell lines", "control PC1 variance", "response PC1 variance",
        "response PC1 vs sensitivity Pearson", "response PC1 vs sensitivity Spearman",
        "median control-vs-drug cosine",
    ],
    "value": [
        summary["cell_lines"], summary["control_pc_explained_variance_ratio"][0],
        summary["response_pc_explained_variance_ratio"][0],
        summary["response_pc1_sensitivity_pearson"],
        summary["response_pc1_sensitivity_spearman"],
        summary["median_control_vs_drug_cosine"],
    ],
}))
assert len(landscape) == len(directions) == 94
"""
        ),
        nbf.v4.new_markdown_cell(
            """## Context & Methods

### Key Assumptions

- 저자 제공 `all_CL_features.rds`의 Disease, sensitivity, mutation 필드를 DepMap ID로 연결했다.
- sensitivity와 mutation은 해석용이며 B0–B4 예측 입력이 아니다.
- PCA는 각 행렬에서 분산이 큰 5,000 genes를 사용한 전체 코호트 탐색 그림이다.
- DMSO 6→24h 변화와 drug response는 DMSO 24h를 반대 부호로 공유하므로 방향 비교는
  독립적인 인과 대조가 아니라 진단적 비교다."""
        ),
        nbf.v4.new_markdown_cell("## Results"),
        nbf.v4.new_code_cell(
            """display(Image(filename=str(root / "results/figures/response_landscape.png")))
display(Image(filename=str(root / "results/figures/response_pc1_sensitivity.png")))
display(Image(filename=str(root / "results/figures/direction_comparison.png")))
"""
        ),
        nbf.v4.new_markdown_cell(
            """## Takeaways

- baseline control PC1은 변동의 12.9%, response PC1은 5.9%를 설명한다.
- response의 주요 축은 외부 약물 민감도와 연결되므로 섭동 신호가 생물학적 차이를 담는다.
- control 시간/source 변화와 약물 반응은 크기가 비슷해도 방향은 같지 않다.
- 다만 공유 DMSO 24h anchor가 음의 유사도를 유도할 수 있어 방향 수치를
  독립 효과로 해석하지 않는다."""
        ),
    ]
    nbf.write(_notebook(cells), root / "notebooks" / "02_response_landscape.ipynb")


def build_baselines(root: Path) -> None:
    cells = [
        nbf.v4.new_markdown_cell(
            """# Held-out cell-line baselines B0–B4

## tl;dr

94개 세포주를 lineage-aware 5-fold로 나누고, 각 outer train 안에서 4-fold로
hyperparameter를 선택했다. B4 direct ridge가 B1 global mean보다 RMSE를 평균
`0.001715` 개선했으며 95% paired cell-line bootstrap CI는 `0.001168–0.002284`다.
방향은 일관되지만 B1 대비 상대 개선은 약 0.53%로 작다."""
        ),
        nbf.v4.new_code_cell(ROOT_CELL),
        nbf.v4.new_code_cell(
            """summary = json.loads((root / "results/logs/baseline_summary.json").read_text())
comparison = pd.read_csv(root / "results/tables/baseline_comparison.csv")
hyperparameters = pd.read_csv(root / "results/tables/baseline_hyperparameters.csv")
display(comparison[[
    "model", "rmse_delta_mean", "pcc_delta_mean", "pcc_context_mean",
    "rmse_gain_vs_b1_mean", "rmse_gain_vs_b1_ci95_low", "rmse_gain_vs_b1_ci95_high",
]])
display(hyperparameters)
assert summary["each_cell_line_predicted_once_per_model"]
assert summary["predictions_rows"] == 470
"""
        ),
        nbf.v4.new_markdown_cell(
            """## Context & Methods

### Model definitions

- B0: no change, `delta=0`
- B1: outer-training lines의 평균 반응
- B2: 같은 lineage의 training 평균; 두 줄 미만이면 B1 fallback
- B3: training-control-only PCA에서 최근접 세포주의 반응
- B4: training-control-only PCA score에서 32,738-gene response로 가는 multi-output ridge

### Leakage controls

- outer test 단위는 cell line이며 각 line은 정확히 한 번 test다.
- top-5,000 control gene 선택과 PCA는 각 outer/inner training partition에서 다시 fit한다.
- B3 차원과 B4 차원/alpha는 outer-train 내부 4-fold RMSE로만 고른다.
- sensitivity, mutation, outer-test treated response는 predictor fit에 사용하지 않는다."""
        ),
        nbf.v4.new_markdown_cell("## Results"),
        nbf.v4.new_code_cell(
            """display(Image(filename=str(root / "results/figures/baseline_performance.png")))
b1 = comparison.set_index("model").loc["B1"]
b4 = comparison.set_index("model").loc["B4"]
relative_gain = b4["rmse_gain_vs_b1_mean"] / b1["rmse_delta_mean"]
print(f"B4 relative RMSE improvement vs B1: {relative_gain:.2%}")
"""
        ),
        nbf.v4.new_markdown_cell(
            """## Takeaways

- B1이 B0보다 낫다: 모든 세포주에 공유되는 trametinib 반응이 크다.
- B2와 B3는 B1보다 나쁘다: lineage 평균이나 단일 최근접 이웃은 반응
  이질성을 안정적으로 설명하지 못한다.
- B4는 B1을 근소하게 이기고 PCC-context는 `0.108`이다. baseline 전사 맥락의 추가 정보는
  검출되지만 현재 저복잡도 모델에서 실용적 크기는 작다.
- 다음 모델은 B4보다 복잡하게 만들기 전에, 평균 반응에 저차원 response component만 더하는
  CCLR이 이 작은 이득을 재현·확대하는지 같은 nested CV에서 확인해야 한다."""
        ),
    ]
    nbf.write(_notebook(cells), root / "notebooks" / "03_baselines.ipynb")


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    (root / "notebooks").mkdir(exist_ok=True)
    build_data_audit(root)
    build_qc(root)
    build_response_landscape(root)
    build_baselines(root)


if __name__ == "__main__":
    main()
