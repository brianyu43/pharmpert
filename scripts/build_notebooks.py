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


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    (root / "notebooks").mkdir(exist_ok=True)
    build_data_audit(root)
    build_qc(root)


if __name__ == "__main__":
    main()
