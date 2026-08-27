"""Build the reader-facing project notebooks with nbformat."""

from __future__ import annotations

import json
from pathlib import Path

import nbformat as nbf
import pandas as pd


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


def build_main_model(root: Path) -> None:
    report = json.loads((root / "results/logs/cclr_summary.json").read_text())
    metrics = report["primary_metrics"]
    gain_b1 = metrics["rmse_gain_vs_b1_mean"]
    gain_b4 = metrics["rmse_gain_vs_b4_mean"]
    b4_low = metrics["rmse_gain_vs_b4_ci95_low"]
    b4_high = metrics["rmse_gain_vs_b4_ci95_high"]
    comparison_phrase = (
        "CCLR이 B4보다 일관되게 낫다"
        if b4_low > 0
        else "CCLR의 B4 대비 차이는 0을 가로질러 우위를 확정하지 못한다"
    )
    cells = [
        nbf.v4.new_markdown_cell(
            f"""# Context-conditioned low-rank response (CCLR)

## tl;dr

94개 세포주의 lineage-aware outer 5-fold와 nested inner 4-fold에서 CCLR을
평가했다. CCLR의 B1 대비 RMSE gain은 `{gain_b1:.6f}`, B4 대비 gain은
`{gain_b4:.6f}` (95% CI `{b4_low:.6f}–{b4_high:.6f}`)다. {comparison_phrase}.
이 결과는 32,738-gene 반응을 training-fold의 공유 response component로 제한했을 때
직접 gene-level ridge보다 일반화가 좋아지는지를 검정한 것이다."""
        ),
        nbf.v4.new_code_cell(ROOT_CELL),
        nbf.v4.new_code_cell(
            """summary = json.loads((root / "results/logs/cclr_summary.json").read_text())
validation = json.loads((root / "results/logs/cclr_validation.json").read_text())
comparison = pd.read_csv(root / "results/tables/model_comparison_w6.csv")
hyperparameters = pd.read_csv(root / "results/tables/cclr_hyperparameters.csv")
component_summary = pd.read_csv(root / "results/tables/cclr_component_summary.csv")
loadings = pd.read_csv(root / "results/tables/cclr_component_top_loadings.csv")
artifact_manifest = pd.read_csv(root / "results/tables/cclr_model_artifacts.csv")
display(comparison.loc[comparison["model"].isin(["B1", "B4", "CCLR"]), [
    "model", "rmse_delta_mean", "pcc_delta_mean", "pcc_context_mean",
    "rmse_gain_vs_b1_mean", "rmse_gain_vs_b1_ci95_low", "rmse_gain_vs_b1_ci95_high",
    "rmse_gain_vs_b4_mean", "rmse_gain_vs_b4_ci95_low", "rmse_gain_vs_b4_ci95_high",
]])
display(hyperparameters)
assert summary["each_cell_line_predicted_once"]
assert summary["predictions_rows"] == 94
assert summary["outer_test_response_used_for_fit"] is False
assert summary["repeat_run_prediction_sha256_match"] is True
assert len(artifact_manifest) == 5
assert validation["status"] == "pass"
"""
        ),
        nbf.v4.new_markdown_cell(
            """## Context & Methods

### Model

각 outer-training fold에서 response를 `Δ ≈ μ + Ws`로 PCA 압축한다. 같은
training control에서 상위 5,000 variable genes와 whitened PCA score `z`를 만들고,
ridge로 `z → s`를 학습한 뒤 `μ + Wŝ`로 32,738-gene response를 복원한다.

### Key Assumptions

- outer test 단위는 cell line이며 각 line은 정확히 한 번만 예측한다.
- control gene filter/PCA와 response PCA는 각 inner/outer training partition에서 다시 fit한다.
- control 차원, response rank, ridge alpha는 inner-CV macro RMSE로만 선택한다.
- response PCA component의 부호와 번호는 fold별로 임의이고 서로 직접 정렬되지 않는다.
- outer-test response는 모델 fit에 쓰지 않고, 예측 완료 후 metric과 component-score 평가에만 쓴다.
- sensitivity와 mutation은 predictor가 아니다."""
        ),
        nbf.v4.new_markdown_cell("## Results"),
        nbf.v4.new_code_cell(
            """display(Image(filename=str(root / "results/figures/cclr_performance.png")))
display(Image(filename=str(root / "results/figures/cclr_components.png")))
"""
        ),
        nbf.v4.new_markdown_cell("### Component inspection"),
        nbf.v4.new_code_cell(
            """display(component_summary.head(20))
display(loadings.loc[
    loadings["loading_rank"].le(10),
    ["outer_fold", "component", "direction", "loading_rank", "gene_symbol", "loading"],
].head(40))
"""
        ),
        nbf.v4.new_markdown_cell(
            f"""## Takeaways

- CCLR의 B1 대비 평균 RMSE gain은 `{gain_b1:.6f}`다. 이는 baseline context가
  공유 저차원 반응 프로그램의 세기를 예측하는지에 대한 주 검정이다.
- CCLR의 B4 대비 paired gain 95% CI는 `{b4_low:.6f}–{b4_high:.6f}`다.
  {comparison_phrase}.
- response rank와 control 차원이 fold마다 달라질 수 있어, component loading은
  fold별 해석 단위로 유지하며 component 번호를 전체 코호트의 고정 pathway로 보지 않는다.
- W6는 모델 비교를 완료하지만 pathway enrichment, rank/feature ablation, lineage·mutation
  확장은 W7의 별도 분석이다."""
        ),
    ]
    nbf.write(_notebook(cells), root / "notebooks" / "04_main_model.ipynb")


def build_ablation(root: Path) -> None:
    comparison = pd.read_csv(root / "results/tables/ablation_metrics.csv").set_index("model")
    pathway = comparison.loc["PATHWAY_D20_R20"]
    lineage = comparison.loc["LINEAGE_D20_R20"]
    mutations = comparison.loc["BRAF_KRAS_D20_R20"]
    cells = [
        nbf.v4.new_markdown_cell(
            f"""# W7 fixed ablation and leakage audit

## tl;dr

동결된 94개 세포주 outer 5-fold에서 16개 비교군을 평가했다. 모든 고정
low-rank 변형은 B1 평균 반응보다 나았지만 B4 direct ridge를 확실히 넘지 못했다.
Pathway panel의 B4 대비 RMSE gain은 `{pathway['rmse_gain_vs_b4_mean']:.6f}`
(95% CI `{pathway['rmse_gain_vs_b4_ci95_low']:.6f}–{pathway['rmse_gain_vs_b4_ci95_high']:.6f}`),
lineage 추가는 `{lineage['rmse_gain_vs_b4_mean']:.6f}`, BRAF/KRAS 추가는
`{mutations['rmse_gain_vs_b4_mean']:.6f}`다. 따라서 W7은 외부 정보가 격차를
줄일 가능성은 보이지만 B4 우위나 새 주 모델을 확정하지 않는다."""
        ),
        nbf.v4.new_code_cell(ROOT_CELL),
        nbf.v4.new_code_cell(
            """summary = json.loads((root / "results/logs/ablation_summary.json").read_text())
comparison = pd.read_csv(root / "results/tables/ablation_metrics.csv")
variants = pd.read_csv(root / "results/tables/ablation_variants.csv")
coverage = pd.read_csv(root / "results/tables/pathway_panel_coverage.csv")
enrichment = pd.read_csv(root / "results/tables/component_pathway_enrichment.csv")
stability = pd.read_csv(root / "results/tables/cclr_subspace_stability.csv")
display(comparison[[
    "model", "rmse_delta_mean", "pcc_context_mean",
    "rmse_gain_vs_b1_mean", "rmse_gain_vs_b1_ci95_low", "rmse_gain_vs_b1_ci95_high",
    "rmse_gain_vs_b4_mean", "rmse_gain_vs_b4_ci95_low", "rmse_gain_vs_b4_ci95_high",
]])
assert summary["cell_lines"] == 94
assert summary["variants"] == 16
assert summary["outer_test_response_used_for_fit"] is False
assert summary["outer_test_used_for_variant_selection"] is False
"""
        ),
        nbf.v4.new_markdown_cell(
            """## Context & Methods

W7은 W6 결과를 소급 재튜닝하는 단계가 아니라 같은 frozen outer test에서 모델의
어떤 제약과 입력이 성능을 바꾸는지 확인하는 진단 분석이다.

### Key Assumptions

- 기준점은 control dimension 20, response rank 20, ridge alpha 100이다.
- control dimension `[5, 10, 20, 30]`과 response rank `[2, 5, 10, 20, 30, 40, 50]`을 모두 보고한다.
- rank 30/40/50은 W6의 rank 20 상한 선택 뒤 사전 동결한 진단 확장이며 W6 주 결과를 바꾸지 않는다.
- MSigDB Hallmark 2026.1.Hs 6개 source set과 기존 8개 immediate-early marker를
  결과 확인 전에 동결했다.
- pathway 입력만 제한하고 response target은 모든 32,738 genes로 유지한다.
- lineage one-hot과 BRAF/KRAS 표준화는 각 outer-training fold에서만 fit한다.
- outer-test response와 sensitivity는 predictor나 variant 선택에 사용하지 않는다.
- 모든 변형을 보고하며 outer-test 평균이 가장 좋은 변형을 새 주 모델로 승격하지 않는다."""
        ),
        nbf.v4.new_markdown_cell("## Data"),
        nbf.v4.new_code_cell(
            """annotations = pd.read_csv(root / "cell_line_annotations.csv")
display(pd.DataFrame({
    "field": ["cell lines", "lineages", "BRAF mutant", "KRAS mutant"],
    "value": [len(annotations), annotations["lineage"].nunique(),
              int(annotations["braf_mut"].sum()), int(annotations["kras_mut"].sum())],
}))
display(coverage)
assert coverage.loc[coverage["row_type"].eq("panel_union"), "defined_symbols"].iloc[0] == 1039
assert coverage.loc[coverage["row_type"].eq("panel_union"), "mapped_symbols"].iloc[0] == 1001
"""
        ),
        nbf.v4.new_markdown_cell("## Results"),
        nbf.v4.new_code_cell(
            """display(Image(filename=str(root / "results/figures/ablation_performance.png")))
display(Image(filename=str(root / "results/figures/complexity_vs_performance.png")))
display(variants[[
    "model", "family", "parameter_count_mean", "rmse_delta_mean", "rmse_gain_vs_b4_mean"
]])
"""
        ),
        nbf.v4.new_markdown_cell("### Pathway enrichment and fold stability"),
        nbf.v4.new_code_cell(
            """significant = (
    enrichment.assign(significant=enrichment["fdr_bh"].lt(0.05))
    .groupby(["collection", "direction"])["significant"].sum().unstack(fill_value=0)
)
display(significant)
display(stability)
display(Image(filename=str(root / "results/figures/component_diagnostics.png")))
print("Mean fold-pair response-subspace overlap:", stability["mean_squared_cosine"].mean())
"""
        ),
        nbf.v4.new_markdown_cell(
            f"""## Takeaways

- control dimension은 10–30 사이에서 차이가 작고, rank 20을 30–50으로 늘려도
  B4를 따라잡지 못한다. W6 상한 선택이 큰 미탐색 이득을 숨겼다는 근거는 없다.
- pathway panel은 평균 격차를 줄였지만 B4 대비 CI
  `{pathway['rmse_gain_vs_b4_ci95_low']:.6f}–{pathway['rmse_gain_vs_b4_ci95_high']:.6f}`가
  0을 포함한다.
- BRAF/KRAS 추가 변형이 고정 low-rank 변형 중 평균 RMSE가 가장 낮지만, B4 대비 CI
  `{mutations['rmse_gain_vs_b4_ci95_low']:.6f}–{mutations['rmse_gain_vs_b4_ci95_high']:.6f}`가
  넓게 0을 포함한다. test 결과로 고른 관찰이므로 새 주 모델로 승격하지 않는다.
- W6 response subspace의 fold-pair mean squared cosine은 평균 약 `0.676`이다.
  상위 방향은 안정적이지만 최약 방향은 거의 정렬되지 않아 component 번호별 강한 해석은 피한다.
- 사전 정의 pathway enrichment는 반복 검출되지만, 예측 정확도 개선과 생물학적
  enrichment를 동일한 주장으로 취급하지 않는다."""
        ),
    ]
    nbf.write(_notebook(cells), root / "notebooks" / "05_ablation.ipynb")


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    (root / "notebooks").mkdir(exist_ok=True)
    build_data_audit(root)
    build_qc(root)
    build_response_landscape(root)
    build_baselines(root)
    build_main_model(root)
    build_ablation(root)


if __name__ == "__main__":
    main()
