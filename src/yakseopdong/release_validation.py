"""Independent validation of the final frozen release surface."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from yakseopdong.pseudobulk import read_vector_parquet
from yakseopdong.release import sha256


def _validate_stage_logs(root: Path) -> dict[str, str]:
    names = ["cclr", "ablation", "temporal", "biology", "robustness", "distribution"]
    statuses: dict[str, str] = {}
    for name in names:
        path = root / "results" / "logs" / f"{name}_validation.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        status = str(payload.get("status", "")).lower()
        if status not in {"pass", "passed"}:
            raise ValueError(f"{name} validator has non-passing status: {status}")
        statuses[name] = status
    return statuses


def _validate_manifest(root: Path, manifest: pd.DataFrame, *, is_table: bool) -> None:
    for row in manifest.itertuples(index=False):
        path = root / str(row.path)
        if not path.is_file():
            raise ValueError(f"manifest path is absent: {row.path}")
        if int(row.bytes) != path.stat().st_size:
            raise ValueError(f"manifest byte count differs for {row.path}")
        if str(row.sha256) != sha256(path):
            raise ValueError(f"manifest SHA-256 differs for {row.path}")
        if is_table and int(row.rows) != len(pd.read_csv(path)):
            raise ValueError(f"manifest row count differs for {row.path}")


def _require_report_sections(root: Path) -> None:
    report = (root / "report" / "final_report.md").read_text(encoding="utf-8")
    supplementary = (root / "report" / "supplementary.md").read_text(encoding="utf-8")
    limitations = (root / "report" / "limitations.md").read_text(encoding="utf-8")
    required = [
        "# Predicting Context-Dependent Transcriptional Responses",
        "## Technical Summary",
        "## Key Findings",
        "## Scope, Data, and Metrics",
        "## Methodology and Model Specification",
        "## Results",
        "## Validation and Reproducibility",
        "## Limitations and Uncertainty",
        "## Conclusions",
        "## Recommended Next Steps",
        "## Further Questions",
    ]
    missing = [heading for heading in required if heading not in report]
    if missing:
        raise ValueError(f"final report is missing sections: {missing}")
    if "작성 전" in report + supplementary + limitations:
        raise ValueError("a report surface is still marked as not written")
    if "0.001715" not in report or "sampling-sensitive" not in report:
        raise ValueError("final report omits the magnitude or sampling-sensitivity claim")
    artifact = json.loads((root / "report" / "artifact.json").read_text(encoding="utf-8"))
    receipt = json.loads(
        (root / "release" / "report_artifact_validation.json").read_text(encoding="utf-8")
    )
    if artifact.get("surface") != "report" or artifact.get("snapshot", {}).get("status") != "ready":
        raise ValueError("bounded report artifact is not release-ready")
    if len(artifact.get("snapshot", {}).get("datasets", {})) != 3:
        raise ValueError("bounded report artifact must expose three reviewed datasets")
    if receipt.get("status") != "passed" or receipt.get("render_status") != "passed":
        raise ValueError("report artifact validation/render receipt is not passing")


def validate_release(root: Path) -> dict[str, object]:
    """Validate final predictions, metrics, manifests, documents, and freeze hashes."""
    stage_statuses = _validate_stage_logs(root)
    config = json.loads((root / "release" / "final_config.json").read_text(encoding="utf-8"))
    freeze = json.loads(
        (root / "results" / "logs" / "release_freeze.json").read_text(encoding="utf-8")
    )

    prediction_meta, predictions = read_vector_parquet(
        root / "results" / "final_predictions.parquet",
        "predicted_delta_log1p_cpm",
    )
    expected_models = {"B0", "B1", "B2", "B3", "B4", "CCLR"}
    model_counts = prediction_meta.groupby("model")["cell_line"].nunique().to_dict()
    if predictions.shape != (564, 32_738):
        raise ValueError("final prediction matrix differs from the 564 x 32738 contract")
    if set(model_counts) != expected_models or set(model_counts.values()) != {94}:
        raise ValueError("every final model must predict exactly 94 unique held-out lines")
    if prediction_meta[["cell_line", "model"]].duplicated().any():
        raise ValueError("final predictions contain duplicate line-model rows")
    if not prediction_meta["external_test"].all():
        raise ValueError("a final prediction is not marked as held out")
    if prediction_meta["treated_response_used_for_fit"].any():
        raise ValueError("a final prediction was marked as using treated test response")

    metrics = pd.read_csv(root / "results" / "final_metrics.csv")
    required_sections = {
        "core_benchmark",
        "temporal_external",
        "robustness",
        "single_cell_distribution",
    }
    if len(metrics) != 64 or not required_sections.issubset(set(metrics["section"])):
        raise ValueError("final metric table differs from its section/row contract")
    lines = pd.read_csv(root / "results" / "final_cell_lines.csv")
    if len(lines) != 94 or lines["cell_line"].nunique() != 94:
        raise ValueError("final cell-line table must contain 94 unique lines")

    figure_manifest = pd.read_csv(root / "results" / "figure_manifest.csv")
    table_manifest = pd.read_csv(root / "results" / "table_manifest.csv")
    if figure_manifest["figure_id"].tolist() != [
        "F1",
        "F2",
        "F3",
        "F4",
        "F5",
        "F6",
        "F7",
        "F8",
        "S1",
        "S2",
    ]:
        raise ValueError("figure manifest must contain ordered F1-F8 and S1-S2")
    if not figure_manifest["visual_inspection"].eq("passed").all():
        raise ValueError("a final figure has not passed visual inspection")
    if table_manifest["table_id"].tolist() != ["T1", "T2", "T3", "T4", "T5", "S1"]:
        raise ValueError("table manifest must contain ordered T1-T5 and S1")
    _validate_manifest(root, figure_manifest, is_table=False)
    _validate_manifest(root, table_manifest, is_table=True)
    _require_report_sections(root)

    hash_contract = {
        "final_predictions_sha256": root / "results" / "final_predictions.parquet",
        "final_metrics_sha256": root / "results" / "final_metrics.csv",
        "final_cell_lines_sha256": root / "results" / "final_cell_lines.csv",
        "figure_manifest_sha256": root / "results" / "figure_manifest.csv",
        "table_manifest_sha256": root / "results" / "table_manifest.csv",
        "final_report_sha256": root / "report" / "final_report.md",
        "supplementary_sha256": root / "report" / "supplementary.md",
        "limitations_sha256": root / "report" / "limitations.md",
        "report_artifact_sha256": root / "report" / "artifact.json",
        "report_artifact_validation_sha256": (
            root / "release" / "report_artifact_validation.json"
        ),
    }
    artifact_hashes = {name: sha256(path) for name, path in hash_contract.items()}
    for name, actual in artifact_hashes.items():
        if str(config.get(name)) != actual or str(freeze.get(name)) != actual:
            raise ValueError(f"release hash contract differs for {name}")
    scalar_contract = {
        "release_version": "1.0.0",
        "evaluation_protocol_version": "1.3",
        "scope_version": "1.1",
        "seed": 20260827,
        "primary_cohort_lines": 94,
        "gene_count": 32_738,
        "prediction_rows": 564,
        "final_metrics_rows": 64,
        "figure_count": 10,
        "table_count": 6,
    }
    for name, expected in scalar_contract.items():
        if config.get(name) != expected or freeze.get(name) != expected:
            raise ValueError(f"release scalar contract differs for {name}")

    report = {
        "stage": "W12_W15_release",
        "status": "passed",
        "stage_validators": stage_statuses,
        "prediction_shape": list(predictions.shape),
        "models": sorted(model_counts),
        "lines_per_model": model_counts,
        "final_metrics_rows": len(metrics),
        "final_cell_lines": lines["cell_line"].nunique(),
        "figure_count": len(figure_manifest),
        "table_count": len(table_manifest),
        "reports_complete": True,
        "treated_test_response_used_for_fit": False,
        "artifact_hashes": artifact_hashes,
    }
    output = root / "results" / "logs" / "release_validation.json"
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report
