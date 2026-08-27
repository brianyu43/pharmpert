"""Build the bounded Data Analytics report manifest from frozen release metrics."""

from __future__ import annotations

import json
from pathlib import Path


def source(identifier: str, label: str, path: str, sql: str, description: str) -> dict:
    return {
        "id": identifier,
        "label": label,
        "path": path,
        "query": {
            "sql": sql,
            "description": description,
            "engine": "duckdb",
            "language": "sql",
            "tables_used": [path],
        },
    }


def build_payload() -> dict:
    metrics_source = source(
        "final_metrics",
        "Frozen final metrics",
        "results/final_metrics.csv",
        "SELECT model, estimate AS rmse, ci95_low, ci95_high, n "
        "FROM read_csv_auto('results/final_metrics.csv') "
        "WHERE section = 'core_benchmark' AND metric = 'rmse_delta' "
        "AND model IN ('B1_MEAN_W5','B4_DIRECT_RIDGE_W5','CCLR_NESTED_W6')",
        "Select the three frozen primary benchmark RMSE rows.",
    )
    metrics_source["query"].update(
        {
            "filters": ["section = core_benchmark", "metric = rmse_delta"],
            "metric_definitions": [
                "rmse = unweighted mean across held-out cell-line gene-wise RMSE"
            ],
        }
    )
    temporal_source = source(
        "temporal_metrics",
        "Frozen temporal metrics",
        "results/final_metrics.csv",
        "SELECT CAST(regexp_extract(analysis, '[0-9]+') AS INTEGER) AS time_hours, "
        "estimate AS gain, ci95_low, ci95_high, n "
        "FROM read_csv_auto('results/final_metrics.csv') "
        "WHERE section = 'temporal_external' AND model = 'B4_FIXED_D20_A100' "
        "AND metric = 'rmse_gain_vs_b1' ORDER BY time_hours",
        "Select B4 gain versus B1 across frozen external time points.",
    )
    temporal_source["query"].update(
        {
            "filters": ["section = temporal_external", "17 non-overlapping external lines"],
            "metric_definitions": ["gain = RMSE(B1) - RMSE(B4); positive is better"],
        }
    )
    robustness_source = source(
        "robustness",
        "Frozen robustness conclusions",
        "results/tables/final_table5_robustness.csv",
        "SELECT conclusion AS test, evidence AS estimate, status "
        "FROM read_csv_auto('results/tables/final_table5_robustness.csv') "
        "ORDER BY conclusion",
        "Read the six frozen robustness conclusions.",
    )
    robustness_source["query"]["metric_definitions"] = [
        "status summarizes whether the primary B4 gain was robust, sensitive, "
        "inconclusive, or limited"
    ]
    title = (
        "Predicting Context-Dependent Transcriptional Responses to MEK Inhibition "
        "across Cancer Cell Lines"
    )
    manifest = {
        "version": 1,
        "surface": "report",
        "title": title,
        "description": (
            "Leakage-safe held-out-cell-line evaluation of trametinib transcriptional "
            "response prediction."
        ),
        "generatedAt": "2026-08-27T00:00:00+09:00",
        "sources": [metrics_source, temporal_source, robustness_source],
        "blocks": [
            {"id": "title", "type": "markdown", "body": f"# {title}"},
            {
                "id": "summary",
                "type": "markdown",
                "body": (
                    "## Executive Summary\nA leakage-safe baseline-state model detects a "
                    "small context-dependent signal beyond the global mean, but the "
                    "absolute gain is sampling-sensitive and does not support a strong "
                    "individualized-prediction claim."
                ),
            },
            {
                "id": "benchmark_heading",
                "type": "markdown",
                "body": (
                    "## Held-out benchmark\nB4 reduced macro RMSE from 0.323750 to "
                    "0.322035: a gain of 0.001715 (95% CI 0.001168–0.002284; about "
                    "0.53%). CCLR remained slightly worse than B4."
                ),
            },
            {"id": "benchmark_block", "type": "chart", "chartId": "benchmark_chart"},
            {
                "id": "temporal_heading",
                "type": "markdown",
                "body": (
                    "## Temporal transfer\nThe fixed 24h B4 model was worse than B1 at 3h "
                    "and 6h, inconclusive at 12h, and modestly better at 24h and 48h."
                ),
            },
            {"id": "temporal_block", "type": "chart", "chartId": "temporal_chart"},
            {
                "id": "robustness_heading",
                "type": "markdown",
                "body": (
                    "## Robustness and limitations\nThe gain survived threshold and gene-filter "
                    "changes, but reversed in all five equal-20-cell repeats. LOLO was "
                    "inconclusive and the split-half full-target noise floor was 0.279028."
                ),
            },
            {"id": "robustness_block", "type": "table", "tableId": "robustness_table"},
            {
                "id": "conclusion",
                "type": "markdown",
                "body": (
                    "## Conclusion\nBaseline transcription contains detectable context "
                    "information, but the measured improvement is too small and "
                    "sampling-sensitive to establish strong personalized prediction."
                ),
            },
        ],
        "charts": [
            {
                "id": "benchmark_chart",
                "title": "Held-out 24h macro RMSE",
                "type": "bar",
                "dataset": "benchmark",
                "source": metrics_source,
                "encodings": {
                    "x": {"field": "model", "type": "nominal"},
                    "y": {"field": "rmse", "type": "quantitative"},
                },
            },
            {
                "id": "temporal_chart",
                "title": "B4 RMSE gain relative to B1 across time",
                "type": "line",
                "dataset": "temporal",
                "source": temporal_source,
                "encodings": {
                    "x": {"field": "time_hours", "type": "quantitative"},
                    "y": {"field": "gain", "type": "quantitative"},
                },
            },
        ],
        "tables": [
            {
                "id": "robustness_table",
                "title": "Robustness conclusion matrix",
                "dataset": "robustness",
                "source": robustness_source,
                "columns": [
                    {"field": "test", "label": "Test"},
                    {"field": "estimate", "label": "Estimate"},
                    {"field": "status", "label": "Status"},
                ],
                "defaultSort": {"field": "test", "direction": "asc"},
            }
        ],
    }
    snapshot = {
        "version": 1,
        "status": "ready",
        "generatedAt": "2026-08-27T00:00:00+09:00",
        "datasets": {
            "benchmark": [
                {
                    "model": "B1 global mean",
                    "rmse": 0.323750,
                    "pcc_delta": 0.420086,
                    "pcc_context": None,
                    "gain_vs_b1": 0,
                    "n_lines": 94,
                },
                {
                    "model": "B4 direct ridge",
                    "rmse": 0.322035,
                    "pcc_delta": 0.428001,
                    "pcc_context": 0.107666,
                    "gain_vs_b1": 0.001715,
                    "n_lines": 94,
                },
                {
                    "model": "CCLR",
                    "rmse": 0.322383,
                    "pcc_delta": 0.425501,
                    "pcc_context": 0.096521,
                    "gain_vs_b1": 0.001367,
                    "n_lines": 94,
                },
            ],
            "temporal": [
                {
                    "time_hours": 3,
                    "gain": -0.007310,
                    "ci_low": -0.010007,
                    "ci_high": -0.004820,
                    "n_lines": 17,
                    "model": "B4 fixed",
                },
                {
                    "time_hours": 6,
                    "gain": -0.004712,
                    "ci_low": -0.007077,
                    "ci_high": -0.002399,
                    "n_lines": 17,
                    "model": "B4 fixed",
                },
                {
                    "time_hours": 12,
                    "gain": -0.001242,
                    "ci_low": -0.003043,
                    "ci_high": 0.000515,
                    "n_lines": 17,
                    "model": "B4 fixed",
                },
                {
                    "time_hours": 24,
                    "gain": 0.002074,
                    "ci_low": 0.000311,
                    "ci_high": 0.004036,
                    "n_lines": 17,
                    "model": "B4 fixed",
                },
                {
                    "time_hours": 48,
                    "gain": 0.003839,
                    "ci_low": 0.002176,
                    "ci_high": 0.005497,
                    "n_lines": 17,
                    "model": "B4 fixed",
                },
            ],
            "robustness": [
                {
                    "test": "Cell-line bootstrap",
                    "estimate": "gain 0.001715; 3/3 CIs > 0",
                    "status": "robust",
                },
                {
                    "test": "Inclusion thresholds",
                    "estimate": "gain 0.001470–0.001892",
                    "status": "robust",
                },
                {
                    "test": "Variable-gene filters",
                    "estimate": "gain 0.001246–0.001637",
                    "status": "robust",
                },
                {
                    "test": "Equal 20 cells",
                    "estimate": "0/5 positive repeats",
                    "status": "sensitive",
                },
                {
                    "test": "Leave one lineage out",
                    "estimate": "0.000514 [-0.000132, 0.001178]",
                    "status": "inconclusive",
                },
                {
                    "test": "Split-half target",
                    "estimate": "floor 0.279028; PCC 0.236373",
                    "status": "limitation",
                },
            ],
        },
    }
    return {
        "surface": "report",
        "manifest": manifest,
        "snapshot": snapshot,
        "sources": [metrics_source, temporal_source, robustness_source],
    }


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    output = root / "report" / "artifact.json"
    output.write_text(json.dumps(build_payload(), indent=2), encoding="utf-8")
    print(output)
