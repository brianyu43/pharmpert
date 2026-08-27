"""Independent validation of W9 biological and prediction-error outputs."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from yakseopdong.landscape import _aligned_rows
from yakseopdong.metrics import spearman_or_nan
from yakseopdong.pathways import benjamini_hochberg
from yakseopdong.pseudobulk import read_vector_parquet


def validate_biology(root: Path) -> dict[str, object]:
    """Recompute W9 effects and enforce interpretation/leakage contracts."""
    associations = pd.read_csv(root / "results" / "tables" / "biological_validation.csv")
    features = pd.read_csv(root / "results" / "tables" / "biological_line_features.csv")
    targets = pd.read_csv(root / "results" / "tables" / "target_gene_validation.csv")
    enrichment = pd.read_csv(
        root / "results" / "tables" / "biological_pathway_enrichment.csv"
    )
    cases = pd.read_csv(root / "results" / "tables" / "prediction_cases.csv")
    classification = pd.read_csv(root / "results" / "tables" / "error_classification.csv")
    if len(features) != 94 or features["cell_line"].nunique() != 94:
        raise ValueError("W9 line features must contain the frozen 94-line cohort")
    if len(associations) != 38 or not associations["interpretation_only"].all():
        raise ValueError("W9 association registry is incomplete or not interpretation-only")
    if len(targets) != 15 or set(cases["case_type"]) != {"best", "median", "worst"}:
        raise ValueError("W9 target genes or prediction cases are incomplete")
    if cases.groupby("case_type").size().ne(5).any():
        raise ValueError("each W9 case group must contain five lines")
    if len(classification) != 94 or int(classification["high_error"].sum()) != 24:
        raise ValueError("W9 error classification differs from the quartile contract")

    finite = associations["p_value"].notna()
    recomputed_fdr = benjamini_hochberg(associations.loc[finite, "p_value"].to_numpy())
    max_fdr_error = float(
        np.max(np.abs(recomputed_fdr - associations.loc[finite, "fdr_bh"].to_numpy()))
    )
    if max_fdr_error > 1e-12:
        raise ValueError("W9 association FDR values do not recompute")
    enrichment_fdr = benjamini_hochberg(enrichment["p_value"].to_numpy())
    max_enrichment_fdr_error = float(np.max(np.abs(enrichment_fdr - enrichment["fdr_bh"])))
    if max_enrichment_fdr_error > 1e-12:
        raise ValueError("W9 enrichment FDR values do not recompute")

    response_meta, response = read_vector_parquet(
        root / "data" / "processed" / "response_24h.parquet", "delta_log1p_cpm"
    )
    lines = features["cell_line"].astype(str).tolist()
    aligned = _aligned_rows(response_meta, response, lines)
    response_rms = np.sqrt(np.mean(np.square(aligned), axis=1))
    max_response_rms_error = float(np.max(np.abs(response_rms - features["response_rms"])))
    if max_response_rms_error > 1e-7:
        raise ValueError("W9 response magnitude was not computed from the frozen response")

    row = associations.loc[
        associations["family"].eq("component_sensitivity")
        & associations["outcome"].eq("response_pc1")
    ].iloc[0]
    recomputed_effect = spearman_or_nan(
        features["trametinib_sensitivity"], features["response_pc1"]
    )
    primary_effect_error = abs(recomputed_effect - float(row.effect))
    if primary_effect_error > 1e-12:
        raise ValueError("W9 primary component association effect does not recompute")

    required_figures = [
        root / "results" / "figures" / "biological_correlates.png",
        root / "results" / "figures" / "prediction_error_drivers.png",
    ]
    if any(not path.exists() or path.stat().st_size < 20_000 for path in required_figures):
        raise ValueError("a required W9 figure is absent or unexpectedly small")
    report = {
        "stage": "W9_biological_interpretation",
        "status": "passed",
        "cell_lines": len(features),
        "association_tests": len(associations),
        "target_genes": len(targets),
        "prediction_cases": len(cases),
        "high_error_lines": int(classification["high_error"].sum()),
        "max_fdr_recomputation_error": max_fdr_error,
        "max_enrichment_fdr_recomputation_error": max_enrichment_fdr_error,
        "max_response_rms_recomputation_error": max_response_rms_error,
        "primary_effect_recomputation_error": primary_effect_error,
        "interpretation_only": True,
        "sensitivity_used_as_predictor": False,
    }
    output = root / "results" / "logs" / "biology_validation.json"
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report
