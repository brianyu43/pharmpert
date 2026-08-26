"""Probe the pertpy McFarland dataset and save auditable metadata."""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np


def _matrix_summary(matrix: Any) -> dict[str, object]:
    summary: dict[str, object] = {
        "type": f"{type(matrix).__module__}.{type(matrix).__name__}",
        "shape": list(matrix.shape),
        "dtype": str(matrix.dtype),
    }
    if hasattr(matrix, "nnz"):
        summary["nnz"] = int(matrix.nnz)
        summary["sparse"] = True
        values = matrix.data
    else:
        summary["sparse"] = False
        values = np.asarray(matrix).ravel()

    sample = np.asarray(values[: min(len(values), 100_000)])
    summary["sample_min"] = float(sample.min()) if sample.size else None
    summary["sample_max"] = float(sample.max()) if sample.size else None
    summary["sample_integer_like"] = (
        bool(np.allclose(sample, np.round(sample))) if sample.size else None
    )
    return summary


def run_probe(root: Path, download: bool) -> dict[str, object]:
    """Load the dataset only with explicit download authorization and summarize it."""
    if not download:
        raise ValueError("dataset loading requires the explicit --download flag")

    import pertpy as pt

    adata = pt.dt.mcfarland_2020()
    obs_summary = {
        column: {
            "dtype": str(adata.obs[column].dtype),
            "n_unique": int(adata.obs[column].nunique(dropna=True)),
            "n_missing": int(adata.obs[column].isna().sum()),
            "examples": [str(value) for value in adata.obs[column].dropna().unique()[:20]],
        }
        for column in adata.obs.columns
    }

    report = {
        "created_at_utc": datetime.now(UTC).isoformat(),
        "shape": list(adata.shape),
        "x": _matrix_summary(adata.X),
        "layers": {key: _matrix_summary(value) for key, value in adata.layers.items()},
        "raw_present": adata.raw is not None,
        "obs": obs_summary,
        "var_columns": list(adata.var.columns),
        "obs_names_unique": bool(adata.obs_names.is_unique),
        "var_names_unique": bool(adata.var_names.is_unique),
    }

    output = root / "data" / "interim" / "mcfarland_probe.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    return report
