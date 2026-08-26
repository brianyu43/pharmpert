"""Runtime and repository smoke checks."""

from __future__ import annotations

import json
import platform
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import yaml

REQUIRED_PATHS = (
    "README.md",
    "pyproject.toml",
    "config/project.yaml",
    "config/cohort.yaml",
    "config/splits.yaml",
    "config/models.yaml",
    "data_manifest.csv",
    "docs/project_blueprint.md",
    "docs/scope.md",
    "docs/evaluation_protocol.md",
    "cell_count_matrix.csv",
)


def package_version(name: str) -> str | None:
    """Return an installed distribution version, or None when unavailable."""
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def repository_root() -> Path:
    """Find the project root from the installed source tree."""
    return Path(__file__).resolve().parents[2]


def smoke_report(root: Path | None = None) -> dict[str, object]:
    """Validate Stage 0 files and return a machine-readable report."""
    project_root = root or repository_root()
    missing = [path for path in REQUIRED_PATHS if not (project_root / path).exists()]

    with (project_root / "config/project.yaml").open(encoding="utf-8") as handle:
        project_config = yaml.safe_load(handle)

    return {
        "ok": not missing,
        "project_root": str(project_root),
        "missing_paths": missing,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "project": project_config["project"],
        "packages": {
            name: package_version(name)
            for name in ("anndata", "numpy", "pandas", "pertpy", "scikit-learn", "scipy")
        },
    }


def print_smoke_report(root: Path | None = None) -> int:
    """Print the smoke report and return a process exit code."""
    report = smoke_report(root)
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0 if report["ok"] else 1
