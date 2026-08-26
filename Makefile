.PHONY: install smoke test lint data-probe core-audit pseudobulk notebooks

install:
	uv sync --all-groups

smoke:
	uv run python -m yakseopdong smoke

test:
	uv run pytest

lint:
	uv run ruff check src tests scripts

data-probe:
	uv run python -m yakseopdong data-probe --download

core-audit:
	uv run python -m yakseopdong core-audit

pseudobulk:
	uv run python -m yakseopdong build-pseudobulk

notebooks:
	uv run python scripts/build_notebooks.py
	uv run jupyter nbconvert --execute --to notebook --inplace notebooks/00_data_audit.ipynb
	uv run jupyter nbconvert --execute --to notebook --inplace notebooks/01_qc.ipynb
