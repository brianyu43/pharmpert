.PHONY: install smoke test lint data-probe core-audit pseudobulk metadata landscape splits baselines cclr validate-cclr ablation validate-ablation temporal validate-temporal biology validate-biology robustness validate-robustness notebooks

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

metadata:
	uv run python -m yakseopdong metadata-audit

landscape:
	uv run python -m yakseopdong landscape

splits:
	uv run python -m yakseopdong build-splits

baselines:
	uv run python -m yakseopdong run-baselines

cclr:
	uv run python -m yakseopdong run-cclr

validate-cclr:
	uv run python -m yakseopdong validate-cclr

ablation:
	uv run python -m yakseopdong run-ablation

validate-ablation:
	uv run python -m yakseopdong validate-ablation

temporal:
	uv run python -m yakseopdong run-temporal

validate-temporal:
	uv run python -m yakseopdong validate-temporal

biology:
	uv run python -m yakseopdong run-biology

validate-biology:
	uv run python -m yakseopdong validate-biology

robustness:
	uv run python -m yakseopdong run-robustness

validate-robustness:
	uv run python -m yakseopdong validate-robustness

notebooks:
	uv run python scripts/build_notebooks.py
	uv run jupyter nbconvert --execute --to notebook --inplace notebooks/00_data_audit.ipynb
	uv run jupyter nbconvert --execute --to notebook --inplace notebooks/01_qc.ipynb
	uv run jupyter nbconvert --execute --to notebook --inplace notebooks/02_response_landscape.ipynb
	uv run jupyter nbconvert --execute --to notebook --inplace notebooks/03_baselines.ipynb
	uv run jupyter nbconvert --execute --to notebook --inplace notebooks/04_main_model.ipynb
	uv run jupyter nbconvert --execute --to notebook --inplace notebooks/05_ablation.ipynb
	uv run jupyter nbconvert --execute --to notebook --inplace notebooks/06_temporal.ipynb
	uv run jupyter nbconvert --execute --to notebook --inplace notebooks/07_biological_interpretation.ipynb
	uv run jupyter nbconvert --execute --to notebook --inplace notebooks/08_robustness.ipynb
