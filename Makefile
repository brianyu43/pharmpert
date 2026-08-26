.PHONY: install smoke test lint data-probe core-audit

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
