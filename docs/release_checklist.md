# Release 1.0.0 Validation Record

Date: 2026-08-27

## Executed gates

| Gate | Command | Result |
|---|---|---|
| locked environment | `uv sync --locked --all-groups` | 165 resolved, 163 checked |
| lint | `uv run ruff check src tests scripts` | passed |
| unit tests | `uv run pytest -q` | 48 passed |
| smoke | `uv run python -m yakseopdong smoke` | passed; Python 3.13.14, macOS arm64 |
| W6 validation | `validate-cclr` | passed; metrics and fold-artifact reconstruction |
| W7 validation | `validate-ablation` | passed; 16 variants, 1,504 line-metric rows |
| W8 validation | `validate-temporal` | passed; 13,713 cells, response reconstruction error 0 |
| W9 validation | `validate-biology` | passed; 38 tests, FDR recomputed |
| W10 validation | `validate-robustness` | passed; 69 summary rows, 470 noise and subsampling rows |
| W11 validation | `validate-distribution` | passed; max metric error `3.55e-15` |
| W12–W15 validation | `validate-release` | passed; 564 × 32,738 predictions, 10 figures, 6 tables |
| notebooks | `make notebooks` | 10/10 executed in place |
| report artifact | MCP `validate_artifact`, then `render_artifact` | passed; 3 datasets, 3 sources |
| diff hygiene | `git diff --check` | passed after final formatting |

Jupyter required local kernel sockets, so notebook execution was authorized outside the restricted filesystem/network sandbox. All communication stayed on the local machine; no external data transfer was part of the run.

## Frozen hashes

- Final predictions: `5ffe460c0a3019fa6369d483fa6dfd2ec771494e472418ab401303c555ee6562`
- Final metrics: `bd1c448709a69b3858f0d91451615e8ed129f26450f574917a7df69297ab281f`
- Final cell lines: `4fd26f2484a62027ac184759b1b71ecb9e23d91a0716aa6d3e73891af78027ed`
- Figure manifest: `8bbba730263adc6f4667aab7ef2bd17389b423c4cc81360b4de4ded9a2886af2`
- Table manifest: `1f2acdc7ac55ec7be1e1aada910394d216b3cd85915eb4ed6b5b0c6257d1ca23`

## Release commits

The analysis parent is `7bd38c61bc4220a35483c47b17168c57aecd3380`. The final content and release-record commits are filled after the two-step commit freeze to avoid a self-referential commit hash.
