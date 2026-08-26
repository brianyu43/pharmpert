"""Sparse pseudobulk construction for the authoritative MIX-seq experiment 3 data."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy import sparse
from scipy.io import mmread

from yakseopdong.core_audit import ARCHIVES, _member_name, load_archive_metadata

SCALE_FACTOR = 1_000_000.0


@dataclass(frozen=True)
class PseudobulkBlock:
    """Cell-line pseudobulk counts and source-QC summaries for one archive."""

    source_key: str
    condition: str
    time_hours: int
    genes: pd.DataFrame
    cell_lines: tuple[str, ...]
    depmap_ids: tuple[str, ...]
    counts: sparse.csr_matrix
    n_cells: np.ndarray
    cell_qc: dict[str, float | int | str]


def read_genes(archive: ZipFile) -> pd.DataFrame:
    """Read the ordered Ensembl ID / gene-symbol table from an archive."""
    with archive.open(_member_name(archive, "/genes.tsv")) as handle:
        genes = pd.read_csv(handle, sep="\t", header=None, names=["gene_id", "gene_symbol"])
    if genes["gene_id"].isna().any() or genes["gene_id"].duplicated().any():
        raise ValueError("gene_id must be complete and unique")
    return genes


def aggregate_normal_cells(
    matrix_gene_by_cell: sparse.spmatrix,
    classifications: pd.DataFrame,
) -> tuple[tuple[str, ...], tuple[str, ...], sparse.csr_matrix, np.ndarray]:
    """Sum raw cell counts by cell line after filtering to normal-quality cells."""
    if matrix_gene_by_cell.shape[1] != len(classifications):
        raise ValueError("matrix cell dimension and classifications row count differ")
    if not np.issubdtype(matrix_gene_by_cell.dtype, np.integer):
        raise ValueError(f"expected integer raw counts, found {matrix_gene_by_cell.dtype}")

    required = {"singlet_ID", "DepMap_ID", "cell_quality"}
    missing = required - set(classifications.columns)
    if missing:
        raise ValueError(f"missing classification columns: {sorted(missing)}")

    normal_mask = classifications["cell_quality"].eq("normal").to_numpy()
    normal = classifications.loc[normal_mask, ["singlet_ID", "DepMap_ID"]].copy()
    if normal.isna().any().any():
        raise ValueError("normal cells contain missing cell-line or DepMap identifiers")

    identity_counts = normal.groupby("singlet_ID")["DepMap_ID"].nunique()
    if (identity_counts != 1).any():
        raise ValueError("a cell line maps to multiple DepMap identifiers")

    cell_lines = tuple(sorted(normal["singlet_ID"].astype(str).unique()))
    line_index = {name: idx for idx, name in enumerate(cell_lines)}
    group_codes = normal["singlet_ID"].astype(str).map(line_index).to_numpy(dtype=np.int64)
    normal_indices = np.flatnonzero(normal_mask)

    matrix_csc = matrix_gene_by_cell.tocsc()
    normal_counts = matrix_csc[:, normal_indices]
    indicator = sparse.csr_matrix(
        (
            np.ones(len(normal_indices), dtype=np.int8),
            (group_codes, np.arange(len(normal_indices), dtype=np.int64)),
        ),
        shape=(len(cell_lines), len(normal_indices)),
    )
    pseudobulk = (indicator @ normal_counts.T).tocsr()
    n_cells = np.bincount(group_codes, minlength=len(cell_lines)).astype(np.int64)

    identities = normal.drop_duplicates("singlet_ID").set_index("singlet_ID")["DepMap_ID"]
    depmap_ids = tuple(identities.loc[list(cell_lines)].astype(str))
    return cell_lines, depmap_ids, pseudobulk, n_cells


def _cell_qc_summary(
    matrix_gene_by_cell: sparse.spmatrix,
    classifications: pd.DataFrame,
    genes: pd.DataFrame,
    source_key: str,
) -> dict[str, float | int | str]:
    normal_indices = np.flatnonzero(classifications["cell_quality"].eq("normal").to_numpy())
    normal_counts = matrix_gene_by_cell.tocsc()[:, normal_indices]
    library_sizes = np.asarray(normal_counts.sum(axis=0)).ravel().astype(np.float64)
    detected_genes = np.diff(normal_counts.indptr).astype(np.float64)
    mitochondrial = genes["gene_symbol"].astype(str).str.upper().str.startswith("MT-").to_numpy()
    mitochondrial_counts = np.asarray(normal_counts[mitochondrial, :].sum(axis=0)).ravel()
    mitochondrial_fraction = np.divide(
        mitochondrial_counts,
        library_sizes,
        out=np.zeros_like(library_sizes),
        where=library_sizes > 0,
    )

    return {
        "source": source_key,
        "normal_cells": int(len(normal_indices)),
        "zero_library_cells": int((library_sizes == 0).sum()),
        "library_size_median": float(np.median(library_sizes)),
        "library_size_p05": float(np.quantile(library_sizes, 0.05)),
        "library_size_p95": float(np.quantile(library_sizes, 0.95)),
        "detected_genes_median": float(np.median(detected_genes)),
        "detected_genes_p05": float(np.quantile(detected_genes, 0.05)),
        "detected_genes_p95": float(np.quantile(detected_genes, 0.95)),
        "mitochondrial_fraction_median": float(np.median(mitochondrial_fraction)),
        "mitochondrial_fraction_p95": float(np.quantile(mitochondrial_fraction, 0.95)),
    }


def load_archive_pseudobulk(
    path: Path,
    source_key: str,
    condition: str,
    time_hours: int,
) -> PseudobulkBlock:
    """Load one compressed Matrix Market archive and aggregate before densifying."""
    classifications, _ = load_archive_metadata(path)
    with ZipFile(path) as archive:
        genes = read_genes(archive)
        with archive.open(_member_name(archive, "/matrix.mtx")) as handle:
            raw_matrix = mmread(handle)

    if raw_matrix.shape != (len(genes), len(classifications)):
        raise ValueError(
            f"matrix shape {raw_matrix.shape} does not match "
            f"{len(genes)} genes x {len(classifications)} cells"
        )
    if not np.issubdtype(raw_matrix.dtype, np.integer):
        raise ValueError(f"matrix is not integer-valued: {raw_matrix.dtype}")
    if raw_matrix.data.size and (raw_matrix.data < 0).any():
        raise ValueError("raw count matrix contains negative values")

    cell_qc = _cell_qc_summary(raw_matrix, classifications, genes, source_key)
    cell_lines, depmap_ids, counts, n_cells = aggregate_normal_cells(
        raw_matrix, classifications
    )
    return PseudobulkBlock(
        source_key=source_key,
        condition=condition,
        time_hours=time_hours,
        genes=genes,
        cell_lines=cell_lines,
        depmap_ids=depmap_ids,
        counts=counts,
        n_cells=n_cells,
        cell_qc=cell_qc,
    )


def align_block(
    block: PseudobulkBlock, target_lines: list[str]
) -> tuple[sparse.csr_matrix, np.ndarray, list[str]]:
    """Align a block to a frozen ordered cell-line cohort."""
    positions = {line: idx for idx, line in enumerate(block.cell_lines)}
    missing = [line for line in target_lines if line not in positions]
    if missing:
        raise ValueError(f"{block.source_key} is missing cell lines: {missing}")
    indices = np.asarray([positions[line] for line in target_lines], dtype=np.int64)
    return (
        block.counts[indices, :].tocsr(),
        block.n_cells[indices],
        [block.depmap_ids[idx] for idx in indices],
    )


def log1p_cpm(counts: sparse.spmatrix | np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Convert cell-line pseudobulk counts to deterministic log1p(CPM)."""
    dense = counts.toarray() if sparse.issparse(counts) else np.asarray(counts)
    dense = dense.astype(np.float32, copy=True)
    library_sizes = dense.sum(axis=1, dtype=np.float64)
    if (library_sizes <= 0).any():
        raise ValueError("pseudobulk rows must have positive library sizes")
    dense *= (SCALE_FACTOR / library_sizes).astype(np.float32)[:, None]
    np.log1p(dense, out=dense)
    return dense, library_sizes


def write_vector_parquet(
    path: Path,
    metadata: pd.DataFrame,
    values: np.ndarray,
    value_name: str,
) -> None:
    """Write one fixed-length float vector per analytical row."""
    if len(metadata) != len(values):
        raise ValueError("metadata and vector row counts differ")
    values = np.asarray(values, dtype=np.float32)
    flat = pa.array(values.reshape(-1), type=pa.float32())
    vector = pa.FixedSizeListArray.from_arrays(flat, list_size=values.shape[1])
    table = pa.Table.from_pandas(metadata.reset_index(drop=True), preserve_index=False)
    table = table.append_column(value_name, vector)
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, path, compression="zstd")


def read_vector_parquet(path: Path, value_name: str) -> tuple[pd.DataFrame, np.ndarray]:
    """Read a vector-valued parquet written by :func:`write_vector_parquet`."""
    table = pq.read_table(path)
    vector_column = table[value_name].combine_chunks()
    values = np.asarray(vector_column.values).reshape(len(table), vector_column.type.list_size)
    metadata = table.drop([value_name]).to_pandas()
    return metadata, values.astype(np.float32, copy=False)


def _metadata_frame(
    cell_lines: list[str],
    depmap_ids: list[str],
    condition: str,
    time_hours: int,
    n_cells: np.ndarray,
    library_sizes: np.ndarray,
    control_source: str,
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "cell_line": cell_lines,
            "depmap_id": depmap_ids,
            "condition": condition,
            "time_hours": time_hours,
            "n_cells": n_cells.astype(np.int64),
            "library_size": library_sizes.astype(np.int64),
            "control_source": control_source,
        }
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_pseudobulk(root: Path) -> dict[str, object]:
    """Build frozen primary/sensitivity pseudobulk matrices and QC evidence."""
    source_specs = {
        "dmso_6h": ("control", 6),
        "dmso_24h": ("control", 24),
        "trametinib_24h": ("trametinib", 24),
    }
    blocks: dict[str, PseudobulkBlock] = {}
    for key, filename in ARCHIVES.items():
        condition, time_hours = source_specs[key]
        blocks[key] = load_archive_pseudobulk(
            root / "data" / "raw" / filename,
            source_key=key,
            condition=condition,
            time_hours=time_hours,
        )

    reference_genes = blocks["dmso_24h"].genes.reset_index(drop=True)
    for key, block in blocks.items():
        if not block.genes.reset_index(drop=True).equals(reference_genes):
            raise ValueError(f"ordered gene metadata differs for {key}")

    cohort = pd.read_csv(root / "cell_count_matrix.csv")
    primary = sorted(
        cohort.loc[cohort["primary_strict_eligible"], "cell_line"].astype(str).tolist()
    )
    sensitivity = sorted(
        cohort.loc[cohort["reproduction_pooled_eligible"], "cell_line"].astype(str).tolist()
    )
    if len(primary) != 94 or len(sensitivity) != 97:
        raise ValueError("frozen cohort counts differ from the audited 94/97 contract")

    processed = root / "data" / "processed"
    processed.mkdir(parents=True, exist_ok=True)
    pq.write_table(
        pa.Table.from_pandas(reference_genes, preserve_index=False),
        processed / "gene_metadata.parquet",
        compression="zstd",
    )

    dmso24_counts, dmso24_n, dmso24_depmap = align_block(blocks["dmso_24h"], primary)
    tram24_counts, tram24_n, tram24_depmap = align_block(blocks["trametinib_24h"], primary)
    if dmso24_depmap != tram24_depmap:
        raise ValueError("DepMap identifiers differ between primary conditions")
    dmso24_log, dmso24_library = log1p_cpm(dmso24_counts)
    tram24_log, tram24_library = log1p_cpm(tram24_counts)
    primary_values = np.vstack([dmso24_log, tram24_log])
    primary_meta = pd.concat(
        [
            _metadata_frame(
                primary,
                dmso24_depmap,
                "control",
                24,
                dmso24_n,
                dmso24_library,
                "DMSO_24hr_expt3",
            ),
            _metadata_frame(
                primary,
                tram24_depmap,
                "trametinib",
                24,
                tram24_n,
                tram24_library,
                "DMSO_24hr_expt3",
            ),
        ],
        ignore_index=True,
    )
    write_vector_parquet(
        processed / "pseudobulk_24h.parquet", primary_meta, primary_values, "log1p_cpm"
    )
    response_meta = primary_meta.iloc[: len(primary)][
        ["cell_line", "depmap_id"]
    ].copy()
    response_meta["treated_source"] = "Trametinib_24hr_expt3"
    response_meta["control_source"] = "DMSO_24hr_expt3"
    response_meta["normalization"] = "log1p_cpm_1e6"
    write_vector_parquet(
        processed / "response_24h.parquet",
        response_meta,
        tram24_log - dmso24_log,
        "delta_log1p_cpm",
    )

    dmso6_counts, dmso6_n, dmso6_depmap = align_block(blocks["dmso_6h"], sensitivity)
    dmso24_all_counts, dmso24_all_n, dmso24_all_depmap = align_block(
        blocks["dmso_24h"], sensitivity
    )
    tram24_all_counts, tram24_all_n, tram24_all_depmap = align_block(
        blocks["trametinib_24h"], sensitivity
    )
    if not (dmso6_depmap == dmso24_all_depmap == tram24_all_depmap):
        raise ValueError("DepMap identifiers differ across sensitivity conditions")
    dmso6_log, dmso6_library = log1p_cpm(dmso6_counts)
    dmso24_all_log, dmso24_all_library = log1p_cpm(dmso24_all_counts)
    control_time_meta = pd.concat(
        [
            _metadata_frame(
                sensitivity,
                dmso6_depmap,
                "control",
                6,
                dmso6_n,
                dmso6_library,
                "DMSO_6hr_expt3",
            ),
            _metadata_frame(
                sensitivity,
                dmso24_all_depmap,
                "control",
                24,
                dmso24_all_n,
                dmso24_all_library,
                "DMSO_24hr_expt3",
            ),
        ],
        ignore_index=True,
    )
    write_vector_parquet(
        processed / "pseudobulk_control_time.parquet",
        control_time_meta,
        np.vstack([dmso6_log, dmso24_all_log]),
        "log1p_cpm",
    )

    pooled_counts = dmso6_counts + dmso24_all_counts
    pooled_n = dmso6_n + dmso24_all_n
    pooled_log, pooled_library = log1p_cpm(pooled_counts)
    tram24_all_log, tram24_all_library = log1p_cpm(tram24_all_counts)
    sensitivity_meta = pd.concat(
        [
            _metadata_frame(
                sensitivity,
                dmso6_depmap,
                "control",
                24,
                pooled_n,
                pooled_library,
                "DMSO_6hr_expt3+DMSO_24hr_expt3",
            ),
            _metadata_frame(
                sensitivity,
                tram24_all_depmap,
                "trametinib",
                24,
                tram24_all_n,
                tram24_all_library,
                "DMSO_6hr_expt3+DMSO_24hr_expt3",
            ),
        ],
        ignore_index=True,
    )
    write_vector_parquet(
        processed / "pseudobulk_pooled_sensitivity.parquet",
        sensitivity_meta,
        np.vstack([pooled_log, tram24_all_log]),
        "log1p_cpm",
    )
    pooled_response_meta = sensitivity_meta.iloc[: len(sensitivity)][
        ["cell_line", "depmap_id"]
    ].copy()
    pooled_response_meta["treated_source"] = "Trametinib_24hr_expt3"
    pooled_response_meta["control_source"] = "DMSO_6hr_expt3+DMSO_24hr_expt3"
    pooled_response_meta["normalization"] = "log1p_cpm_1e6"
    write_vector_parquet(
        processed / "response_pooled_sensitivity.parquet",
        pooled_response_meta,
        tram24_all_log - pooled_log,
        "delta_log1p_cpm",
    )

    from yakseopdong.qc import run_qc

    qc_summary = run_qc(
        root=root,
        genes=reference_genes,
        primary_lines=primary,
        sensitivity_lines=sensitivity,
        dmso6_counts=dmso6_counts,
        dmso24_counts=dmso24_all_counts,
        dmso6_log=dmso6_log,
        dmso24_log=dmso24_all_log,
        tram24_log=tram24_all_log,
        cell_qc=pd.DataFrame([block.cell_qc for block in blocks.values()]),
    )

    output_paths = sorted(processed.glob("*.parquet"))
    output_contract = {
        "gene_metadata.parquet": (32_738, 0, ""),
        "pseudobulk_24h.parquet": (188, 32_738, "log1p_cpm"),
        "response_24h.parquet": (94, 32_738, "delta_log1p_cpm"),
        "pseudobulk_control_time.parquet": (194, 32_738, "log1p_cpm"),
        "pseudobulk_pooled_sensitivity.parquet": (194, 32_738, "log1p_cpm"),
        "response_pooled_sensitivity.parquet": (97, 32_738, "delta_log1p_cpm"),
    }
    processed_manifest = pd.DataFrame(
        [
            {
                "path": str(path.relative_to(root)),
                "row_count": output_contract[path.name][0],
                "vector_length": output_contract[path.name][1],
                "value_column": output_contract[path.name][2],
                "bytes": path.stat().st_size,
                "sha256": _sha256(path),
                "genes_tsv_sha256": (
                    "8778dd78085086be6dcafa8334f95f7cf76b25526e8ef2b63186db5127a1492c"
                ),
            }
            for path in output_paths
        ]
    )
    processed_manifest.to_csv(root / "processed_manifest.csv", index=False)
    report = {
        "normalization": "log1p(CPM), scale_factor=1e6",
        "gene_count": int(len(reference_genes)),
        "primary_lines": len(primary),
        "sensitivity_lines": len(sensitivity),
        "primary_rows": int(len(primary_meta)),
        "primary_response_rows": int(len(response_meta)),
        "raw_matrices_integer_nonnegative": True,
        "cell_qc": [block.cell_qc for block in blocks.values()],
        "qc": qc_summary,
        "outputs": {
            row["path"]: {"bytes": row["bytes"], "sha256": row["sha256"]}
            for row in processed_manifest.to_dict(orient="records")
        },
    }
    log_path = root / "results" / "logs" / "pseudobulk_summary.json"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    return report
