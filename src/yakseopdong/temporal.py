"""W8 time-course pseudobulk, descriptive trajectories, and external transfer."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import anndata as ad
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import yaml
from matplotlib.ticker import MaxNLocator
from numpy.typing import NDArray
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge

from yakseopdong.ablation import predict_fixed_low_rank
from yakseopdong.landscape import _aligned_rows
from yakseopdong.metrics import (
    bootstrap_mean_interval,
    pearson_or_nan,
    rmse,
    spearman_or_nan,
)
from yakseopdong.models import fit_control_embedding, fit_response_embedding
from yakseopdong.pathways import load_geneset_config
from yakseopdong.pseudobulk import (
    _sha256,
    log1p_cpm,
    read_vector_parquet,
    write_vector_parquet,
)

TIMEPOINTS = (3, 6, 12, 24, 48)
TAG_PATTERN = re.compile(r"^(DMSO|Tram|Untreated)_(3|6|12|24|48)hr$")
TRANSFER_MODELS = ("B1", "B4_FIXED_D20_A100", "CCLR_FIXED_D20_R20_A100")
SEED = 20260827


def parse_timecourse_tag(value: object) -> tuple[str, int] | None:
    """Parse an assigned hash tag into the frozen analytical condition and hour."""
    match = TAG_PATTERN.fullmatch(str(value))
    if match is None:
        return None
    source, hour_text = match.groups()
    hour = int(hour_text)
    if source == "Untreated" and hour != 48:
        return None
    condition = {"DMSO": "control", "Tram": "trametinib", "Untreated": "untreated"}[
        source
    ]
    return condition, hour


def eligibility_from_counts(
    counts: pd.DataFrame, thresholds: tuple[int, ...] = (10, 20, 30)
) -> pd.DataFrame:
    """Return one row per line with all-time matched-condition eligibility flags."""
    required = counts.loc[counts["condition"].isin(["control", "trametinib"])].copy()
    group_count = required.groupby("depmap_id", observed=True).size()
    minimum = required.groupby("depmap_id", observed=True)["n_cells"].min()
    identity = (
        counts[["cell_line", "depmap_id", "disease"]]
        .drop_duplicates("depmap_id")
        .set_index("depmap_id")
    )
    result = identity.join(group_count.rename("matched_group_count")).join(
        minimum.rename("min_cells_across_matched_groups")
    )
    for threshold in thresholds:
        result[f"eligible_t{threshold}"] = (
            result["matched_group_count"].eq(10)
            & result["min_cells_across_matched_groups"].ge(threshold)
        )
    return result.reset_index().sort_values("cell_line", ignore_index=True)


def aggregate_grouped_rows(
    matrix: Any,
    selected_indices: NDArray[np.int64],
    group_codes: NDArray[np.int64],
    n_groups: int,
    chunk_size: int = 512,
) -> NDArray[np.int64]:
    """Aggregate selected rows of a backed sparse matrix by group in gene chunks."""
    if len(selected_indices) != len(group_codes):
        raise ValueError("selected row indices and group codes differ in length")
    if len(selected_indices) == 0 or n_groups <= 0:
        raise ValueError("time-course aggregation requires selected rows and groups")
    indicator = sparse.csr_matrix(
        (
            np.ones(len(group_codes), dtype=np.int8),
            (group_codes, np.arange(len(group_codes), dtype=np.int64)),
        ),
        shape=(n_groups, len(group_codes)),
    )
    output = np.zeros((n_groups, int(matrix.shape[1])), dtype=np.int64)
    for start in range(0, int(matrix.shape[1]), chunk_size):
        stop = min(start + chunk_size, int(matrix.shape[1]))
        block = matrix[selected_indices, start:stop]
        if not sparse.issparse(block):
            block = sparse.csr_matrix(np.asarray(block))
        if block.data.size:
            if (block.data < 0).any() or not np.allclose(block.data, np.round(block.data)):
                raise ValueError("time-course X is not non-negative integer-like raw counts")
        output[:, start:stop] = (indicator @ block).toarray().astype(np.int64)
    return output


def aggregate_h5ad_csc(
    path: Path,
    selected_indices: NDArray[np.int64],
    group_codes: NDArray[np.int64],
    n_groups: int,
    chunk_size: int = 256,
) -> NDArray[np.int64]:
    """Aggregate rows directly from the on-disk H5AD CSC representation.

    AnnData's backed sparse row slicing materializes the entire 711-million-entry
    matrix for this file. Reading each contiguous CSC gene block avoids that
    behavior while preserving exact raw-count sums.
    """
    if len(selected_indices) != len(group_codes):
        raise ValueError("selected row indices and group codes differ in length")
    with h5py.File(path, "r") as handle:
        group = handle["X"]
        if str(group.attrs.get("encoding-type")) != "csc_matrix":
            raise ValueError("time-course H5AD X must be stored as CSC")
        n_rows, n_genes = (int(value) for value in group.attrs["shape"])
        indptr = np.asarray(group["indptr"], dtype=np.int64)
        row_to_group = np.full(n_rows, -1, dtype=np.int32)
        row_to_group[selected_indices] = group_codes.astype(np.int32)
        output = np.zeros((n_groups, n_genes), dtype=np.int64)
        for start in range(0, n_genes, chunk_size):
            stop = min(start + chunk_size, n_genes)
            pointer = indptr[start : stop + 1]
            first, last = int(pointer[0]), int(pointer[-1])
            row_indices = np.asarray(group["indices"][first:last], dtype=np.int32)
            mapped_groups = row_to_group[row_indices]
            keep = mapped_groups >= 0
            if not keep.any():
                continue
            lengths = np.diff(pointer)
            local_columns = np.repeat(
                np.arange(stop - start, dtype=np.int32), lengths
            )[keep]
            data = np.asarray(group["data"][first:last])[keep]
            if (data < 0).any() or not np.allclose(data, np.round(data)):
                raise ValueError("time-course X is not non-negative integer-like raw counts")
            block = sparse.coo_matrix(
                (
                    data.astype(np.int64),
                    (mapped_groups[keep], local_columns),
                ),
                shape=(n_groups, stop - start),
            ).toarray()
            output[:, start:stop] = block
    return output


def _load_protocol(root: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    with (root / "config" / "cohort.yaml").open(encoding="utf-8") as handle:
        cohort = yaml.safe_load(handle)["timecourse"]
    with (root / "config" / "models.yaml").open(encoding="utf-8") as handle:
        models = yaml.safe_load(handle)["w8_temporal"]
    return cohort, models


def _timecourse_obs(adata: ad.AnnData) -> tuple[pd.DataFrame, NDArray[np.int64]]:
    obs = adata.obs.reset_index(drop=False).rename(columns={"index": "cell_barcode"})
    # Cast categorical tags to object before returning tuple/None values; pandas
    # otherwise attempts to construct a categorical MultiIndex from the mapper.
    parsed = obs["hash_tag"].astype(object).map(parse_timecourse_tag)
    selected = obs["cell_quality"].eq("normal") & parsed.notna()
    frame = obs.loc[
        selected,
        ["cell_barcode", "singlet_ID", "DepMap_ID", "disease", "hash_tag", "ncounts"],
    ].copy()
    frame[["condition", "time_hours"]] = pd.DataFrame(
        parsed.loc[selected].tolist(), index=frame.index
    )
    frame = frame.rename(
        columns={"singlet_ID": "cell_line", "DepMap_ID": "depmap_id"}
    )
    if frame[["cell_line", "depmap_id", "condition", "time_hours"]].isna().any().any():
        raise ValueError("assigned normal time-course cells contain missing identities")
    if frame.groupby("cell_line", observed=True)["depmap_id"].nunique().ne(1).any():
        raise ValueError("a time-course cell line maps to multiple DepMap IDs")
    return frame, frame.index.to_numpy(dtype=np.int64)


def _write_temporal_manifest(root: Path, paths: list[tuple[Path, int, int, str]]) -> None:
    manifest_path = root / "processed_manifest.csv"
    manifest = pd.read_csv(manifest_path)
    relative_paths = {str(path.relative_to(root)) for path, *_ in paths}
    manifest = manifest.loc[~manifest["path"].isin(relative_paths)].copy()
    gene_hash = str(manifest["genes_tsv_sha256"].dropna().iloc[0])
    temporal = pd.DataFrame(
        [
            {
                "path": str(path.relative_to(root)),
                "row_count": rows,
                "vector_length": vector_length,
                "value_column": value_column,
                "bytes": path.stat().st_size,
                "sha256": _sha256(path),
                "genes_tsv_sha256": gene_hash,
            }
            for path, rows, vector_length, value_column in paths
        ]
    )
    pd.concat([manifest, temporal], ignore_index=True).sort_values("path").to_csv(
        manifest_path, index=False
    )


def build_timecourse_pseudobulk(root: Path) -> dict[str, object]:
    """Build time-matched pseudobulk response matrices from the integrated H5AD."""
    cohort_config, _ = _load_protocol(root)
    source = root / str(cohort_config["source"])
    adata = ad.read_h5ad(source, backed="r")
    try:
        frame, selected_indices = _timecourse_obs(adata)
        grouping_columns = ["cell_line", "depmap_id", "disease", "condition", "time_hours"]
        groups = (
            frame.groupby(grouping_columns, observed=True, sort=True)
            .size()
            .rename("n_cells")
            .reset_index()
            .sort_values(grouping_columns, ignore_index=True)
        )
        frame = frame.merge(
            groups.assign(group_code=np.arange(len(groups), dtype=np.int64)),
            on=grouping_columns,
            how="left",
            validate="many_to_one",
        )
        # The merge preserves left-row order, hence this aligns to selected_indices.
        group_codes = frame["group_code"].to_numpy(dtype=np.int64)
        genes = adata.var.reset_index(drop=False).rename(columns={"index": "h5ad_var_name"})
    finally:
        adata.file.close()
    counts = aggregate_h5ad_csc(source, selected_indices, group_codes, len(groups))

    official_genes = pq.read_table(
        root / "data" / "processed" / "gene_metadata.parquet"
    ).to_pandas()
    if len(genes) != len(official_genes):
        raise ValueError("H5AD and official gene tables differ in length")
    if "ensembl_id" not in genes or not genes["ensembl_id"].astype(str).equals(
        official_genes["gene_id"].astype(str)
    ):
        raise ValueError("H5AD Ensembl order differs from the official processed gene order")

    normalized, library_sizes = log1p_cpm(counts)
    groups["library_size"] = library_sizes.astype(np.int64)
    groups["source_tag"] = [
        ("DMSO" if condition == "control" else "Tram" if condition == "trametinib" else "Untreated")
        + f"_{hour}hr"
        for condition, hour in zip(groups["condition"], groups["time_hours"], strict=True)
    ]
    cohort = eligibility_from_counts(groups)
    groups = groups.merge(
        cohort[["depmap_id", "eligible_t10", "eligible_t20", "eligible_t30"]],
        on="depmap_id",
        validate="many_to_one",
    )

    processed = root / "data" / "processed"
    pseudobulk_path = processed / "pseudobulk_timecourse.parquet"
    write_vector_parquet(pseudobulk_path, groups, normalized, "log1p_cpm")

    response_rows: list[dict[str, object]] = []
    response_vectors: list[NDArray[np.float32]] = []
    group_positions = {
        (str(row.depmap_id), str(row.condition), int(row.time_hours)): index
        for index, row in groups.iterrows()
    }
    for _, line in cohort.iterrows():
        for hour in TIMEPOINTS:
            control_position = group_positions[(str(line.depmap_id), "control", hour)]
            treated_position = group_positions[(str(line.depmap_id), "trametinib", hour)]
            response_rows.append(
                {
                    "cell_line": line.cell_line,
                    "depmap_id": line.depmap_id,
                    "disease": line.disease,
                    "time_hours": hour,
                    "n_control_cells": int(groups.loc[control_position, "n_cells"]),
                    "n_treated_cells": int(groups.loc[treated_position, "n_cells"]),
                    "eligible_t10": bool(line.eligible_t10),
                    "eligible_t20": bool(line.eligible_t20),
                    "eligible_t30": bool(line.eligible_t30),
                    "control_source": f"DMSO_{hour}hr",
                    "treated_source": f"Tram_{hour}hr",
                    "normalization": "log1p_cpm_1e6",
                }
            )
            response_vectors.append(normalized[treated_position] - normalized[control_position])
    response_path = processed / "response_timecourse.parquet"
    write_vector_parquet(
        response_path,
        pd.DataFrame(response_rows),
        np.asarray(response_vectors),
        "delta_log1p_cpm",
    )

    groups.to_csv(root / "results" / "tables" / "timecourse_cell_counts.csv", index=False)
    cohort.to_csv(root / "results" / "tables" / "timecourse_cohort.csv", index=False)
    _write_temporal_manifest(
        root,
        [
            (pseudobulk_path, len(groups), normalized.shape[1], "log1p_cpm"),
            (response_path, len(response_rows), normalized.shape[1], "delta_log1p_cpm"),
        ],
    )
    return {
        "source": str(source.relative_to(root)),
        "source_sha256": _sha256(source),
        "normal_assigned_cells": int(len(frame)),
        "pseudobulk_groups": int(len(groups)),
        "response_rows": int(len(response_rows)),
        "cell_lines_all": int(len(cohort)),
        "cell_lines_t10": int(cohort["eligible_t10"].sum()),
        "cell_lines_t20": int(cohort["eligible_t20"].sum()),
        "cell_lines_t30": int(cohort["eligible_t30"].sum()),
        "genes": int(normalized.shape[1]),
    }


def _pathway_scores(
    genes: pd.DataFrame,
    metadata: pd.DataFrame,
    response: NDArray[np.floating],
    genesets: dict[str, Any],
) -> pd.DataFrame:
    symbols = genes["gene_symbol"].astype(str)
    rows: list[pd.DataFrame] = []
    for name, collection in genesets["collections"].items():
        indices = np.flatnonzero(symbols.isin(set(collection["genes"])).to_numpy())
        if not len(indices):
            continue
        block = metadata[["cell_line", "depmap_id", "disease", "time_hours"]].copy()
        block["pathway"] = name
        block["mapped_genes"] = len(indices)
        block["mean_delta_log1p_cpm"] = np.asarray(response[:, indices]).mean(axis=1)
        rows.append(block)
    return pd.concat(rows, ignore_index=True)


def _paired_bootstrap_difference(
    values: NDArray[np.floating], seed: int = SEED, n_bootstrap: int = 2_000
) -> tuple[float, float, float]:
    array = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    draws = rng.choice(array, size=(n_bootstrap, len(array)), replace=True).mean(axis=1)
    low, high = np.quantile(draws, [0.025, 0.975])
    return float(array.mean()), float(low), float(high)


def _descriptive_temporal_analysis(root: Path) -> dict[str, object]:
    response_meta, response = read_vector_parquet(
        root / "data" / "processed" / "response_timecourse.parquet", "delta_log1p_cpm"
    )
    genes = pq.read_table(root / "data" / "processed" / "gene_metadata.parquet").to_pandas()
    primary_mask = response_meta["eligible_t10"].astype(bool).to_numpy()
    primary_meta = response_meta.loc[primary_mask].reset_index(drop=True)
    primary = np.asarray(response[primary_mask], dtype=np.float64)
    if len(primary_meta) != 110 or primary_meta["depmap_id"].nunique() != 22:
        raise ValueError("W8 primary cohort must contain 22 lines x 5 time points")

    variances = primary.var(axis=0)
    eligible = np.flatnonzero(variances > 1e-10)
    selected = eligible[np.argsort(-variances[eligible], kind="stable")[:5_000]]
    pca = PCA(n_components=5, svd_solver="randomized", random_state=SEED)
    scores = pca.fit_transform(primary[:, selected])
    score_table = primary_meta[["cell_line", "depmap_id", "disease", "time_hours"]].copy()
    for component in range(scores.shape[1]):
        score_table[f"response_pc{component + 1}"] = scores[:, component]
    score_table.to_csv(root / "results" / "tables" / "temporal_component_scores.csv", index=False)

    loading_rows: list[dict[str, object]] = []
    for component, loadings in enumerate(pca.components_, start=1):
        order = np.argsort(loadings, kind="stable")
        for direction, positions in (("negative", order[:25]), ("positive", order[-25:][::-1])):
            for rank, position in enumerate(positions, start=1):
                gene_index = int(selected[position])
                loading_rows.append(
                    {
                        "component": component,
                        "direction": direction,
                        "rank": rank,
                        "gene_id": genes.loc[gene_index, "gene_id"],
                        "gene_symbol": genes.loc[gene_index, "gene_symbol"],
                        "loading": float(loadings[position]),
                    }
                )
    pd.DataFrame(loading_rows).to_csv(
        root / "results" / "tables" / "temporal_component_loadings.csv", index=False
    )

    mean_by_time = {
        hour: primary[primary_meta["time_hours"].eq(hour).to_numpy()].mean(axis=0)
        for hour in TIMEPOINTS
    }
    heterogeneity_rows: list[dict[str, object]] = []
    for row_index, row in primary_meta.iterrows():
        vector = primary[row_index]
        global_mean = mean_by_time[int(row.time_hours)]
        heterogeneity_rows.append(
            {
                "cell_line": row.cell_line,
                "depmap_id": row.depmap_id,
                "disease": row.disease,
                "time_hours": int(row.time_hours),
                "response_rms": float(np.sqrt(np.mean(np.square(vector)))),
                "rmse_to_time_mean": rmse(vector, global_mean),
            }
        )
    heterogeneity = pd.DataFrame(heterogeneity_rows)
    heterogeneity.to_csv(
        root / "results" / "tables" / "temporal_heterogeneity_by_line.csv", index=False
    )
    heterogeneity_summary_rows: list[dict[str, object]] = []
    for hour in TIMEPOINTS:
        subset = primary[primary_meta["time_hours"].eq(hour).to_numpy()]
        cross_line_gene_sd_rms = float(np.sqrt(np.mean(np.var(subset, axis=0, ddof=1))))
        line_values = heterogeneity.loc[
            heterogeneity["time_hours"].eq(hour), "rmse_to_time_mean"
        ].to_numpy()
        mean, low, high = bootstrap_mean_interval(line_values, seed=SEED + hour)
        heterogeneity_summary_rows.append(
            {
                "time_hours": hour,
                "cross_line_gene_sd_rms": cross_line_gene_sd_rms,
                "mean_line_rmse_to_time_mean": mean,
                "ci95_low": low,
                "ci95_high": high,
            }
        )
    heterogeneity_summary = pd.DataFrame(heterogeneity_summary_rows)
    heterogeneity_summary.to_csv(
        root / "results" / "tables" / "temporal_heterogeneity_summary.csv", index=False
    )

    per_line_period = (
        heterogeneity.assign(
            period=np.where(
                heterogeneity["time_hours"].isin([3, 6]),
                "early",
                np.where(heterogeneity["time_hours"].isin([24, 48]), "late", "middle"),
            )
        )
        .loc[lambda x: x["period"].ne("middle")]
        .groupby(["cell_line", "period"], observed=True)["rmse_to_time_mean"]
        .mean()
        .unstack()
    )
    early_late_difference = per_line_period["late"] - per_line_period["early"]
    early_late_mean, early_late_low, early_late_high = _paired_bootstrap_difference(
        early_late_difference.to_numpy()
    )

    genesets = load_geneset_config(root / "config" / "genesets.yaml")
    pathway = _pathway_scores(genes, primary_meta, primary, genesets)
    pathway.to_csv(root / "results" / "tables" / "timecourse_pathway_scores.csv", index=False)
    pathway_summary_rows: list[dict[str, object]] = []
    for (name, hour), block in pathway.groupby(["pathway", "time_hours"], observed=True):
        mean, low, high = bootstrap_mean_interval(
            block["mean_delta_log1p_cpm"].to_numpy(), seed=SEED + int(hour)
        )
        pathway_summary_rows.append(
            {
                "pathway": name,
                "time_hours": int(hour),
                "mapped_genes": int(block["mapped_genes"].iloc[0]),
                "mean_delta_log1p_cpm": mean,
                "ci95_low": low,
                "ci95_high": high,
                "cross_line_sd": float(block["mean_delta_log1p_cpm"].std(ddof=1)),
            }
        )
    pathway_summary = pd.DataFrame(pathway_summary_rows)
    pathway_summary.to_csv(
        root / "results" / "tables" / "timecourse_pathway_summary.csv", index=False
    )

    return {
        "primary_lines": 22,
        "primary_rows": 110,
        "pca_variable_genes": int(len(selected)),
        "pca_explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
        "early_late_rmse_difference_mean": early_late_mean,
        "early_late_rmse_difference_ci95": [early_late_low, early_late_high],
    }


def _temporal_control_qc(root: Path) -> dict[str, object]:
    metadata, values = read_vector_parquet(
        root / "data" / "processed" / "pseudobulk_timecourse.parquet", "log1p_cpm"
    )
    lines = sorted(metadata.loc[metadata["eligible_t10"], "cell_line"].unique())
    dmso = _aligned_rows(metadata, values, lines, condition="control", time_hours=48)
    untreated = _aligned_rows(metadata, values, lines, condition="untreated", time_hours=48)
    rows = [
        {
            "cell_line": line,
            "rmse_dmso_vs_untreated": rmse(dmso[index], untreated[index]),
            "pcc_dmso_vs_untreated": pearson_or_nan(dmso[index], untreated[index]),
        }
        for index, line in enumerate(lines)
    ]
    table = pd.DataFrame(rows)
    table.to_csv(root / "results" / "tables" / "timecourse_control_48h_qc.csv", index=False)
    return {
        "lines": len(lines),
        "median_rmse": float(table["rmse_dmso_vs_untreated"].median()),
        "median_pcc": float(table["pcc_dmso_vs_untreated"].median()),
        "pooling_performed": False,
    }


def _fixed_transfer_predictions(root: Path) -> dict[str, object]:
    primary_meta, primary = read_vector_parquet(
        root / "data" / "processed" / "pseudobulk_24h.parquet", "log1p_cpm"
    )
    response_meta, train_response = read_vector_parquet(
        root / "data" / "processed" / "response_24h.parquet", "delta_log1p_cpm"
    )
    train_lines = sorted(response_meta["cell_line"].astype(str).tolist())
    train_control = _aligned_rows(
        primary_meta, primary, train_lines, condition="control", time_hours=24
    )
    train_response = _aligned_rows(response_meta, train_response, train_lines)
    train_depmap = set(response_meta["depmap_id"].astype(str))

    temporal_meta, temporal_values = read_vector_parquet(
        root / "data" / "processed" / "pseudobulk_timecourse.parquet", "log1p_cpm"
    )
    temporal_response_meta, temporal_response = read_vector_parquet(
        root / "data" / "processed" / "response_timecourse.parquet", "delta_log1p_cpm"
    )
    external_meta = (
        temporal_response_meta.loc[
            temporal_response_meta["eligible_t10"]
            & ~temporal_response_meta["depmap_id"].astype(str).isin(train_depmap)
        ]
        .drop_duplicates("depmap_id")
        .sort_values("cell_line", ignore_index=True)
    )
    external_lines = external_meta["cell_line"].astype(str).tolist()
    if len(external_lines) < 10:
        raise ValueError("too few genuinely external W8 lines for temporal transfer")

    control_embedding = fit_control_embedding(
        train_control,
        max_components=20,
        max_variable_genes=5_000,
        min_mean_log1p_cpm=0.1,
        seed=SEED,
    )
    train_scores = control_embedding.transform(train_control)[:, :20]
    direct = Ridge(alpha=100.0, fit_intercept=True, solver="svd")
    direct.fit(train_scores, train_response)
    response_embedding = fit_response_embedding(train_response, max_components=20, seed=SEED + 1)
    mean_response = train_response.mean(axis=0).astype(np.float32)

    prediction_meta_rows: list[dict[str, object]] = []
    prediction_vectors: list[NDArray[np.float32]] = []
    metric_rows: list[dict[str, object]] = []
    for hour in TIMEPOINTS:
        test_control = _aligned_rows(
            temporal_meta,
            temporal_values,
            external_lines,
            condition="control",
            time_hours=hour,
        )
        observed = _aligned_rows(
            temporal_response_meta, temporal_response, external_lines, time_hours=hour
        )
        test_scores = control_embedding.transform(test_control)[:, :20]
        model_predictions = {
            "B1": np.repeat(mean_response[None, :], len(external_lines), axis=0),
            "B4_FIXED_D20_A100": np.asarray(direct.predict(test_scores), dtype=np.float32),
            "CCLR_FIXED_D20_R20_A100": predict_fixed_low_rank(
                train_response,
                train_scores,
                test_scores,
                response_embedding,
                response_rank=20,
                alpha=100.0,
            ),
        }
        for model, predictions in model_predictions.items():
            for index, cell_line in enumerate(external_lines):
                row = external_meta.iloc[index]
                prediction_meta_rows.append(
                    {
                        "cell_line": cell_line,
                        "depmap_id": row.depmap_id,
                        "disease": row.disease,
                        "time_hours": hour,
                        "model": model,
                        "training_lines": len(train_lines),
                        "external_test": True,
                        "temporal_response_used_for_fit": False,
                    }
                )
                prediction_vectors.append(predictions[index])
                metric_rows.append(
                    {
                        "cell_line": cell_line,
                        "depmap_id": row.depmap_id,
                        "disease": row.disease,
                        "time_hours": hour,
                        "model": model,
                        "rmse_delta": rmse(observed[index], predictions[index]),
                        "pcc_delta": pearson_or_nan(observed[index], predictions[index]),
                        "spearman_delta": spearman_or_nan(observed[index], predictions[index]),
                        "rmse_gain_vs_b1": rmse(observed[index], mean_response)
                        - rmse(observed[index], predictions[index]),
                        "external_test": True,
                        "temporal_response_used_for_fit": False,
                    }
                )
    prediction_path = root / "results" / "predictions" / "temporal_transfer_predictions.parquet"
    write_vector_parquet(
        prediction_path,
        pd.DataFrame(prediction_meta_rows),
        np.asarray(prediction_vectors),
        "predicted_delta_log1p_cpm",
    )
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(
        root / "results" / "tables" / "temporal_transfer_metrics_by_line.csv",
        index=False,
    )
    summary_rows: list[dict[str, object]] = []
    for (hour, model), block in metrics.groupby(["time_hours", "model"], observed=True):
        for metric in ("rmse_delta", "pcc_delta", "spearman_delta", "rmse_gain_vs_b1"):
            mean, low, high = bootstrap_mean_interval(
                block[metric].to_numpy(), seed=SEED + int(hour) + TRANSFER_MODELS.index(model)
            )
            summary_rows.append(
                {
                    "time_hours": int(hour),
                    "model": model,
                    "metric": metric,
                    "macro_mean": mean,
                    "ci95_low": low,
                    "ci95_high": high,
                    "n_cell_lines": len(block),
                }
            )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(root / "results" / "tables" / "temporal_transfer_summary.csv", index=False)
    best_hour = int(
        summary.loc[
            summary["model"].eq("B4_FIXED_D20_A100")
            & summary["metric"].eq("rmse_delta")
        ].sort_values("macro_mean").iloc[0]["time_hours"]
    )
    return {
        "training_lines": len(train_lines),
        "eligible_temporal_lines": int(
            temporal_response_meta.loc[
                temporal_response_meta["eligible_t10"], "depmap_id"
            ].nunique()
        ),
        "external_test_lines": len(external_lines),
        "overlapping_lines_excluded": int(
            temporal_response_meta.loc[
                temporal_response_meta["eligible_t10"]
                & temporal_response_meta["depmap_id"].astype(str).isin(train_depmap),
                "depmap_id",
            ].nunique()
        ),
        "best_b4_transfer_hour_by_rmse": best_hour,
        "temporal_response_used_for_fit": False,
        "prediction_sha256": _sha256(prediction_path),
    }


def _write_temporal_figures(root: Path) -> None:
    """Write the two W8 chart contracts as static, report-ready figures."""
    scores = pd.read_csv(root / "results" / "tables" / "temporal_component_scores.csv")
    pathways = pd.read_csv(root / "results" / "tables" / "timecourse_pathway_summary.csv")
    heterogeneity = pd.read_csv(root / "results" / "tables" / "temporal_heterogeneity_summary.csv")
    transfer = pd.read_csv(root / "results" / "tables" / "temporal_transfer_summary.csv")
    figure_dir = root / "results" / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    colors = ["#136F63", "#D1495B", "#00798C", "#EDAE49", "#6A4C93"]

    # Contract: discrete-time point/interval panels; x=time, y=mean score; color=component/pathway.
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
    for index in range(3):
        column = f"response_pc{index + 1}"
        means = scores.groupby("time_hours")[column].mean().reindex(TIMEPOINTS)
        errors = scores.groupby("time_hours")[column].sem().reindex(TIMEPOINTS) * 1.96
        axes[0].errorbar(
            TIMEPOINTS,
            means,
            yerr=errors,
            fmt="o",
            capsize=3,
            label=f"PC{index + 1}",
            color=colors[index],
        )
    axes[0].axhline(0, color="#9CA3AF", linewidth=0.8)
    axes[0].set(
        title="Shared response components",
        xlabel="Hours",
        ylabel="Mean PCA score ± 1.96 SE",
    )
    axes[0].legend(frameon=False)
    selected_pathways = [
        "immediate_early_response",
        "MAPK_KRAS_signaling",
        "E2F_targets",
        "G2M_checkpoint",
    ]
    for index, pathway in enumerate(selected_pathways):
        block = (
            pathways.loc[pathways["pathway"].eq(pathway)]
            .set_index("time_hours")
            .reindex(TIMEPOINTS)
        )
        axes[1].errorbar(
            TIMEPOINTS,
            block["mean_delta_log1p_cpm"],
            yerr=np.vstack(
                [
                    block["mean_delta_log1p_cpm"] - block["ci95_low"],
                    block["ci95_high"] - block["mean_delta_log1p_cpm"],
                ]
            ),
            fmt="o",
            capsize=3,
            label=pathway.replace("_", " "),
            color=colors[index],
        )
    axes[1].axhline(0, color="#9CA3AF", linewidth=0.8)
    axes[1].set(
        title="Predefined pathway responses",
        xlabel="Hours",
        ylabel="Mean Δ log1p(CPM), 95% CI",
    )
    axes[1].legend(frameon=False, fontsize=8)
    for axis in axes:
        axis.set_xticks(TIMEPOINTS)
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", color="#E5E7EB", linewidth=0.7)
    fig.suptitle(
        "Trametinib response across five discrete time points",
        x=0.02,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    fig.savefig(figure_dir / "temporal_component_trajectories.png", dpi=180)
    plt.close(fig)

    # Contract: left shows cross-line dispersion; right compares fixed 24h-model transfer.
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
    axes[0].errorbar(
        heterogeneity["time_hours"],
        heterogeneity["mean_line_rmse_to_time_mean"],
        yerr=np.vstack(
            [
                heterogeneity["mean_line_rmse_to_time_mean"] - heterogeneity["ci95_low"],
                heterogeneity["ci95_high"] - heterogeneity["mean_line_rmse_to_time_mean"],
            ]
        ),
        fmt="o",
        capsize=3,
        color=colors[0],
        label="Line RMSE to time mean",
    )
    axes[0].scatter(
        heterogeneity["time_hours"],
        heterogeneity["cross_line_gene_sd_rms"],
        marker="s",
        color=colors[1],
        label="Gene-wise SD RMS",
    )
    axes[0].set(
        title="Cross-line response heterogeneity",
        xlabel="Hours",
        ylabel="Dispersion (Δ log1p(CPM))",
    )
    axes[0].legend(frameon=False, fontsize=8)
    rmse_summary = transfer.loc[transfer["metric"].eq("rmse_delta")]
    for index, model in enumerate(TRANSFER_MODELS):
        block = (
            rmse_summary.loc[rmse_summary["model"].eq(model)]
            .set_index("time_hours")
            .reindex(TIMEPOINTS)
        )
        axes[1].errorbar(
            TIMEPOINTS,
            block["macro_mean"],
            yerr=np.vstack(
                [block["macro_mean"] - block["ci95_low"], block["ci95_high"] - block["macro_mean"]],
            ),
            fmt="o",
            capsize=3,
            color=colors[index],
            label=model,
        )
    axes[1].set(
        title="Frozen 24h model on external lines",
        xlabel="Hours",
        ylabel="Macro RMSE, 95% CI",
    )
    axes[1].legend(frameon=False, fontsize=8)
    for axis in axes:
        axis.set_xticks(TIMEPOINTS)
        axis.yaxis.set_major_locator(MaxNLocator(nbins=6))
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", color="#E5E7EB", linewidth=0.7)
    fig.suptitle(
        "Temporal heterogeneity and model-transfer limits",
        x=0.02,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    fig.savefig(figure_dir / "temporal_heterogeneity_transfer.png", dpi=180)
    plt.close(fig)


def run_temporal(root: Path) -> dict[str, object]:
    """Execute W8 from raw integrated counts through figures and audit log."""
    build = build_timecourse_pseudobulk(root)
    descriptive = _descriptive_temporal_analysis(root)
    control_qc = _temporal_control_qc(root)
    transfer = _fixed_transfer_predictions(root)
    _write_temporal_figures(root)
    report = {
        "stage": "W8_temporal",
        "status": "complete_pending_independent_validation",
        "cohort": build,
        "descriptive": descriptive,
        "control_48h_qc": control_qc,
        "external_transfer": transfer,
        "interpretation_contract": {
            "pca": "whole_timecourse_descriptive_only",
            "primary_delta": "Tram_time_minus_DMSO_time",
            "untreated_48h": "control_QC_only_not_pooled",
            "single_cell_trajectory_claim": False,
        },
    }
    log_path = root / "results" / "logs" / "temporal_summary.json"
    log_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    return report
