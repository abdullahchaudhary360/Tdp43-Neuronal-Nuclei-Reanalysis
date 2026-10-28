#!/usr/bin/env python3
"""Validate the supplied matrix and run donor-aware DESeq2-family models.

Contrast direction is nuclear TDP-43 loss minus nuclear TDP-43 retained.
Thus, negative log2FoldChange values indicate lower expression upon TDP-43 loss.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import pydeseq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

BIOMARKERS = {
    "ENSG00000104435": "STMN2",
    "ENSG00000130477": "UNC13A",
    "ENSG00000120948": "TARDBP",
}
MODELS = {
    "status_only": "~nuclear_status",
    "diagnosis_adjusted": "~diagnosis + nuclear_status",
    "donor_aware": "~donor + nuclear_status",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_inputs(counts_path: Path, summary_path: Path, metadata_path: Path):
    summary = pd.read_csv(summary_path, sep="\t", index_col=0)
    order = [
        col.split("/")[-1].replace("_Aligned.sortedByCoord.out.bam", "")
        for col in summary.columns
    ]
    counts = pd.read_csv(
        counts_path, sep="\t", header=None, names=["gene_id", *order]
    ).set_index("gene_id")
    summary.columns = order
    metadata = pd.read_csv(metadata_path).set_index("Sample").loc[order].copy()
    metadata["nuclear_status"] = metadata["Nuclear_TDP43"].map(
        {"Retained": "Positive", "Loss": "Negative"}
    )
    metadata["diagnosis"] = metadata["Diagnosis"]
    metadata["donor"] = metadata["Donor"]
    if metadata[["nuclear_status", "diagnosis", "donor"]].isna().any().any():
        raise ValueError("Metadata contain missing model variables.")
    return counts, summary, metadata


def validate_inputs(counts: pd.DataFrame, summary: pd.DataFrame) -> dict:
    if counts.index.duplicated().any():
        raise ValueError("Duplicate Ensembl gene identifiers detected.")
    if counts.isna().any().any():
        raise ValueError("Missing values detected in count matrix.")
    if (counts < 0).any().any() or not np.allclose(counts.to_numpy() % 1, 0):
        raise ValueError("Counts must be non-negative integers.")
    library_sums = counts.sum(axis=0)
    assigned = summary.loc["Assigned", counts.columns]
    if not np.array_equal(library_sums.to_numpy(), assigned.to_numpy()):
        raise ValueError("Library sums do not match featureCounts Assigned counts.")
    return {
        "genes": int(counts.shape[0]),
        "libraries": int(counts.shape[1]),
        "all_zero_genes": int((counts.sum(axis=1) == 0).sum()),
        "assigned_min": int(library_sums.min()),
        "assigned_median": float(library_sums.median()),
        "assigned_max": int(library_sums.max()),
        "matrix_matches_featureCounts_Assigned": True,
    }


def fit_model(name: str, design: str, counts: pd.DataFrame, metadata: pd.DataFrame,
              outdir: Path, n_cpus: int) -> dict:
    print(f"Fitting {name}: {design}", flush=True)
    dds = DeseqDataSet(
        counts=counts.T,
        metadata=metadata,
        design=design,
        refit_cooks=True,
        n_cpus=n_cpus,
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(
        dds,
        contrast=["nuclear_status", "Negative", "Positive"],
        alpha=0.05,
        cooks_filter=True,
        independent_filter=True,
        n_cpus=n_cpus,
        quiet=True,
    )
    stats.summary()
    results = stats.results_df.copy().sort_values("padj")
    results.index.name = "gene_id"
    results.to_csv(outdir / f"{name}_results.tsv", sep="\t")
    if name == "donor_aware":
        dds.obs["size_factors"].to_csv(outdir / "size_factors.tsv", sep="\t")
    significant = results[(results["padj"] < 0.05) & results["log2FoldChange"].notna()]
    summary = {
        "design": design,
        "contrast": "TDP43_loss_minus_TDP43_retained",
        "genes_with_adjusted_p": int(results["padj"].notna().sum()),
        "significant_FDR_0_05": int(len(significant)),
        "downregulated_in_TDP43_loss": int((significant["log2FoldChange"] < 0).sum()),
        "upregulated_in_TDP43_loss": int((significant["log2FoldChange"] > 0).sum()),
    }
    for gene_id, symbol in BIOMARKERS.items():
        row = results.loc[gene_id]
        summary[symbol] = {
            "baseMean": float(row["baseMean"]),
            "log2FoldChange": float(row["log2FoldChange"]),
            "lfcSE": float(row["lfcSE"]),
            "pvalue": float(row["pvalue"]),
            "padj": float(row["padj"]),
        }
    (outdir / f"{name}_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--models", nargs="+", choices=list(MODELS), default=list(MODELS))
    parser.add_argument("--min-total-count", type=int, default=10)
    parser.add_argument("--n-cpus", type=int, default=1)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    counts, summary, metadata = load_inputs(args.counts, args.summary, args.metadata)
    report = validate_inputs(counts, summary)
    report.update({
        "donors": int(metadata["donor"].nunique()),
        "complete_pair_donors": int(
            (metadata.groupby("donor")["nuclear_status"].nunique() == 2).sum()
        ),
        "minimum_total_count_per_gene": args.min_total_count,
        "pydeseq2_version": pydeseq2.__version__,
        "python_version": platform.python_version(),
        "input_sha256": {
            args.counts.name: sha256(args.counts),
            args.summary.name: sha256(args.summary),
            args.metadata.name: sha256(args.metadata),
        },
    })
    filtered = counts.loc[counts.sum(axis=1) >= args.min_total_count]
    report["genes_after_filter"] = int(filtered.shape[0])
    (args.outdir / "data_validation.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    if args.validate_only:
        return

    model_summaries = {}
    for name in args.models:
        model_summaries[name] = fit_model(
            name, MODELS[name], filtered, metadata, args.outdir, args.n_cpus
        )
    (args.outdir / "model_comparison.json").write_text(
        json.dumps(model_summaries, indent=2)
    )


if __name__ == "__main__":
    main()
