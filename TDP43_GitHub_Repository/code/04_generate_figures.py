#!/usr/bin/env python3
"""Reproduce all manuscript figures from the supplied count matrix and results.

The contrast is nuclear TDP-43 loss minus nuclear TDP-43 retained. Negative
log2 fold changes therefore indicate lower expression in TDP-43-loss nuclei.
This script performs descriptive normalization/PCA and visualization only; it
does not estimate cryptic-exon inclusion, PSI, or junction-level splicing.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA
from PIL import Image

MARKERS = {
    "ENSG00000104435": "STMN2",
    "ENSG00000130477": "UNC13A",
    "ENSG00000120948": "TARDBP",
}
MODEL_LABELS = {
    "status_only": "Status only",
    "diagnosis_adjusted": "Diagnosis adjusted",
    "donor_aware": "Donor aware",
}
STATUS_COLORS = {"Retained": "#2878B5", "Loss": "#C44E52"}
DIAG_MARKERS = {"FTD": "o", "FTD-ALS": "s"}
UP = "#C44E52"
DOWN = "#2878B5"
GREY = "#A7ADB4"
DARK = "#263238"
ACCENT = "#4C78A8"


def setup_style() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "axes.linewidth": 0.7,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "figure.dpi": 120,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def panel_label(ax, label: str) -> None:
    ax.text(-0.14, 1.06, label, transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="top", ha="left", color=DARK)


def save_figure(fig: plt.Figure, outbase: Path, tiff: bool = True) -> None:
    outbase.parent.mkdir(parents=True, exist_ok=True)
    png = outbase.with_suffix(".png")
    fig.savefig(png, dpi=400, facecolor="white")
    fig.savefig(outbase.with_suffix(".pdf"), facecolor="white")
    if tiff:
        with Image.open(png) as im:
            im.convert("RGB").save(outbase.with_suffix(".tiff"), compression="tiff_lzw", dpi=(400, 400))
    plt.close(fig)


def load_inputs(data_dir: Path, metadata_dir: Path, tables_dir: Path):
    summary = pd.read_csv(data_dir / "raw_counts.txt.summary", sep="\t", index_col=0)
    order = [c.split("/")[-1].replace("_Aligned.sortedByCoord.out.bam", "") for c in summary.columns]
    summary.columns = order
    counts = pd.read_csv(data_dir / "counts_matrix.txt", sep="\t", header=None,
                         names=["gene_id", *order]).set_index("gene_id")
    meta = pd.read_csv(metadata_dir / "metadata_verified.csv").set_index("Sample").loc[order].copy()
    results = {}
    for name in MODEL_LABELS:
        # Supplied results are CSV; a fresh run of 02 writes TSV. Accept both
        # so the complete validation → model → annotation → figure workflow
        # can be reproduced without changing file contents by hand.
        csv_path = tables_dir / f"{name}_results.csv"
        tsv_path = tables_dir / f"{name}_results.tsv"
        if csv_path.exists():
            results[name] = pd.read_csv(csv_path)
        elif tsv_path.exists():
            results[name] = pd.read_csv(tsv_path, sep="\t")
        else:
            raise FileNotFoundError(f"Missing {name} results: expected {csv_path.name} or {tsv_path.name}")
    annotated = pd.read_csv(tables_dir / "donor_aware_results_annotated.tsv", sep="\t")
    return counts, summary, meta, results, annotated


def median_ratio_normalize(counts: pd.DataFrame):
    filtered = counts.loc[counts.sum(axis=1) >= 10].copy()
    # Standard median-of-ratios size factors: use genes with nonzero counts
    # in every library so each geometric mean is finite and comparable.
    with np.errstate(divide="ignore"):
        logged = np.log(filtered.astype(float))
    geometric = logged.mean(axis=1)
    valid = np.isfinite(geometric)
    size_factors = np.exp(logged.loc[valid].sub(geometric.loc[valid], axis=0).median(axis=0))
    normalized = filtered.div(size_factors, axis=1)
    return filtered, normalized, size_factors


def model_summary(results: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for model, df in results.items():
        sig = df[(df.padj < 0.05) & df.log2FoldChange.notna()]
        row = {
            "model": model,
            "design": {
                "status_only": "~ nuclear_status",
                "diagnosis_adjusted": "~ diagnosis + nuclear_status",
                "donor_aware": "~ donor + nuclear_status",
            }[model],
            "genes_with_adjusted_p": int(df.padj.notna().sum()),
            "significant_FDR_0_05": int(len(sig)),
            "lower_in_TDP43_loss": int((sig.log2FoldChange < 0).sum()),
            "higher_in_TDP43_loss": int((sig.log2FoldChange > 0).sum()),
            "FDR_and_abs_LFC_at_least_1": int((sig.log2FoldChange.abs() >= 1).sum()),
        }
        for gid, symbol in MARKERS.items():
            r = df.loc[df.gene_id == gid].iloc[0]
            row[f"{symbol}_log2FC"] = r.log2FoldChange
            row[f"{symbol}_lfcSE"] = r.lfcSE
            row[f"{symbol}_padj"] = r.padj
        rows.append(row)
    return pd.DataFrame(rows)


def figure1(outdir: Path, meta: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(10.5, 5.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.05, 1], hspace=0.5, wspace=0.35)
    ax = fig.add_subplot(gs[0, :])
    ax.axis("off")
    steps = [
        ("GSE126542", "Frontal-cortex nuclei\npublic RNA-seq"),
        ("Cohort audit", "NeuN+ libraries\nTDP-43 labels verified"),
        ("Count validation", "62,754 genes\n12 libraries"),
        ("Primary model", "~ donor + status\n5 complete pairs"),
        ("Interpretation", "Gene-level replication\nSplicing not measured"),
    ]
    xs = np.linspace(0.02, 0.82, len(steps))
    for i, ((title, subtitle), x) in enumerate(zip(steps, xs)):
        box = FancyBboxPatch((x, 0.22), 0.15, 0.56, boxstyle="round,pad=0.012,rounding_size=0.015",
                             transform=ax.transAxes, facecolor="#F3F6F8", edgecolor="#78909C", linewidth=0.8)
        ax.add_patch(box)
        ax.text(x + 0.075, 0.62, title, transform=ax.transAxes, ha="center", va="center",
                fontsize=9, fontweight="bold", color=DARK)
        ax.text(x + 0.075, 0.40, subtitle, transform=ax.transAxes, ha="center", va="center",
                fontsize=7.2, color="#455A64")
        if i < len(steps) - 1:
            ax.add_patch(FancyArrowPatch((x + 0.153, 0.50), (xs[i + 1] - 0.005, 0.50),
                                         transform=ax.transAxes, arrowstyle="-|>", mutation_scale=10,
                                         linewidth=0.8, color="#78909C"))
    ax.text(0.0, 1.02, "A", transform=ax.transAxes, fontsize=12, fontweight="bold", color=DARK)

    axb = fig.add_subplot(gs[1, 0:2])
    axb.axis("off")
    donors = ["S1", "S2", "S3", "S4", "S5", "S6", "S7"]
    ypos = {d: len(donors) - 1 - i for i, d in enumerate(donors)}
    for donor in donors:
        sub = meta[meta.Donor == donor]
        y = ypos[donor]
        if len(sub) == 2:
            axb.plot([0.28, 0.72], [y, y], color="#B0BEC5", lw=1.0, zorder=1)
        for status, x in [("Retained", 0.28), ("Loss", 0.72)]:
            ss = sub[sub.Nuclear_TDP43 == status]
            if len(ss):
                diag = ss.iloc[0].Diagnosis
                axb.scatter(x, y, s=90, marker=DIAG_MARKERS[diag], color=STATUS_COLORS[status],
                            edgecolor="white", linewidth=0.6, zorder=3)
            else:
                axb.scatter(x, y, s=65, marker="o", facecolor="white", edgecolor="#CFD8DC", linewidth=0.8)
        axb.text(0.05, y, donor, va="center", ha="left", fontsize=8, fontweight="bold")
        axb.text(0.93, y, sub.iloc[0].Diagnosis, va="center", ha="left", fontsize=7, color="#546E7A")
    axb.text(0.28, len(donors) - 0.25, "TDP-43 retained", ha="center", fontsize=8, color=DOWN, fontweight="bold")
    axb.text(0.72, len(donors) - 0.25, "TDP-43 loss", ha="center", fontsize=8, color=UP, fontweight="bold")
    axb.set_xlim(0, 1.12)
    axb.set_ylim(-0.8, len(donors))
    axb.set_title("Seven donors; five complete matched pairs", loc="left", pad=12, fontweight="bold")
    axb.text(-0.08, 1.08, "B", transform=axb.transAxes, fontsize=12, fontweight="bold", color=DARK)

    axc = fig.add_subplot(gs[1, 2])
    axc.axis("off")
    metrics = [
        ("12", "NeuN+ libraries"),
        ("6 + 6", "retained + loss"),
        ("44,225", "genes after filtering"),
        ("3", "statistical models"),
    ]
    for i, (value, label) in enumerate(metrics):
        y = 0.88 - i * 0.25
        axc.text(0.02, y, value, fontsize=15, fontweight="bold", color=ACCENT, va="center")
        axc.text(0.02, y - 0.09, label, fontsize=7.5, color="#546E7A", va="center")
    axc.text(-0.12, 1.08, "C", transform=axc.transAxes, fontsize=12, fontweight="bold", color=DARK)
    save_figure(fig, outdir / "Figure_1_study_design")


def figure2(outdir: Path, summary: pd.DataFrame, meta: pd.DataFrame,
            normalized: pd.DataFrame, tables_out: Path) -> None:
    order = meta.index.tolist()
    lognorm = np.log2(normalized + 1)
    variable = lognorm.var(axis=1).nlargest(2000).index
    x = lognorm.loc[variable].T
    x = x - x.mean(axis=0)
    pca = PCA(n_components=2, random_state=0)
    coords = pca.fit_transform(x)
    pcdf = pd.DataFrame({"Sample": order, "PC1": coords[:, 0], "PC2": coords[:, 1],
                         "Donor": meta.Donor.values, "Diagnosis": meta.Diagnosis.values,
                         "Nuclear_TDP43": meta.Nuclear_TDP43.values})
    pcdf.to_csv(tables_out / "Table_S6_PCA_coordinates.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5), gridspec_kw={"width_ratios": [1.05, 1.25]})
    ax = axes[0]
    assigned = summary.loc["Assigned", order] / 1e6
    colors = [STATUS_COLORS[s] for s in meta.Nuclear_TDP43]
    bars = ax.bar(np.arange(len(order)), assigned, color=colors, width=0.72)
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels([s.replace("SRR857", "…") for s in order], rotation=45, ha="right")
    ax.set_ylabel("featureCounts-assigned reads (millions)")
    ax.set_xlabel("Library (SRR857…)")
    ax.set_title("Count depth was consistent across libraries", loc="left", fontweight="bold")
    ax.axhline(assigned.median(), color="#607D8B", lw=0.8, ls="--")
    ax.text(len(order) - 0.4, assigned.median() + 0.5, f"median {assigned.median():.1f}M",
            ha="right", fontsize=7, color="#546E7A")
    ax.set_ylim(0, max(assigned) * 1.16)
    panel_label(ax, "A")

    ax = axes[1]
    for donor, sub in pcdf.groupby("Donor"):
        if len(sub) == 2:
            ax.plot(sub.PC1, sub.PC2, color="#B0BEC5", lw=0.8, zorder=1)
    for status in ["Retained", "Loss"]:
        for diag in ["FTD", "FTD-ALS"]:
            sub = pcdf[(pcdf.Nuclear_TDP43 == status) & (pcdf.Diagnosis == diag)]
            ax.scatter(sub.PC1, sub.PC2, s=58, marker=DIAG_MARKERS[diag],
                       color=STATUS_COLORS[status], edgecolor="white", linewidth=0.5,
                       label=f"{status}, {diag}", zorder=3)
    for _, r in pcdf.iterrows():
        ax.text(r.PC1, r.PC2 + 2.0, r.Donor, fontsize=6.5, ha="center", color="#455A64")
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}% variance)")
    ax.set_title("PCA of 2,000 most variable genes", loc="left", fontweight="bold")
    ax.legend(frameon=False, ncol=2, loc="best", handletextpad=0.2, columnspacing=0.8)
    ax.text(0.01, 0.01, "Lines connect matched donors", transform=ax.transAxes, fontsize=6.5, color="#607D8B")
    panel_label(ax, "B")
    fig.subplots_adjust(wspace=0.32)
    save_figure(fig, outdir / "Figure_2_QC_and_PCA")


def volcano(ax, df: pd.DataFrame, annotated: pd.DataFrame | None = None, title: str = "") -> None:
    d = df.dropna(subset=["log2FoldChange", "padj"]).copy()
    d["neglog10"] = -np.log10(d.padj.clip(lower=1e-300))
    sig = d.padj < 0.05
    strong = sig & (d.log2FoldChange.abs() >= 1)
    ax.scatter(d.loc[~sig, "log2FoldChange"], d.loc[~sig, "neglog10"], s=2.5, c=GREY, alpha=0.28, rasterized=True)
    ax.scatter(d.loc[sig & ~strong, "log2FoldChange"], d.loc[sig & ~strong, "neglog10"], s=3,
               c="#D4A72C", alpha=0.45, rasterized=True)
    ax.scatter(d.loc[strong & (d.log2FoldChange < 0), "log2FoldChange"],
               d.loc[strong & (d.log2FoldChange < 0), "neglog10"], s=4, c=DOWN, alpha=0.55, rasterized=True)
    ax.scatter(d.loc[strong & (d.log2FoldChange > 0), "log2FoldChange"],
               d.loc[strong & (d.log2FoldChange > 0), "neglog10"], s=4, c=UP, alpha=0.55, rasterized=True)
    ax.axhline(-math.log10(0.05), color="#78909C", lw=0.7, ls="--")
    ax.axvline(-1, color="#B0BEC5", lw=0.6, ls=":")
    ax.axvline(1, color="#B0BEC5", lw=0.6, ls=":")
    ax.set_xlabel("log2 fold change (TDP-43 loss / retained)")
    ax.set_ylabel("−log10 adjusted P")
    ax.set_title(title, loc="left", fontweight="bold")
    if annotated is not None:
        ann = annotated.drop_duplicates("gene_id").set_index("gene_id")
        labels = set()
        eligible = annotated[(annotated.padj < 0.05) & (annotated.log2FoldChange.abs() >= 1) &
                             (annotated.gene_type == "protein-coding")].nsmallest(6, "padj")
        labels.update(eligible.gene_id)
        labels.update(MARKERS)
        offsets = {}
        for i, gid in enumerate(labels):
            if gid not in ann.index:
                continue
            r = ann.loc[gid]
            if not isinstance(r, pd.Series):
                r = r.iloc[0]
            symbol = r.gene_symbol if pd.notna(r.gene_symbol) else MARKERS.get(gid, gid)
            x, y = r.log2FoldChange, -math.log10(max(r.padj, 1e-300))
            dx = 10 if x >= 0 else -10
            dy = offsets.get(gid, 8 + (i % 3) * 7)
            ax.annotate(symbol, (x, y), xytext=(dx, dy), textcoords="offset points",
                        ha="left" if dx > 0 else "right", fontsize=6.7, fontweight="bold" if gid in MARKERS else "normal",
                        arrowprops=dict(arrowstyle="-", lw=0.45, color="#78909C"), color=DARK)


def figure3(outdir: Path, results: dict[str, pd.DataFrame], annotated: pd.DataFrame,
            summaries: pd.DataFrame, tables_out: Path) -> None:
    fig = plt.figure(figsize=(11, 8.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 0.85], hspace=0.34, wspace=0.28)
    ax = fig.add_subplot(gs[0, :])
    volcano(ax, results["donor_aware"], annotated,
            "Donor-aware differential expression: 3,442 genes at FDR < 0.05")
    panel_label(ax, "A")

    ax = fig.add_subplot(gs[1, 0])
    ss = summaries.set_index("model").loc[list(MODEL_LABELS)]
    x = np.arange(len(ss))
    ax.bar(x, ss.lower_in_TDP43_loss, color=DOWN, width=0.62, label="Lower in TDP-43 loss")
    ax.bar(x, ss.higher_in_TDP43_loss, bottom=ss.lower_in_TDP43_loss, color=UP, width=0.62,
           label="Higher in TDP-43 loss")
    for i, row in enumerate(ss.itertuples()):
        ax.text(i, row.significant_FDR_0_05 + max(ss.significant_FDR_0_05) * 0.02,
                f"{row.significant_FDR_0_05:,}", ha="center", fontsize=7, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([MODEL_LABELS[m] for m in ss.index], rotation=12, ha="right")
    ax.set_ylabel("Genes at FDR < 0.05")
    ax.set_ylim(0, max(ss.significant_FDR_0_05) * 1.18)
    ax.set_title("Model-dependent global significance", loc="left", fontweight="bold")
    ax.legend(frameon=False, loc="upper left")
    panel_label(ax, "B")

    ax = fig.add_subplot(gs[1, 1])
    a = results["donor_aware"][["gene_id", "log2FoldChange"]].rename(columns={"log2FoldChange": "donor"})
    b = results["diagnosis_adjusted"][["gene_id", "log2FoldChange"]].rename(columns={"log2FoldChange": "diagnosis"})
    merged = a.merge(b, on="gene_id").replace([np.inf, -np.inf], np.nan).dropna()
    r, _ = pearsonr(merged.donor, merged.diagnosis)
    rho, _ = spearmanr(merged.donor, merged.diagnosis)
    ax.scatter(merged.diagnosis, merged.donor, s=2.2, color="#78909C", alpha=0.18, rasterized=True)
    lim = max(abs(merged.diagnosis).max(), abs(merged.donor).max()) * 1.03
    ax.plot([-lim, lim], [-lim, lim], color=DARK, lw=0.7, ls="--")
    ax.axhline(0, color="#CFD8DC", lw=0.5)
    ax.axvline(0, color="#CFD8DC", lw=0.5)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("log2 fold change: diagnosis adjusted")
    ax.set_ylabel("log2 fold change: donor aware")
    ax.set_title("Effect sizes were highly concordant", loc="left", fontweight="bold")
    ax.text(0.03, 0.96, f"Pearson r = {r:.3f}\nSpearman ρ = {rho:.3f}", transform=ax.transAxes,
            va="top", fontsize=7.5, color=DARK)
    merged.to_csv(tables_out / "Table_S7_effect_size_concordance.tsv", sep="\t", index=False)
    panel_label(ax, "C")
    save_figure(fig, outdir / "Figure_3_differential_expression")


def figure4(outdir: Path, normalized: pd.DataFrame, meta: pd.DataFrame,
            results: dict[str, pd.DataFrame], tables_out: Path) -> None:
    fig = plt.figure(figsize=(11, 7.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.05, 0.9], hspace=0.45, wspace=0.38)
    donor_counts = meta.Donor.value_counts()
    paired_donors = donor_counts[donor_counts == 2].index.tolist()
    marker_rows = []
    for i, (gid, symbol) in enumerate(MARKERS.items()):
        ax = fig.add_subplot(gs[0, i])
        tab = pd.DataFrame({"Sample": meta.index, "Donor": meta.Donor.values,
                            "Diagnosis": meta.Diagnosis.values, "Nuclear_TDP43": meta.Nuclear_TDP43.values,
                            "Normalized_count": normalized.loc[gid, meta.index].values,
                            "Gene": symbol, "Ensembl_ID": gid})
        marker_rows.append(tab)
        for donor in paired_donors:
            sub = tab[tab.Donor == donor].set_index("Nuclear_TDP43")
            ax.plot([0, 1], [sub.loc["Retained", "Normalized_count"], sub.loc["Loss", "Normalized_count"]],
                    color="#B0BEC5", lw=0.9, zorder=1)
        for status, x in [("Retained", 0), ("Loss", 1)]:
            for diag in ["FTD", "FTD-ALS"]:
                sub = tab[(tab.Nuclear_TDP43 == status) & (tab.Diagnosis == diag)]
                ax.scatter(np.full(len(sub), x), sub.Normalized_count, s=44, marker=DIAG_MARKERS[diag],
                           color=STATUS_COLORS[status], edgecolor="white", linewidth=0.5, zorder=3)
        r = results["donor_aware"].loc[results["donor_aware"].gene_id == gid].iloc[0]
        ax.set_title(f"{symbol}\nLFC = {r.log2FoldChange:.2f}; FDR = {r.padj:.2e}", fontweight="bold", pad=8)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Retained", "Loss"])
        ax.set_ylabel("Normalized count" if i == 0 else "")
        ax.set_xlim(-0.3, 1.3)
        ax.set_ylim(bottom=0)
        panel_label(ax, chr(ord("A") + i))

    ax = fig.add_subplot(gs[1, :])
    y_positions = []
    labels = []
    y = 0
    for symbol in MARKERS.values():
        for model in ["status_only", "diagnosis_adjusted", "donor_aware"]:
            r = results[model].loc[results[model].gene_id == next(k for k, v in MARKERS.items() if v == symbol)].iloc[0]
            y_positions.append(y)
            labels.append(f"{symbol} | {MODEL_LABELS[model]}")
            color = DARK if model == "donor_aware" else "#78909C"
            ax.errorbar(r.log2FoldChange, y, xerr=1.96 * r.lfcSE, fmt="o", ms=5 if model == "donor_aware" else 4,
                        color=color, ecolor=color, elinewidth=1.0, capsize=2)
            y += 1
        y += 0.65
    ax.axvline(0, color="#455A64", lw=0.7)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=6.8)
    ax.invert_yaxis()
    ax.set_xlabel("log2 fold change (TDP-43 loss / retained); mean ± 95% Wald interval")
    ax.set_title("Prespecified marker effects across models", loc="left", fontweight="bold")
    ax.grid(axis="x", color="#ECEFF1", lw=0.6)
    panel_label(ax, "D")
    pd.concat(marker_rows).to_csv(tables_out / "Table_S8_marker_normalized_counts.tsv", sep="\t", index=False)
    save_figure(fig, outdir / "Figure_4_marker_expression")


def figure5(outdir: Path, normalized: pd.DataFrame, meta: pd.DataFrame, annotated: pd.DataFrame,
            tables_out: Path) -> None:
    eligible = annotated[(annotated.padj < 0.05) & (annotated.log2FoldChange.abs() >= 1) &
                         (annotated.baseMean >= 100) & (annotated.gene_type == "protein-coding") &
                         annotated.gene_symbol.notna()].nsmallest(30, "padj").copy()
    genes = eligible.gene_id.tolist()
    lognorm = np.log2(normalized.loc[genes, meta.index] + 1)
    z = lognorm.sub(lognorm.mean(axis=1), axis=0).div(lognorm.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    row_link = linkage(pdist(z.values, metric="correlation"), method="average")
    dend = dendrogram(row_link, no_plot=True)
    order_rows = dend["leaves"]
    z = z.iloc[order_rows]
    labels = eligible.set_index("gene_id").loc[z.index, "gene_symbol"].tolist()
    eligible[["gene_id", "gene_symbol", "baseMean", "log2FoldChange", "padj"]].to_csv(
        tables_out / "Table_S9_heatmap_gene_selection.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(10.5, 8.5))
    gs = fig.add_gridspec(2, 2, width_ratios=[0.15, 1], height_ratios=[0.07, 1],
                          wspace=0.03, hspace=0.05)
    axd = fig.add_subplot(gs[1, 0])
    dendrogram(row_link, orientation="left", ax=axd, color_threshold=None, above_threshold_color="#78909C")
    axd.axis("off")
    ax = fig.add_subplot(gs[1, 1])
    cmap = LinearSegmentedColormap.from_list("expression", ["#2166AC", "#F7F7F7", "#B2182B"])
    im = ax.imshow(z.values, aspect="auto", cmap=cmap, vmin=-2.5, vmax=2.5, interpolation="nearest")
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels, fontsize=6.5, fontstyle="italic")
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_xticks(np.arange(len(meta)))
    ax.set_xticklabels(meta.Donor + "\n" + meta.Diagnosis, fontsize=6.3)
    ax.set_xlabel("Libraries ordered by donor and nuclear TDP-43 status")
    ax.set_title("Leading protein-coding differential genes", loc="left", fontweight="bold", pad=12)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    axa = fig.add_subplot(gs[0, 1])
    axa.axis("off")
    for i, (_, r) in enumerate(meta.iterrows()):
        axa.add_patch(plt.Rectangle((i, 0.48), 1, 0.45, color=STATUS_COLORS[r.Nuclear_TDP43], linewidth=0))
        axa.add_patch(plt.Rectangle((i, 0.02), 1, 0.38, color="#90A4AE" if r.Diagnosis == "FTD" else "#37474F", linewidth=0))
    axa.set_xlim(0, len(meta))
    axa.set_ylim(0, 1)
    axa.text(-0.5, 0.70, "TDP-43", ha="right", va="center", fontsize=6.5)
    axa.text(-0.5, 0.20, "Diagnosis", ha="right", va="center", fontsize=6.5)
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.17)
    cbar.set_label("Row Z score", fontsize=7)
    cbar.ax.tick_params(labelsize=6)
    fig.text(0.02, 0.98, "A", fontsize=12, fontweight="bold", color=DARK, va="top")
    save_figure(fig, outdir / "Figure_5_leading_gene_heatmap")


def figure_s1(outdir: Path, summary: pd.DataFrame, meta: pd.DataFrame) -> None:
    """Generate the complete featureCounts assignment audit."""
    order = meta.index.tolist()
    categories = ["Assigned", "Unassigned_MultiMapping", "Unassigned_NoFeatures", "Unassigned_Ambiguity"]
    labels = ["Assigned", "Multimapping", "No feature", "Ambiguous"]
    colors = ["#4C78A8", "#F2A65A", "#B0BEC5", "#8E6C8A"]
    vals = summary.loc[categories, order].astype(float)
    frac = vals.div(vals.sum(axis=0), axis=1)
    fig, ax = plt.subplots(figsize=(10.5, 4.5))
    bottom = np.zeros(len(order))
    for cat, label, color in zip(categories, labels, colors):
        ax.bar(np.arange(len(order)), frac.loc[cat] * 100, bottom=bottom, label=label, color=color, width=0.72)
        bottom += frac.loc[cat].values * 100
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels([s.replace("SRR857", "…") for s in order], rotation=45, ha="right")
    ax.set_ylabel("Reads in featureCounts summary (%)")
    ax.set_xlabel("Library (SRR857…)")
    ax.set_ylim(0, 100)
    ax.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.12))
    ax.set_title("featureCounts assignment categories; these are not STAR alignment rates", loc="left", fontweight="bold")
    save_figure(fig, outdir / "Figure_S1_featureCounts_assignment")


def figure_s2(outdir: Path, results: dict[str, pd.DataFrame]) -> None:
    """Generate sensitivity-model volcano plots."""
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5))
    for ax, model, label in zip(axes, ["status_only", "diagnosis_adjusted"], ["A", "B"]):
        volcano(ax, results[model], None, MODEL_LABELS[model])
        panel_label(ax, label)
    fig.subplots_adjust(wspace=0.28)
    save_figure(fig, outdir / "Figure_S2_sensitivity_volcanoes")


def figure_s3(outdir: Path) -> None:
    """Generate the explicit evidence boundary for cryptic-exon interpretation."""
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    ax.axis("off")
    boxes = [
        (0.02, "Nuclear TDP-43 loss", "Measured source label\nNeuN+ sorted nuclei"),
        (0.27, "Whole-gene expression", "Measured here\nfeatureCounts + PyDESeq2"),
        (0.52, "Cryptic-exon inclusion / PSI", "Not measured here\nRequires junction-level evidence"),
        (0.77, "RNA or protein consequence", "Not measured here\nRequires isoform/protein assays"),
    ]
    for i, (x, title, subtitle) in enumerate(boxes):
        measured = i < 2
        box = FancyBboxPatch((x, 0.40), 0.20, 0.35, boxstyle="round,pad=0.015,rounding_size=0.02",
                             transform=ax.transAxes, facecolor="#EDF4F8" if measured else "#FBEFEF",
                             edgecolor="#4C78A8" if measured else "#C44E52", linewidth=1.0)
        ax.add_patch(box)
        ax.text(x + 0.10, 0.64, title, transform=ax.transAxes, ha="center", va="center",
                fontsize=8.2, fontweight="bold", color=DARK)
        ax.text(x + 0.10, 0.49, subtitle, transform=ax.transAxes, ha="center", va="center",
                fontsize=7, color="#546E7A")
        if i < 3:
            ax.add_patch(FancyArrowPatch((x + 0.205, 0.575), (boxes[i + 1][0] - 0.006, 0.575),
                                         transform=ax.transAxes, arrowstyle="-|>", mutation_scale=11,
                                         color="#78909C", linewidth=0.9))
    ax.text(0.02, 0.20, "Recommended follow-up", transform=ax.transAxes, fontsize=8.5, fontweight="bold", color=DARK)
    ax.text(0.02, 0.11,
            "Reprocess all eligible sorted libraries and quantify splice junctions with an annotation-aware method.\nValidate UNC13A and STMN2 events by RT–PCR or long-read RNA analysis; test protein-level consequences independently.",
            transform=ax.transAxes, fontsize=7.5, color="#455A64", va="top")
    ax.set_title("Mechanistic model and evidence boundary", loc="left", fontweight="bold", pad=12)
    save_figure(fig, outdir / "Figure_S3_cryptic_exon_evidence_boundary")


def supplementary_figures(outdir: Path, summary: pd.DataFrame, meta: pd.DataFrame,
                          results: dict[str, pd.DataFrame]) -> None:
    figure_s1(outdir, summary, meta)
    figure_s2(outdir, results)
    figure_s3(outdir)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-dir", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.package_dir.resolve()
    data_dir, metadata_dir, tables_dir = root / "data", root / "metadata", root / "tables"
    figures_dir = root / "figures"
    figures_dir.mkdir(exist_ok=True)
    tables_dir.mkdir(exist_ok=True)
    setup_style()

    counts, summary, meta, results, annotated = load_inputs(data_dir, metadata_dir, tables_dir)
    _, normalized, size_factors = median_ratio_normalize(counts)
    size_factors.rename("size_factor").to_csv(tables_dir / "Table_S10_size_factors.tsv", sep="\t")
    summaries = model_summary(results)
    summaries.to_csv(tables_dir / "Table_S5_model_comparison.tsv", sep="\t", index=False)

    figure1(figures_dir, meta)
    figure2(figures_dir, summary, meta, normalized, tables_dir)
    figure3(figures_dir, results, annotated, summaries, tables_dir)
    figure4(figures_dir, normalized, meta, results, tables_dir)
    figure5(figures_dir, normalized, meta, annotated, tables_dir)
    supplementary_figures(figures_dir, summary, meta, results)

    stmn2 = results["donor_aware"].loc[results["donor_aware"].gene_id == "ENSG00000104435"].iloc[0]
    unc13a = results["donor_aware"].loc[results["donor_aware"].gene_id == "ENSG00000130477"].iloc[0]
    report = {
        "contrast": "nuclear_TDP43_loss_minus_retained",
        "libraries": int(len(meta)),
        "donors": int(meta.Donor.nunique()),
        "complete_pairs": int((meta.groupby("Donor").Nuclear_TDP43.nunique() == 2).sum()),
        "STMN2_percent_lower_in_loss": float((1 - 2 ** stmn2.log2FoldChange) * 100),
        "UNC13A_percent_change_in_loss": float((2 ** unc13a.log2FoldChange - 1) * 100),
        "figures_generated": 8,
        "splicing_or_PSI_estimated": False,
    }
    (tables_dir / "figure_generation_summary.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
