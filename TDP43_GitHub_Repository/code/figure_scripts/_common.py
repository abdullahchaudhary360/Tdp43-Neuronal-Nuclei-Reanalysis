"""Shared loader for the one-script-per-figure entry points."""
from __future__ import annotations
import argparse
import importlib.util
from pathlib import Path


def load_figure_module():
    module_path = Path(__file__).resolve().parents[1] / '04_generate_figures.py'
    spec = importlib.util.spec_from_file_location('tdp43_figures', module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run(figure_id: str) -> None:
    parser = argparse.ArgumentParser(description=f'Regenerate {figure_id} from the verified package inputs.')
    parser.add_argument('--package-dir', type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    root = args.package_dir.resolve()
    fg = load_figure_module()
    fg.setup_style()
    counts, summary, meta, results, annotated = fg.load_inputs(root/'data', root/'metadata', root/'tables')
    _, normalized, _ = fg.median_ratio_normalize(counts)
    summaries = fg.model_summary(results)
    outdir = root/'figures'
    tables = root/'tables'
    outdir.mkdir(exist_ok=True)
    actions = {
        'Figure_1_study_design': lambda: fg.figure1(outdir, meta),
        'Figure_2_QC_and_PCA': lambda: fg.figure2(outdir, summary, meta, normalized, tables),
        'Figure_3_differential_expression': lambda: fg.figure3(outdir, results, annotated, summaries, tables),
        'Figure_4_marker_expression': lambda: fg.figure4(outdir, normalized, meta, results, tables),
        'Figure_5_leading_gene_heatmap': lambda: fg.figure5(outdir, normalized, meta, annotated, tables),
        'Figure_S1_featureCounts_assignment': lambda: fg.figure_s1(outdir, summary, meta),
        'Figure_S2_sensitivity_volcanoes': lambda: fg.figure_s2(outdir, results),
        'Figure_S3_cryptic_exon_evidence_boundary': lambda: fg.figure_s3(outdir),
    }
    if figure_id not in actions:
        raise ValueError(f'Unknown figure: {figure_id}')
    actions[figure_id]()
    print(f'Regenerated {figure_id} in {outdir}')
