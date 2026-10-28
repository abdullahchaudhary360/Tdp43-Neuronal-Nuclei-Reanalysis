#!/usr/bin/env python3
"""Annotate Ensembl gene identifiers through the MyGene.info query service.

The script annotates all FDR-significant genes from the donor-aware model,
all three prespecified marker genes, and the leading genes by adjusted P value.
It preserves unmapped identifiers rather than silently dropping them.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import pandas as pd
import requests

BIOMARKERS = {
    "ENSG00000104435": "STMN2",
    "ENSG00000130477": "UNC13A",
    "ENSG00000120948": "TARDBP",
}


def query_mygene(gene_ids: list[str], batch_size: int = 500) -> pd.DataFrame:
    rows: list[dict] = []
    endpoint = "https://mygene.info/v3/query"
    for start in range(0, len(gene_ids), batch_size):
        batch = gene_ids[start:start + batch_size]
        payload = {
            "q": ",".join(batch),
            "scopes": "ensembl.gene",
            "fields": "symbol,name,type_of_gene",
            "species": "human",
            "size": 1,
        }
        last_error = None
        for attempt in range(3):
            try:
                response = requests.post(endpoint, data=payload, timeout=120)
                response.raise_for_status()
                rows.extend(response.json())
                last_error = None
                break
            except (requests.RequestException, ValueError) as exc:
                last_error = exc
                time.sleep(2 ** attempt)
        if last_error is not None:
            raise RuntimeError(f"MyGene.info query failed for batch {start}: {last_error}")
        print(f"Annotated/query attempted: {min(start + batch_size, len(gene_ids))}/{len(gene_ids)}", flush=True)
        time.sleep(0.2)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True,
                        help="Donor-aware result table with gene_id and padj columns")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--top", type=int, default=50)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    results = pd.read_csv(args.results, sep=None, engine="python")
    significant = results.loc[results["padj"] < 0.05, "gene_id"].tolist()
    leading = results.nsmallest(args.top, "padj")["gene_id"].tolist()
    targets = list(dict.fromkeys(significant + leading + list(BIOMARKERS)))
    annotation = query_mygene(targets)

    if annotation.empty:
        raise RuntimeError("MyGene.info returned no annotation records.")
    # querymany-style responses can contain duplicate hits. Keep the first hit
    # per queried Ensembl ID and retain an audit table of all returned records.
    annotation.to_csv(args.outdir / "mygene_raw_annotation.tsv", sep="\t", index=False)
    first = annotation.drop_duplicates("query", keep="first").copy()
    first = first.rename(columns={"query": "gene_id", "symbol": "gene_symbol",
                                  "name": "gene_name", "type_of_gene": "gene_type"})
    keep = ["gene_id", "gene_symbol", "gene_name", "gene_type", "notfound"]
    for col in keep:
        if col not in first:
            first[col] = pd.NA
    merged = results.merge(first[keep], on="gene_id", how="left")
    merged.to_csv(args.outdir / "donor_aware_results_annotated.tsv", sep="\t", index=False)

    sig = merged.loc[merged["padj"] < 0.05].sort_values("padj")
    sig.to_csv(args.outdir / "Table_S2_donor_aware_significant_genes.tsv", sep="\t", index=False)
    merged.nsmallest(args.top, "padj").to_csv(
        args.outdir / "Table_S3_top_genes.tsv", sep="\t", index=False)
    missing = merged.loc[merged["gene_symbol"].isna(), "gene_id"].drop_duplicates()
    missing.to_csv(args.outdir / "unmapped_gene_ids.txt", sep="\t", index=False, header=False)
    print(f"Queried {len(targets)} genes; {len(sig)} FDR-significant rows exported; "
          f"{len(missing)} identifiers lack a symbol in the returned annotation.")


if __name__ == "__main__":
    main()
