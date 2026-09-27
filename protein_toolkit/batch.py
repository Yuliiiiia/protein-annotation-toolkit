"""
batch.py — Analyze several PDB structures in one run and produce a
side-by-side comparison table (Markdown and/or CSV).

Usage:
    python -m protein_toolkit.batch 1CRN 1UBQ 6VXX
    python -m protein_toolkit.batch 1CRN 1UBQ 6VXX --domains
    python -m protein_toolkit.batch 1CRN 1UBQ 6VXX --csv comparison.csv

A single failing PDB ID (bad accession, network hiccup, malformed
structure) does not abort the batch — its row just records the error and
every other structure is still analyzed and reported.
"""

import csv as csv_module
from pathlib import Path

from .fetch import fetch_structure
from .analyze import analyze_structure
from .domains import annotate_domains

FIELDNAMES = [
    "pdb_id", "chains", "residues", "hetero_groups", "primary_chain",
    "mw_da", "pI", "gravy", "instability",
    "helix_pct", "sheet_pct", "disulfides", "radius_of_gyration",
    "pfam_domains", "error",
]


def analyze_one(pdb_id: str, out_dir: str = "data", include_domains: bool = False) -> dict:
    """
    Fetch + analyze a single PDB ID, returning one flat summary row.
    Never raises: any failure is captured in the row's "error" field so a
    batch of many IDs isn't derailed by one bad entry.
    """
    row = {name: None for name in FIELDNAMES}
    row["pdb_id"] = pdb_id.strip().upper()

    try:
        path = fetch_structure(pdb_id, out_dir=out_dir)
        results = analyze_structure(path)
    except Exception as exc:
        row["error"] = str(exc)
        return row

    comp = results["composition"]

    # Summarize at the structure level using its longest chain, so a
    # multi-chain assembly still collapses to one readable comparison row.
    primary_chain_id, primary_props = None, {}
    for chain_id, props in results["chain_properties"].items():
        if props and props.get("length", 0) > primary_props.get("length", 0):
            primary_chain_id, primary_props = chain_id, props

    ss = results["secondary_structure"]

    row.update({
        "chains": results["num_chains"],
        "residues": comp["amino_acid_residues"],
        "hetero_groups": comp["hetero_groups"],
        "primary_chain": primary_chain_id,
        "mw_da": primary_props.get("molecular_weight_da"),
        "pI": primary_props.get("isoelectric_point"),
        "gravy": primary_props.get("gravy_hydrophobicity"),
        "instability": primary_props.get("instability_index"),
        "helix_pct": ss.get("helix_percent") if ss.get("available") else None,
        "sheet_pct": ss.get("sheet_percent") if ss.get("available") else None,
        "disulfides": len(results["disulfide_candidates"]),
        "radius_of_gyration": results["radius_of_gyration_angstrom"],
    })

    if include_domains:
        dom = annotate_domains(pdb_id)
        if dom.get("available"):
            hits = sorted({h["pfam_id"] for hits in dom["chains"].values() for h in hits})
            row["pfam_domains"] = ", ".join(hits) if hits else "-"
        # if lookup failed, leave as None -> rendered as "-" downstream

    return row


def analyze_batch(pdb_ids: list, out_dir: str = "data", include_domains: bool = False) -> list:
    """Run analyze_one() across a list of PDB IDs, in order."""
    return [
        analyze_one(pdb_id, out_dir=out_dir, include_domains=include_domains)
        for pdb_id in pdb_ids
    ]


def render_markdown_table(rows: list, include_domains: bool = False) -> str:
    headers = [
        "PDB ID", "Chains", "Residues", "MW (Da)", "pI", "GRAVY",
        "Instability", "Helix %", "Sheet %", "Disulfides", "Rg (Å)",
    ]
    if include_domains:
        headers.append("Pfam domains")

    lines = [
        "| " + " | ".join(headers) + " |",
        "|" + "|".join(["---"] * len(headers)) + "|",
    ]

    for row in rows:
        if row.get("error"):
            cells = [row["pdb_id"], f"**ERROR:** {row['error']}"] + [""] * (len(headers) - 2)
        else:
            def fmt(value):
                return "-" if value is None else str(value)

            cells = [
                row["pdb_id"], fmt(row["chains"]), fmt(row["residues"]),
                fmt(row["mw_da"]), fmt(row["pI"]), fmt(row["gravy"]),
                fmt(row["instability"]), fmt(row["helix_pct"]), fmt(row["sheet_pct"]),
                fmt(row["disulfides"]), fmt(row["radius_of_gyration"]),
            ]
            if include_domains:
                cells.append(fmt(row["pfam_domains"]))
        lines.append("| " + " | ".join(cells) + " |")

    return "\n".join(lines)


def write_csv(rows: list, path: str, include_domains: bool = False) -> None:
    fieldnames = [f for f in FIELDNAMES if include_domains or f != "pfam_domains"]
    with open(path, "w", newline="") as f:
        writer = csv_module.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Analyze multiple PDB structures and produce a comparison table"
    )
    parser.add_argument("pdb_ids", nargs="+", help="Two or more PDB IDs, e.g. 1CRN 1UBQ 6VXX")
    parser.add_argument("--out-dir", default="data", help="Directory for downloaded files")
    parser.add_argument(
        "--domains", action="store_true",
        help="Also look up Pfam domains per structure (requires network, adds time per ID)"
    )
    parser.add_argument("--csv", help="Also write the comparison table to this CSV file")
    parser.add_argument("--save", help="Save the Markdown table to this file instead of printing")
    args = parser.parse_args()

    import sys
    print(f"Analyzing {len(args.pdb_ids)} structures...", file=sys.stderr)
    rows = analyze_batch(args.pdb_ids, out_dir=args.out_dir, include_domains=args.domains)

    failed = [r["pdb_id"] for r in rows if r.get("error")]
    if failed:
        print(f"Warning: {len(failed)} failed: {', '.join(failed)}", file=sys.stderr)

    table = render_markdown_table(rows, include_domains=args.domains)

    if args.csv:
        write_csv(rows, args.csv, include_domains=args.domains)
        print(f"CSV written to {args.csv}", file=sys.stderr)

    if args.save:
        with open(args.save, "w") as f:
            f.write(table)
        print(f"Markdown table saved to {args.save}", file=sys.stderr)
    else:
        print(table)


if __name__ == "__main__":
    main()
