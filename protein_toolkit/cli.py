"""
cli.py — Command-line entry point.

Usage:
    python -m protein_toolkit.cli 1CRN
    python -m protein_toolkit.cli 1CRN --file data/1crn.cif   # skip download
    python -m protein_toolkit.cli 1CRN --json                 # raw JSON output
    python -m protein_toolkit.cli 1CRN --domains               # + Pfam domain annotation
"""

import argparse
import json
import sys

from .fetch import fetch_structure
from .analyze import analyze_structure
from .domains import annotate_domains
from .report import render_markdown


def main():
    parser = argparse.ArgumentParser(description="Protein structure viewer/annotator")
    parser.add_argument("pdb_id", help="4-character PDB ID, e.g. 1CRN")
    parser.add_argument("--file", help="Use a local structure file instead of downloading")
    parser.add_argument("--out-dir", default="data", help="Directory for downloaded files")
    parser.add_argument("--json", action="store_true", help="Print raw JSON instead of Markdown")
    parser.add_argument("--save", help="Save the report to this file path")
    parser.add_argument(
        "--domains", action="store_true",
        help="Look up Pfam domains via PDBe/InterPro (requires network, adds a few seconds)"
    )
    args = parser.parse_args()

    if args.file:
        structure_path = args.file
    else:
        print(f"Fetching {args.pdb_id} from RCSB PDB...", file=sys.stderr)
        structure_path = fetch_structure(args.pdb_id, out_dir=args.out_dir)

    print("Analyzing structure...", file=sys.stderr)
    results = analyze_structure(structure_path)

    if args.domains:
        print("Looking up Pfam domains...", file=sys.stderr)
        results["domains"] = annotate_domains(args.pdb_id)

    output = json.dumps(results, indent=2) if args.json else render_markdown(results)

    if args.save:
        with open(args.save, "w") as f:
            f.write(output)
        print(f"Report saved to {args.save}", file=sys.stderr)
    else:
        print(output)


if __name__ == "__main__":
    main()
