"""
fetch.py — Download protein structure files from the RCSB Protein Data Bank.
"""

from pathlib import Path
from Bio.PDB import PDBList


def fetch_structure(pdb_id: str, out_dir: str = "data") -> Path:
    """
    Download a structure file (mmCIF format) for the given PDB ID.

    Raises FileNotFoundError with a clear message if the ID doesn't exist
    or the download otherwise fails. Biopython's PDBList doesn't raise on
    a 404 itself -- it just prints a message and returns None -- so
    without this check, callers further down the pipeline (like
    batch.analyze_one's error handling) would get a confusing, unrelated
    TypeError from trying to treat None as a file path instead of a clean,
    catchable error about the real problem.
    """
    pdb_id = pdb_id.strip().upper()
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    pdbl = PDBList()
    file_path = pdbl.retrieve_pdb_file(
        pdb_id, pdir=out_dir, file_format="mmCif"
    )
    if not file_path or not Path(file_path).exists():
        raise FileNotFoundError(
            f"Could not download structure '{pdb_id}' from RCSB PDB "
            "-- check that the PDB ID is valid and currently in the archive"
        )
    return Path(file_path)


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python fetch.py <PDB_ID>")
        sys.exit(1)

    path = fetch_structure(sys.argv[1])
    print(f"Downloaded structure to: {path}")
