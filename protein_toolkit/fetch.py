"""
fetch.py — Download protein structure files from the RCSB Protein Data Bank.
"""

from pathlib import Path
from Bio.PDB import PDBList


def fetch_structure(pdb_id: str, out_dir: str = "data") -> Path:
    """
    Download a structure file (mmCIF format) for the given PDB ID.

    Parameters
    ----------
    pdb_id : str
        4-character PDB accession code, e.g. "1CRN".
    out_dir : str
        Directory to save the downloaded file into.

    Returns
    -------
    Path to the downloaded structure file.
    """
    pdb_id = pdb_id.strip().upper()
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    pdbl = PDBList()
    file_path = pdbl.retrieve_pdb_file(
        pdb_id, pdir=out_dir, file_format="mmCif"
    )
    return Path(file_path)


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python fetch.py <PDB_ID>")
        sys.exit(1)

    path = fetch_structure(sys.argv[1])
    print(f"Downloaded structure to: {path}")
