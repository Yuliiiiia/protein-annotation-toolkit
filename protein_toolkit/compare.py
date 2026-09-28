"""
compare.py — Structural alignment and RMSD between two protein structures.

Comparing two structures directly by residue *number* breaks whenever the
numbering differs (different starting residue, insertions, a construct
with extra tag residues, etc.) — which is the normal case, not the
exception, when comparing two independently-deposited PDB entries. So this
module first does a sequence alignment (Biopython's pairwise aligner,
BLOSUM62) to find which residues actually correspond to each other, then
only superimposes CA atoms at those matched, non-gap positions using
Bio.PDB.Superimposer, which finds the rotation+translation that minimizes
RMSD (the standard, correct way to compare two structures regardless of
how they happen to be oriented/positioned in space).
"""

from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import PPBuilder, Superimposer


def _get_chain_residues(structure, chain_id: str):
    """
    Return (sequence_string, [Residue, ...]) for one chain, in the same
    order PPBuilder emits them -- so seq[i] always corresponds to
    residues[i]. Concatenates all peptide segments in the chain (a chain
    with an internal break still yields one aligned sequence/residue list).
    """
    model = structure[0]
    if chain_id not in model:
        available = ", ".join(c.id for c in model)
        raise ValueError(f"Chain '{chain_id}' not found (available: {available})")

    chain = model[chain_id]
    ppb = PPBuilder()
    seq_parts, residues = [], []
    for pp in ppb.build_peptides(chain):
        seq_parts.append(str(pp.get_sequence()))
        residues.extend(list(pp))

    return "".join(seq_parts), residues


def _sequence_alignment_pairs(seq1: str, seq2: str):
    """
    Global pairwise alignment (BLOSUM62) between two sequences, returning
    a list of (index_in_seq1, index_in_seq2) for every aligned column that
    isn't a gap in either sequence. Uses the best-scoring alignment.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(seq1, seq2)[0]
    blocks1, blocks2 = alignment.aligned  # parallel arrays of (start, end) ungapped blocks

    pairs = []
    for (s1, e1), (s2, e2) in zip(blocks1, blocks2):
        for offset in range(e1 - s1):
            pairs.append((s1 + offset, s2 + offset))
    return pairs


def compare_structures(
    structure1, structure2,
    chain1: str = "A", chain2: str = "A",
    min_aligned_residues: int = 3,
) -> dict:
    """
    Align chain1 of structure1 against chain2 of structure2 by sequence,
    then compute the RMSD of their CA atoms at matched positions after
    optimal superposition.

    structure1 / structure2 are already-parsed Bio.PDB Structure objects
    (as returned by analyze.load_structure), not file paths -- callers
    that only have file paths should load them first.

    Returns a dict with rmsd_angstrom, residues_aligned, and
    sequence_identity_percent, or raises ValueError if there aren't
    enough alignable CA atoms (e.g. completely unrelated sequences, or a
    chain ID that doesn't exist).
    """
    seq1, residues1 = _get_chain_residues(structure1, chain1)
    seq2, residues2 = _get_chain_residues(structure2, chain2)

    if not seq1 or not seq2:
        raise ValueError("One or both chains have no standard-residue sequence to align")

    pairs = _sequence_alignment_pairs(seq1, seq2)

    atoms1, atoms2 = [], []
    identical = 0
    for i1, i2 in pairs:
        r1, r2 = residues1[i1], residues2[i2]
        if "CA" not in r1 or "CA" not in r2:
            continue
        atoms1.append(r1["CA"])
        atoms2.append(r2["CA"])
        if seq1[i1] == seq2[i2]:
            identical += 1

    if len(atoms1) < min_aligned_residues:
        raise ValueError(
            f"Only {len(atoms1)} alignable CA atoms found between "
            f"{chain1} and {chain2} (need at least {min_aligned_residues}) "
            "-- these chains may not be homologous, or chain IDs may be wrong"
        )

    superimposer = Superimposer()
    superimposer.set_atoms(atoms1, atoms2)  # computes optimal rotation+translation internally

    return {
        "rmsd_angstrom": round(float(superimposer.rms), 3),
        "residues_aligned": len(atoms1),
        "sequence_identity_percent": round(100 * identical / len(atoms1), 1),
        "chain1": chain1,
        "chain2": chain2,
        "seq1_length": len(seq1),
        "seq2_length": len(seq2),
    }


def compare_structure_files(
    file1: str, file2: str,
    chain1: str = "A", chain2: str = "A",
) -> dict:
    """Convenience wrapper: load two structure files by path, then compare."""
    from .analyze import load_structure  # local import avoids a circular import at module load

    structure1 = load_structure(file1)
    structure2 = load_structure(file2)
    return compare_structures(structure1, structure2, chain1=chain1, chain2=chain2)


def main():
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description="Compute RMSD between two protein structures (sequence-aligned, CA-based)"
    )
    parser.add_argument("pdb_id_or_file_1")
    parser.add_argument("pdb_id_or_file_2")
    parser.add_argument("--chain1", default="A", help="Chain ID in the first structure (default: A)")
    parser.add_argument("--chain2", default="A", help="Chain ID in the second structure (default: A)")
    parser.add_argument("--out-dir", default="data", help="Directory for downloaded files")
    args = parser.parse_args()

    from pathlib import Path
    from .fetch import fetch_structure

    def resolve(arg):
        # Treat it as a local file if it exists / looks like a path; otherwise fetch by PDB ID.
        if Path(arg).exists():
            return arg
        return fetch_structure(arg, out_dir=args.out_dir)

    path1 = resolve(args.pdb_id_or_file_1)
    path2 = resolve(args.pdb_id_or_file_2)

    try:
        result = compare_structure_files(path1, path2, chain1=args.chain1, chain2=args.chain2)
    except ValueError as exc:
        print(f"Comparison failed: {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"RMSD: {result['rmsd_angstrom']} Å")
    print(f"Aligned residues: {result['residues_aligned']}")
    print(f"Sequence identity: {result['sequence_identity_percent']}%")


if __name__ == "__main__":
    main()
