import copy
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from protein_toolkit import compare
from protein_toolkit.analyze import load_structure

FIXTURE = Path(__file__).parent / "1a8o_fixture.cif"


def _random_rotation_matrix(seed=42):
    """A real, non-trivial 3D rotation matrix (via QR decomposition of a
    random matrix), used to prove RMSD-after-superposition actually
    corrects for arbitrary orientation rather than assuming identity."""
    rng = np.random.default_rng(seed)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    return q.astype(np.float32)


def _apply_rigid_transform(structure, rotation, translation):
    """Return a deep copy of structure with every atom coordinate
    rotated and translated by a fixed rigid-body transform."""
    moved = copy.deepcopy(structure)
    for atom in moved.get_atoms():
        atom.coord = rotation @ atom.coord + translation
    return moved


def test_compare_structure_against_itself_gives_near_zero_rmsd():
    structure = load_structure(FIXTURE)
    result = compare.compare_structures(structure, structure, chain1="A", chain2="A")

    assert result["rmsd_angstrom"] < 0.01
    assert result["sequence_identity_percent"] == 100.0
    assert result["residues_aligned"] > 0


def test_rmsd_is_invariant_to_rigid_body_transform():
    """
    The real point of using Superimposer instead of comparing raw
    coordinates: two copies of the same structure that have been
    arbitrarily rotated and translated relative to each other should
    still report ~0 RMSD, because Superimposer finds the optimal
    alignment first. This is the actual behavior that makes RMSD
    comparison meaningful for two independently-deposited PDB entries,
    which are never in the same coordinate frame to begin with.
    """
    structure = load_structure(FIXTURE)
    rotation = _random_rotation_matrix()
    translation = np.array([37.0, -12.0, 8.5], dtype=np.float32)
    transformed = _apply_rigid_transform(structure, rotation, translation)

    result = compare.compare_structures(structure, transformed, chain1="A", chain2="A")

    assert result["rmsd_angstrom"] < 0.01
    assert result["sequence_identity_percent"] == 100.0


def test_rmsd_reflects_real_coordinate_perturbation():
    """
    A structure compared against a copy with real per-atom noise added
    should report a small but genuinely non-zero RMSD -- confirming the
    metric is sensitive to actual structural difference, not just
    silently returning 0 regardless of input.
    """
    structure = load_structure(FIXTURE)
    noisy = copy.deepcopy(structure)
    rng = np.random.default_rng(7)
    for atom in noisy.get_atoms():
        atom.coord = atom.coord + rng.normal(scale=0.5, size=3).astype(np.float32)

    result = compare.compare_structures(structure, noisy, chain1="A", chain2="A")

    assert 0.05 < result["rmsd_angstrom"] < 2.0


def test_compare_raises_on_missing_chain():
    structure = load_structure(FIXTURE)
    try:
        compare.compare_structures(structure, structure, chain1="Z", chain2="A")
        assert False, "expected ValueError for missing chain"
    except ValueError as exc:
        assert "Z" in str(exc)


def test_compare_raises_on_too_few_aligned_residues():
    structure = load_structure(FIXTURE)
    seq1, residues1 = compare._get_chain_residues(structure, "A")
    # Build a structure/sequence pair with essentially nothing in common
    # by comparing against a trivially short synthetic sequence/residue set.
    try:
        compare.compare_structures(structure, structure, chain1="A", chain2="A",
                                    min_aligned_residues=10_000)
        assert False, "expected ValueError when requiring more aligned residues than exist"
    except ValueError as exc:
        assert "alignable CA atoms" in str(exc)


def test_sequence_alignment_pairs_handles_identical_sequences():
    pairs = compare._sequence_alignment_pairs("ACDEFG", "ACDEFG")
    assert pairs == [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5)]


def test_sequence_alignment_pairs_handles_insertion():
    # seq2 has an extra residue inserted in the middle relative to seq1
    pairs = compare._sequence_alignment_pairs("ACDFG", "ACDEFG")
    aligned_seq1_chars = [p[0] for p in pairs]
    # every position in the shorter sequence should still find a match
    assert len(pairs) >= 4
    assert 0 in aligned_seq1_chars and 4 in aligned_seq1_chars
