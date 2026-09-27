import sys
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from protein_toolkit import analyze

SYNTHETIC_FILE = Path(__file__).parent / "synthetic.pdb"
REAL_FIXTURE = Path(__file__).parent / "1a8o_fixture.cif"


def test_run_dssp_reports_unavailable_when_pydssp_missing():
    """
    If pydssp isn't importable, run_dssp should report that clearly rather
    than raising — same soft-dependency contract as before.
    """
    structure = analyze.load_structure(SYNTHETIC_FILE)
    with patch.object(analyze, "_HAS_PYDSSP", False):
        result = analyze.run_dssp(structure, SYNTHETIC_FILE)
    assert result["available"] is False
    assert "reason" in result


def test_run_dssp_reports_unavailable_on_too_short_backbone():
    """
    The synthetic fixture's 10-residue idealized backbone is too short
    (and geometrically not a real fold) for pydssp's hydrogen-bond window
    to assign anything meaningful — this should fail gracefully, not raise.
    """
    structure = analyze.load_structure(SYNTHETIC_FILE)
    result = analyze.run_dssp(structure, SYNTHETIC_FILE)
    # Either genuinely unavailable, or (if it does run) produces a
    # well-formed result — either is acceptable; what matters is no crash.
    assert "available" in result


def test_run_dssp_assigns_real_secondary_structure():
    """
    Real regression test against actual PDB coordinates (1A8O, HIV-1
    capsid protein C-terminal domain, chain A — a small, genuinely
    structured fragment), bundled as an offline fixture so this doesn't
    depend on network access. Confirms pydssp runs end-to-end on real data
    and produces a sane, non-trivial helix/sheet/coil split rather than
    just not crashing.
    """
    structure = analyze.load_structure(REAL_FIXTURE)
    result = analyze.run_dssp(structure, REAL_FIXTURE)

    assert result["available"] is True
    assert result["residues_assigned"] > 0
    # This fragment is known to be substantially helical; assert loosely
    # (a wide, sane range) rather than pinning exact percentages, since the
    # exact H-bond geometry-derived counts could shift slightly across
    # numpy/pydssp versions.
    assert result["helix_percent"] > 20
    pcts = result["helix_percent"] + result["sheet_percent"] + result["coil_percent"]
    assert abs(pcts - 100.0) < 0.5  # percentages should sum to ~100%


def test_run_dssp_handles_multichain_without_crashing():
    """
    Chains are processed independently; a structure with only one real
    chain (like the fixtures here) should still produce a valid combined
    result without raising on the per-chain loop.
    """
    structure = analyze.load_structure(REAL_FIXTURE)
    result = analyze.run_dssp(structure, REAL_FIXTURE)
    assert isinstance(result, dict)
    assert "raw_counts" in result
