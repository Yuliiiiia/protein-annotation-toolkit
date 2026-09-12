import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from protein_toolkit import analyze

TEST_FILE = Path(__file__).parent / "synthetic.pdb"


def test_run_dssp_reports_unavailable_when_binary_missing():
    """
    With no mkdssp binary present at all, both the direct mmCIF/PDB attempt
    and the PDB-fallback attempt should fail, and the function should report
    that clearly rather than raising.
    """
    structure = analyze.load_structure(TEST_FILE)
    result = analyze.run_dssp(structure, TEST_FILE)
    assert result["available"] is False
    assert "reason" in result


def test_run_dssp_falls_back_to_pdb_conversion_on_cif_failure():
    """
    Simulate the real-world failure mode: DSSP raising on a .cif file
    (as newer mkdssp builds do without chemical-component-dictionary data
    configured) — the second attempt, on a converted temporary .pdb file,
    should be what actually gets called and should succeed.
    """
    structure = analyze.load_structure(TEST_FILE)
    fake_cif_path = "fake_structure.cif"

    calls = []

    def fake_dssp_constructor(model, path):
        calls.append(path)
        if str(path).endswith(".cif"):
            raise Exception("mmcif_pdbx dictionary error: Is a directory")
        mock_dssp = MagicMock()
        mock_dssp.keys.return_value = [("A", 1), ("A", 2), ("A", 3)]
        mock_dssp.__getitem__.side_effect = lambda k: (None, None, "H")
        return mock_dssp

    with patch("protein_toolkit.analyze.DSSP", side_effect=fake_dssp_constructor):
        result = analyze.run_dssp(structure, fake_cif_path)

    assert result["available"] is True
    assert result["helix_percent"] == 100.0
    # first call used the original .cif path, second used a converted .pdb file
    assert calls[0] == fake_cif_path
    assert calls[1].endswith(".pdb")
