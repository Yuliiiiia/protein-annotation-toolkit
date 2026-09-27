import sys
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from protein_toolkit import batch


FAKE_RESULTS_OK = {
    "structure_id": "1crn",
    "num_chains": 1,
    "chain_sequences": {"A": "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"},
    "chain_properties": {
        "A": {
            "length": 46, "molecular_weight_da": 4736.43,
            "isoelectric_point": 5.73, "gravy_hydrophobicity": 0.37,
            "aromaticity": 0.065, "instability_index": 31.89,
            "amino_acid_percent": {},
        }
    },
    "composition": {
        "amino_acid_residues": 46, "water_molecules": 0,
        "hetero_groups": 0, "hetero_group_types": {},
    },
    "disulfide_candidates": [
        {"residue_1": "A:CYS3", "residue_2": "A:CYS40", "distance_angstrom": 2.0},
        {"residue_1": "A:CYS4", "residue_2": "A:CYS32", "distance_angstrom": 2.03},
    ],
    "radius_of_gyration_angstrom": 9.67,
    "secondary_structure": {
        "available": True, "residues_assigned": 46,
        "helix_percent": 43.5, "sheet_percent": 8.7, "coil_percent": 47.8,
        "raw_counts": {"H": 20, "E": 4, "-": 22},
    },
}


def test_analyze_one_success():
    with patch.object(batch, "fetch_structure", return_value="data/1crn.cif"), \
         patch.object(batch, "analyze_structure", return_value=FAKE_RESULTS_OK):
        row = batch.analyze_one("1crn")

    assert row["pdb_id"] == "1CRN"
    assert row["error"] is None
    assert row["chains"] == 1
    assert row["residues"] == 46
    assert row["primary_chain"] == "A"
    assert row["mw_da"] == 4736.43
    assert row["helix_pct"] == 43.5
    assert row["disulfides"] == 2


def test_analyze_one_picks_longest_chain_as_primary():
    multi_chain_results = dict(FAKE_RESULTS_OK)
    multi_chain_results["chain_properties"] = {
        "A": {"length": 10, "molecular_weight_da": 1000.0, "isoelectric_point": 6.0,
              "gravy_hydrophobicity": 0.1, "aromaticity": 0.0, "instability_index": 20.0,
              "amino_acid_percent": {}},
        "B": {"length": 200, "molecular_weight_da": 22000.0, "isoelectric_point": 7.0,
              "gravy_hydrophobicity": -0.2, "aromaticity": 0.1, "instability_index": 35.0,
              "amino_acid_percent": {}},
    }
    with patch.object(batch, "fetch_structure", return_value="data/x.cif"), \
         patch.object(batch, "analyze_structure", return_value=multi_chain_results):
        row = batch.analyze_one("XXXX")

    assert row["primary_chain"] == "B"
    assert row["mw_da"] == 22000.0


def test_analyze_one_handles_fetch_failure_gracefully():
    """
    A bad PDB ID or network failure should produce a row with an error
    message, not raise -- this is what lets the rest of a batch continue.
    """
    with patch.object(batch, "fetch_structure", side_effect=Exception("404: unknown PDB ID")):
        row = batch.analyze_one("ZZZZ")

    assert row["pdb_id"] == "ZZZZ"
    assert row["error"] == "404: unknown PDB ID"
    assert row["chains"] is None


def test_analyze_one_handles_analyze_failure_gracefully():
    with patch.object(batch, "fetch_structure", return_value="data/x.cif"), \
         patch.object(batch, "analyze_structure", side_effect=Exception("corrupt structure file")):
        row = batch.analyze_one("XXXX")

    assert row["error"] == "corrupt structure file"


def test_analyze_one_with_domains():
    fake_domains = {"available": True, "chains": {"A": [
        {"pfam_id": "PF00321", "name": "Plant thionin", "pdb_start": 2, "pdb_end": 46,
         "uniprot_accession": "P01542"},
    ]}}
    with patch.object(batch, "fetch_structure", return_value="data/1crn.cif"), \
         patch.object(batch, "analyze_structure", return_value=FAKE_RESULTS_OK), \
         patch.object(batch, "annotate_domains", return_value=fake_domains):
        row = batch.analyze_one("1crn", include_domains=True)

    assert row["pfam_domains"] == "PF00321"


def test_analyze_batch_runs_all_ids_even_with_one_failure():
    def fake_fetch(pdb_id, out_dir="data"):
        if pdb_id == "BAD1":
            raise Exception("not found")
        return f"data/{pdb_id}.cif"

    with patch.object(batch, "fetch_structure", side_effect=fake_fetch), \
         patch.object(batch, "analyze_structure", return_value=FAKE_RESULTS_OK):
        rows = batch.analyze_batch(["1CRN", "BAD1", "1UBQ"])

    assert len(rows) == 3
    assert rows[0]["error"] is None
    assert rows[1]["error"] == "not found"
    assert rows[2]["error"] is None


def test_render_markdown_table_basic():
    with patch.object(batch, "fetch_structure", return_value="data/1crn.cif"), \
         patch.object(batch, "analyze_structure", return_value=FAKE_RESULTS_OK):
        row = batch.analyze_one("1crn")

    table = batch.render_markdown_table([row])
    assert "| PDB ID |" in table
    assert "1CRN" in table
    assert "43.5" in table  # helix percent made it into the table


def test_render_markdown_table_shows_errors():
    error_row = {name: None for name in batch.FIELDNAMES}
    error_row["pdb_id"] = "ZZZZ"
    error_row["error"] = "404: unknown PDB ID"

    table = batch.render_markdown_table([error_row])
    assert "ZZZZ" in table
    assert "ERROR" in table
    assert "404" in table


def test_write_csv(tmp_path):
    with patch.object(batch, "fetch_structure", return_value="data/1crn.cif"), \
         patch.object(batch, "analyze_structure", return_value=FAKE_RESULTS_OK):
        row = batch.analyze_one("1crn")

    csv_path = tmp_path / "out.csv"
    batch.write_csv([row], str(csv_path))

    content = csv_path.read_text()
    assert "1CRN" in content
    assert "pdb_id" in content.splitlines()[0]  # header row present
