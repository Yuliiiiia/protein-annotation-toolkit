import sys
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from protein_toolkit import webapp

FAKE_RESULTS = {
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
    ],
    "radius_of_gyration_angstrom": 9.67,
    "secondary_structure": {
        "available": True, "residues_assigned": 46,
        "helix_percent": 43.5, "sheet_percent": 8.7, "coil_percent": 47.8,
        "raw_counts": {"H": 20, "E": 4, "-": 22},
    },
}


def _client():
    webapp.app.config["TESTING"] = True
    return webapp.app.test_client()


def test_index_page_loads():
    client = _client()
    resp = client.get("/")
    assert resp.status_code == 200
    assert b"PDB ID" in resp.data


def test_structure_lookup_form_redirects_to_structure_view():
    client = _client()
    resp = client.post("/structure", data={"pdb_id": "1crn"})
    assert resp.status_code == 302
    assert "/structure/1crn" in resp.headers["Location"]


def test_structure_lookup_with_empty_id_redirects_home_with_flash():
    client = _client()
    resp = client.post("/structure", data={"pdb_id": ""}, follow_redirects=True)
    assert resp.status_code == 200
    assert b"Enter a PDB ID" in resp.data


def test_structure_view_renders_report_and_viewer():
    client = _client()
    with patch.object(webapp, "fetch_structure", return_value="tests/1a8o_fixture.cif"), \
         patch.object(webapp, "analyze_structure", return_value=FAKE_RESULTS):
        resp = client.get("/structure/1crn")

    assert resp.status_code == 200
    body = resp.data.decode()
    # report content made it into the page
    assert "4736.43" in body       # molecular weight
    assert "43.5" in body          # helix percent
    assert "CYS3" in body          # disulfide candidate
    # the 3D viewer script and structure payload are present
    assert "3Dmol" in body
    assert "$3Dmol.createViewer" in body


def test_structure_view_includes_domains_when_requested():
    fake_domains = {"available": True, "chains": {"A": [
        {"pfam_id": "PF00321", "name": "Plant thionin", "pdb_start": 2, "pdb_end": 46,
         "uniprot_accession": "P01542"},
    ]}}
    client = _client()
    with patch.object(webapp, "fetch_structure", return_value="tests/1a8o_fixture.cif"), \
         patch.object(webapp, "analyze_structure", return_value=FAKE_RESULTS), \
         patch.object(webapp, "annotate_domains", return_value=fake_domains):
        resp = client.get("/structure/1crn?domains=1")

    body = resp.data.decode()
    assert "PF00321" in body
    assert "Plant thionin" in body


def test_structure_view_omits_domains_section_by_default():
    client = _client()
    with patch.object(webapp, "fetch_structure", return_value="tests/1a8o_fixture.cif"), \
         patch.object(webapp, "analyze_structure", return_value=FAKE_RESULTS):
        resp = client.get("/structure/1crn")  # no ?domains=1

    body = resp.data.decode()
    assert "Pfam domains" not in body


def test_structure_view_handles_fetch_failure_gracefully():
    client = _client()
    with patch.object(webapp, "fetch_structure", side_effect=FileNotFoundError("bad ID")):
        resp = client.get("/structure/ZZZZ", follow_redirects=True)

    assert resp.status_code == 200
    assert b"bad ID" in resp.data
