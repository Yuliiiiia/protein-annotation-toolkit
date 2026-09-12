import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from protein_toolkit import domains


def _mock_response(json_data, status_ok=True):
    mock = MagicMock()
    mock.json.return_value = json_data
    if not status_ok:
        mock.raise_for_status.side_effect = Exception("HTTP error")
    else:
        mock.raise_for_status.return_value = None
    return mock


PDBE_MAPPING_RESPONSE = {
    "1abc": {
        "UniProt": {
            "P00519": {
                "identifier": "ABL1_HUMAN",
                "name": "Tyrosine-protein kinase ABL1",
                "mappings": [
                    {
                        "chain_id": "A",
                        "struct_asym_id": "A",
                        "unp_start": 10,
                        "unp_end": 300,
                        "start": {"author_residue_number": 1},
                        "end": {"author_residue_number": 291},
                    }
                ],
            }
        }
    }
}

INTERPRO_PFAM_RESPONSE = {
    "results": [
        {
            "metadata": {"accession": "PF00069", "name": "Protein kinase domain"},
            "proteins": [
                {
                    "accession": "P00519",
                    "entry_protein_locations": [
                        {"fragments": [{"start": 15, "end": 270}]}
                    ],
                }
            ],
        },
        {
            # a domain entirely outside the structure's covered UniProt range
            # (unp 10-300) -- should be dropped by the mapping step
            "metadata": {"accession": "PF99999", "name": "Unrelated domain"},
            "proteins": [
                {
                    "accession": "P00519",
                    "entry_protein_locations": [
                        {"fragments": [{"start": 500, "end": 520}]}
                    ],
                }
            ],
        },
    ]
}


def test_get_uniprot_mappings_success():
    with patch("protein_toolkit.domains.requests.get") as mock_get:
        mock_get.return_value = _mock_response(PDBE_MAPPING_RESPONSE)
        result = domains.get_uniprot_mappings("1abc")

    assert "A" in result
    segment = result["A"][0]
    assert segment["uniprot_accession"] == "P00519"
    assert segment["unp_start"] == 10
    assert segment["pdb_start"] == 1


def test_get_uniprot_mappings_network_failure():
    with patch("protein_toolkit.domains.requests.get", side_effect=Exception("timeout")):
        result = domains.get_uniprot_mappings("1abc")

    assert "_error" in result


def test_get_pfam_domains_parses_fragments():
    with patch("protein_toolkit.domains.requests.get") as mock_get:
        mock_get.return_value = _mock_response(INTERPRO_PFAM_RESPONSE)
        hits = domains.get_pfam_domains("P00519")

    assert len(hits) == 2
    assert hits[0]["pfam_id"] == "PF00069"
    assert hits[0]["unp_start"] == 15
    assert hits[0]["unp_end"] == 270


def test_map_unp_range_to_pdb_offset():
    segment = {"unp_start": 10, "unp_end": 300, "pdb_start": 1, "pdb_end": 291}
    mapped = domains._map_unp_range_to_pdb(15, 270, segment)
    assert mapped == (6, 261)  # offset = 1 - 10 = -9


def test_map_unp_range_to_pdb_no_overlap():
    segment = {"unp_start": 10, "unp_end": 300, "pdb_start": 1, "pdb_end": 291}
    mapped = domains._map_unp_range_to_pdb(500, 520, segment)
    assert mapped is None


def test_annotate_domains_full_pipeline():
    with patch("protein_toolkit.domains.requests.get") as mock_get:
        def side_effect(url, params=None, timeout=None):
            if "mappings/uniprot" in url:
                return _mock_response(PDBE_MAPPING_RESPONSE)
            return _mock_response(INTERPRO_PFAM_RESPONSE)

        mock_get.side_effect = side_effect
        result = domains.annotate_domains("1abc")

    assert result["available"] is True
    assert "A" in result["chains"]
    hits = result["chains"]["A"]
    # only the PF00069 hit overlaps the mapped UniProt range (10-300);
    # PF99999 (500-520) falls outside it and should be dropped
    assert len(hits) == 1
    assert hits[0]["pfam_id"] == "PF00069"
    assert hits[0]["pdb_start"] == 6
    assert hits[0]["pdb_end"] == 261


def test_annotate_domains_handles_missing_entry():
    with patch("protein_toolkit.domains.requests.get") as mock_get:
        mock_get.return_value = _mock_response({})
        result = domains.annotate_domains("9zzz")

    assert result["available"] is False
    assert "reason" in result
