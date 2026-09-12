import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from protein_toolkit import analyze
from tests.make_synthetic_pdb import build_synthetic_structure
from Bio.PDB import PDBIO

TEST_FILE = Path(__file__).parent / "synthetic.pdb"


def setup_module(module):
    if not TEST_FILE.exists():
        structure = build_synthetic_structure()
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(TEST_FILE))


def test_load_structure():
    structure = analyze.load_structure(TEST_FILE)
    assert structure is not None
    chains = list(structure[0])
    assert len(chains) == 1


def test_get_chain_sequences():
    structure = analyze.load_structure(TEST_FILE)
    seqs = analyze.get_chain_sequences(structure)
    assert "A" in seqs
    assert seqs["A"] == "MACGLCVSKW"
    assert len(seqs["A"]) == 10


def test_composition():
    structure = analyze.load_structure(TEST_FILE)
    comp = analyze.summarize_composition(structure)
    assert comp["amino_acid_residues"] == 10
    assert comp["water_molecules"] == 0
    assert comp["hetero_groups"] == 0


def test_physicochemical_properties():
    props = analyze.compute_physicochemical_properties("MACGLCVSKW")
    assert props["length"] == 10
    assert props["molecular_weight_da"] > 0
    assert 0 < props["isoelectric_point"] < 14


def test_disulfide_candidates():
    structure = analyze.load_structure(TEST_FILE)
    bonds = analyze.find_disulfide_candidates(structure, max_distance=2.5)
    # the synthetic structure places CYS3 and CYS6 SG atoms ~0.3 A apart
    assert len(bonds) == 1
    assert "CYS3" in bonds[0]["residue_1"] or "CYS3" in bonds[0]["residue_2"]


def test_radius_of_gyration():
    structure = analyze.load_structure(TEST_FILE)
    rg = analyze.compute_radius_of_gyration(structure)
    assert rg > 0


def test_full_pipeline_runs():
    results = analyze.analyze_structure(TEST_FILE)
    assert results["num_chains"] == 1
    assert "A" in results["chain_sequences"]
    assert results["composition"]["amino_acid_residues"] == 10
    assert isinstance(results["disulfide_candidates"], list)
