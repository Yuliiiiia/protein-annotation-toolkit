"""
analyze.py — Structural and sequence analysis of protein structures.

Given a PDB/mmCIF file, this module extracts:
  - chains, residue counts, hetero-atoms (ligands/waters)
  - amino acid sequence per chain
  - molecular weight, isoelectric point, GRAVY (hydrophobicity)
  - secondary structure (via DSSP, if the `mkdssp` binary is installed)
  - candidate disulfide bonds (close SG-SG cysteine pairs)
  - radius of gyration (a rough measure of compactness)
"""

from pathlib import Path
from statistics import mean

from Bio.PDB import MMCIFParser, PDBParser, PPBuilder, is_aa
from Bio.SeqUtils.ProtParam import ProteinAnalysis

try:
    from Bio.PDB import DSSP
    _HAS_DSSP = True
except ImportError:
    _HAS_DSSP = False


def _get_parser(file_path: Path):
    suffix = file_path.suffix.lower()
    if suffix in (".cif", ".mmcif"):
        return MMCIFParser(QUIET=True)
    return PDBParser(QUIET=True)


def load_structure(file_path: str):
    """Parse a PDB or mmCIF file and return a Bio.PDB Structure object."""
    file_path = Path(file_path)
    parser = _get_parser(file_path)
    structure_id = file_path.stem
    return parser.get_structure(structure_id, str(file_path))


def get_chain_sequences(structure) -> dict:
    """Return {chain_id: amino_acid_sequence} for each polypeptide chain."""
    ppb = PPBuilder()
    sequences = {}
    for model in structure:
        for chain in model:
            seq_parts = [str(pp.get_sequence()) for pp in ppb.build_peptides(chain)]
            seq = "".join(seq_parts)
            if seq:
                sequences[chain.id] = seq
        break  # only first model (NMR files can have many)
    return sequences


def summarize_composition(structure) -> dict:
    """Count residue types, hetero groups, and waters across the structure."""
    aa_count, hetero_count, water_count = 0, 0, 0
    hetero_names = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                het_flag = residue.id[0]
                if het_flag == " " and is_aa(residue, standard=True):
                    aa_count += 1
                elif het_flag == "W":
                    water_count += 1
                elif het_flag.startswith("H_"):
                    hetero_count += 1
                    name = residue.resname
                    hetero_names[name] = hetero_names.get(name, 0) + 1
        break
    return {
        "amino_acid_residues": aa_count,
        "water_molecules": water_count,
        "hetero_groups": hetero_count,
        "hetero_group_types": hetero_names,
    }


def compute_physicochemical_properties(sequence: str) -> dict:
    """Compute molecular weight, isoelectric point, GRAVY, and aromaticity."""
    if not sequence:
        return {}
    analysis = ProteinAnalysis(sequence.replace("X", ""))
    return {
        "length": len(sequence),
        "molecular_weight_da": round(analysis.molecular_weight(), 2),
        "isoelectric_point": round(analysis.isoelectric_point(), 2),
        "gravy_hydrophobicity": round(analysis.gravy(), 3),
        "aromaticity": round(analysis.aromaticity(), 3),
        "instability_index": round(analysis.instability_index(), 2),
        "amino_acid_percent": {
            k: round(v * 100, 1) for k, v in analysis.amino_acids_percent.items() if v > 0
        },
    }


def find_disulfide_candidates(structure, max_distance: float = 2.5) -> list:
    """
    Find pairs of cysteine SG atoms within max_distance (Angstroms) of each
    other — candidate disulfide bonds. A true disulfide bond is ~2.05 A.
    """
    sg_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname == "CYS" and "SG" in residue:
                    sg_atoms.append((chain.id, residue.id[1], residue["SG"]))
        break

    bonds = []
    for i in range(len(sg_atoms)):
        for j in range(i + 1, len(sg_atoms)):
            chain_i, resnum_i, atom_i = sg_atoms[i]
            chain_j, resnum_j, atom_j = sg_atoms[j]
            distance = float(atom_i - atom_j)
            if distance <= max_distance:
                bonds.append({
                    "residue_1": f"{chain_i}:CYS{resnum_i}",
                    "residue_2": f"{chain_j}:CYS{resnum_j}",
                    "distance_angstrom": round(distance, 2),
                })
    return bonds


def compute_radius_of_gyration(structure) -> float:
    """Rough compactness measure: RMS distance of all atoms from centroid."""
    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coords.append(atom.get_coord())
        break
    if not coords:
        return 0.0
    cx = mean(c[0] for c in coords)
    cy = mean(c[1] for c in coords)
    cz = mean(c[2] for c in coords)
    sq_dists = [
        (c[0] - cx) ** 2 + (c[1] - cy) ** 2 + (c[2] - cz) ** 2 for c in coords
    ]
    rg = (sum(sq_dists) / len(sq_dists)) ** 0.5
    return round(float(rg), 2)


def run_dssp(structure, file_path: str) -> dict:
    """
    Run DSSP (if the mkdssp binary is available) to get per-residue secondary
    structure, then summarize into overall percentages.

    Newer mkdssp builds (4.x, e.g. the conda-forge package) parse mmCIF
    files using libcifpp, which expects local chemical-component-dictionary
    data that often isn't configured out of the box. That makes DSSP fail
    on .cif input with a "mmcif_pdbx ... Is a directory" style error even
    though DSSP itself works fine. Classic PDB-format input avoids that
    dictionary lookup entirely, so if the first attempt fails and the input
    was mmCIF, we convert the already-parsed structure to a temporary
    legacy .pdb file and retry once before giving up.

    Returns a dict with "available": False and a "reason" if DSSP is not
    installed or both attempts fail — this is a soft dependency throughout.
    """
    if not _HAS_DSSP:
        return {"available": False, "reason": "Bio.PDB.DSSP not importable"}

    model = structure[0]
    last_error = None

    try:
        dssp = DSSP(model, str(file_path))
        return _summarize_dssp(dssp)
    except Exception as exc:
        last_error = exc

    if Path(file_path).suffix.lower() in (".cif", ".mmcif"):
        tmp_path = None
        try:
            import tempfile
            from Bio.PDB import PDBIO

            with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
                tmp_path = tmp.name
            io = PDBIO()
            io.set_structure(structure)
            io.save(tmp_path)

            dssp = DSSP(model, tmp_path)
            return _summarize_dssp(dssp)
        except Exception as exc:
            last_error = exc
        finally:
            if tmp_path:
                Path(tmp_path).unlink(missing_ok=True)

    return {"available": False, "reason": f"DSSP failed on both mmCIF and PDB input: {last_error}"}


def _summarize_dssp(dssp) -> dict:
    ss_counts = {}
    for key in dssp.keys():
        ss = dssp[key][2]
        ss_counts[ss] = ss_counts.get(ss, 0) + 1

    total = sum(ss_counts.values()) or 1
    # DSSP codes: H/G/I = helix, E/B = sheet/strand, rest = coil/turn/loop
    helix = sum(ss_counts.get(c, 0) for c in "HGI")
    sheet = sum(ss_counts.get(c, 0) for c in "EB")
    coil = total - helix - sheet

    return {
        "available": True,
        "residues_assigned": total,
        "helix_percent": round(100 * helix / total, 1),
        "sheet_percent": round(100 * sheet / total, 1),
        "coil_percent": round(100 * coil / total, 1),
        "raw_counts": ss_counts,
    }


def analyze_structure(file_path: str) -> dict:
    """
    Run the full analysis pipeline on a structure file and return a single
    result dict combining composition, per-chain properties, disulfides,
    radius of gyration, and secondary structure (if DSSP available).
    """
    structure = load_structure(file_path)
    sequences = get_chain_sequences(structure)

    chain_properties = {
        chain_id: compute_physicochemical_properties(seq)
        for chain_id, seq in sequences.items()
    }

    return {
        "structure_id": structure.id,
        "num_chains": len(sequences),
        "chain_sequences": sequences,
        "chain_properties": chain_properties,
        "composition": summarize_composition(structure),
        "disulfide_candidates": find_disulfide_candidates(structure),
        "radius_of_gyration_angstrom": compute_radius_of_gyration(structure),
        "secondary_structure": run_dssp(structure, file_path),
    }
