"""
Builds a small synthetic (non-biological, idealized-geometry) PDB file
purely so the analysis pipeline can be unit-tested without network access
to the real RCSB PDB. Not a real protein structure.
"""

from Bio.PDB import Structure, Model, Chain, Residue, Atom, PDBIO

THREE_LETTER = ["MET", "ALA", "CYS", "GLY", "LEU", "CYS", "VAL", "SER", "LYS", "TRP"]


def build_synthetic_structure():
    structure = Structure.Structure("SYN1")
    model = Model.Model(0)
    structure.add(model)
    chain = Chain.Chain("A")
    model.add(chain)

    bond_len = 1.5
    for i, resname in enumerate(THREE_LETTER):
        res_id = (" ", i + 1, " ")
        residue = Residue.Residue(res_id, resname, "")
        x = i * bond_len * 2
        n_coord = [x, 0.0, 0.0]
        ca_coord = [x + 1.0, 0.3, 0.0]
        c_coord = [x + 2.0, 0.0, 0.0]
        o_coord = [x + 2.0, 1.2, 0.0]

        for name, coord, element in [
            ("N", n_coord, "N"),
            ("CA", ca_coord, "C"),
            ("C", c_coord, "C"),
            ("O", o_coord, "O"),
        ]:
            atom = Atom.Atom(name, coord, 20.0, 1.0, " ", name, i * 10, element)
            residue.add(atom)

        if resname == "CYS":
            # place SG atoms at fixed, nearby coordinates (independent of
            # sequence position) across the two CYS residues, so the
            # disulfide-candidate finder has something real to detect
            sg_coord = [50.0, 50.0, 2.0 if i == 2 else 2.3]
            sg_atom = Atom.Atom("SG", sg_coord, 20.0, 1.0, " ", "SG", i * 10 + 5, "S")
            residue.add(sg_atom)

        chain.add(residue)

    return structure


if __name__ == "__main__":
    structure = build_synthetic_structure()
    io = PDBIO()
    io.set_structure(structure)
    io.save("tests/synthetic.pdb")
    print("Wrote tests/synthetic.pdb")
