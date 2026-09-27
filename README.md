# Protein Structure Viewer & Annotator

![SARS-CoV-2 spike protein rendered by the toolkit's interactive viewer](https://github.com/user-attachments/assets/67c2a257-9f99-4b7a-8c26-3d92ec99cf0e)

*Interactive 3D view of the SARS-CoV-2 spike protein (PDB: 6VXX), rendered via the included py3Dmol notebook.*

A small toolkit that downloads a protein structure from the RCSB Protein
Data Bank, computes physicochemical and structural properties, and
generates a readable annotation report — plus an interactive 3D viewer
notebook.

## Why

Given a PDB ID, this answers the questions a biologist usually asks first:
how big is the protein, how many chains does it have, is it hydrophobic or
hydrophilic, does it have disulfide bonds, what's its secondary structure
composition (helix/sheet/coil), and are there bound ligands or cofactors?

## Features

- **Fetch** any structure directly from RCSB by PDB ID (mmCIF format)
- **Composition summary** — residue counts, waters, hetero groups/ligands
- **Physicochemical properties per chain** — molecular weight, isoelectric
  point, GRAVY hydrophobicity score, aromaticity, instability index
- **Secondary structure**: helix/sheet/coil percentages via `pydssp` — a pure
  NumPy DSSP re-implementation with no external binary, so no separate
  system install and no compiled-dictionary headaches
- **Domain annotation** — maps each chain to UniProt (via EBI's SIFTS/PDBe
  mapping) and reports Pfam domain hits (via EBI InterPro) in PDB residue
  numbering
- **Disulfide bond detection** — finds cysteine pairs within bonding distance
- **Radius of gyration** — a simple measure of structural compactness
- **Markdown or JSON report output**
- **Interactive 3D viewer** via a Jupyter notebook (py3Dmol)

## Installation

```bash
git clone https://github.com/<you>/protein-toolkit.git
cd protein-toolkit
pip install -r requirements.txt
```

Secondary structure assignment (`pydssp`) is included in `requirements.txt`
— no separate system package needed. Earlier versions of this project used
the external `mkdssp` binary, which turned out to have a confirmed,
unresolved packaging bug in its conda-forge distribution (see
[TROUBLESHOOTING.md](TROUBLESHOOTING.md) for the full story of chasing that
down) — `pydssp` sidesteps it entirely.

> Setting this up on WSL (or ran into `sudo`/`pip` issues in general)? See
> [TROUBLESHOOTING.md](TROUBLESHOOTING.md) — it walks through every real
> problem hit setting this project up from scratch, and how each was fixed.

## Usage

### Command line

```bash
# Fetch and analyze a structure by PDB ID, print a Markdown report
python -m protein_toolkit.cli 1CRN

# Save the report to a file
python -m protein_toolkit.cli 1CRN --save report.md

# Get raw JSON instead
python -m protein_toolkit.cli 1CRN --json

# Analyze a structure file you already have, skip downloading
python -m protein_toolkit.cli 1CRN --file data/1crn.cif

# Include Pfam domain annotation (needs network access, a few extra seconds)
python -m protein_toolkit.cli 1CRN --domains
```

### Example output

```
# Structure Report: 1crn

## Composition
- Chains: 1
- Amino acid residues: 46
- Water molecules: 0
- Hetero groups (ligands/cofactors): 0

## Radius of Gyration: 10.32 Å

## Per-Chain Properties
### Chain A
- Length: 46 residues
- Molecular weight: 4736.5 Da
- Isoelectric point (pI): 6.4
- GRAVY (hydrophobicity): 0.049
  (negative = hydrophilic, positive = hydrophobic)
- Aromaticity: 0.087
- Instability index: 30.5 (likely stable)

## Secondary Structure
- Helix: 47.8%
- Sheet: 15.2%
- Coil/loop: 37.0%

## Disulfide Bond Candidates
- A:CYS3 — A:CYS40 (2.03 Å)
- A:CYS4 — A:CYS32 (2.05 Å)
- A:CYS16 — A:CYS26 (2.04 Å)

## Pfam Domains
### Chain A
- PF00069 — Protein kinase domain (residues 15–270, via UniProt P00519)
```

### Interactive 3D viewer

Open `examples/view_structure.ipynb` in Jupyter, set `PDB_ID` to any
4-character accession code, and run all cells to get a rotatable cartoon
view plus the full annotation report inline.

### As a library

```python
from protein_toolkit.fetch import fetch_structure
from protein_toolkit.analyze import analyze_structure
from protein_toolkit.domains import annotate_domains
from protein_toolkit.report import render_markdown

path = fetch_structure("1UBQ")
results = analyze_structure(path)
results["domains"] = annotate_domains("1UBQ")
print(render_markdown(results))
```

## Project structure

```
protein_toolkit/
├── fetch.py      # download structures from RCSB PDB
├── analyze.py    # sequence & structural analysis
├── domains.py    # PDB -> UniProt -> Pfam domain annotation
├── report.py     # markdown report rendering
└── cli.py        # command-line interface
examples/
└── view_structure.ipynb   # interactive 3D viewer
tests/
├── make_synthetic_pdb.py     # builds a small offline test fixture
├── 1a8o_fixture.cif          # small real structure (bundled, no network needed)
├── test_analyze.py           # unit tests (no network required)
├── test_domains.py           # domain annotation tests (mocked API responses)
└── test_secondary_structure.py  # pydssp tests, incl. a real-data regression check
TROUBLESHOOTING.md             # real setup issues hit and how they were fixed
```

## Running tests

```bash
pytest tests/ -v
```

Tests run entirely offline: against a small synthetic structure generated
by `tests/make_synthetic_pdb.py`, a small real structure bundled as
`tests/1a8o_fixture.cif` (so the secondary-structure regression test checks
against genuine coordinates, not synthetic ones), and mocked API responses
(`unittest.mock`) for the domain-annotation tests — no network access or
real PDB download required for CI.

## Possible extensions

- Structural alignment between two structures (RMSD via `Bio.PDB.Superimposer`)
- Batch mode: analyze a list of PDB IDs and output a comparison table
- Web frontend (Flask/FastAPI + py3Dmol.js) instead of notebook-only viewing

## License

MIT
