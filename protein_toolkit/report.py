"""
report.py — Render analysis results as a readable Markdown report.
"""


def render_markdown(results: dict) -> str:
    lines = [f"# Structure Report: {results['structure_id']}", ""]

    comp = results["composition"]
    lines += [
        "## Composition",
        f"- Chains: {results['num_chains']}",
        f"- Amino acid residues: {comp['amino_acid_residues']}",
        f"- Water molecules: {comp['water_molecules']}",
        f"- Hetero groups (ligands/cofactors): {comp['hetero_groups']}",
    ]
    if comp["hetero_group_types"]:
        het_list = ", ".join(f"{k} ({v})" for k, v in comp["hetero_group_types"].items())
        lines.append(f"  - Types: {het_list}")
    lines.append("")

    lines.append(f"## Radius of Gyration: {results['radius_of_gyration_angstrom']} Å")
    lines.append("")

    lines.append("## Per-Chain Properties")
    for chain_id, props in results["chain_properties"].items():
        if not props:
            continue
        lines += [
            f"### Chain {chain_id}",
            f"- Length: {props['length']} residues",
            f"- Molecular weight: {props['molecular_weight_da']} Da",
            f"- Isoelectric point (pI): {props['isoelectric_point']}",
            f"- GRAVY (hydrophobicity): {props['gravy_hydrophobicity']}",
            f"  (negative = hydrophilic, positive = hydrophobic)",
            f"- Aromaticity: {props['aromaticity']}",
            f"- Instability index: {props['instability_index']}"
            f" ({'likely unstable' if props['instability_index'] > 40 else 'likely stable'})",
            "",
        ]

    ss = results["secondary_structure"]
    lines.append("## Secondary Structure")
    if ss.get("available"):
        lines += [
            f"- Helix: {ss['helix_percent']}%",
            f"- Sheet: {ss['sheet_percent']}%",
            f"- Coil/loop: {ss['coil_percent']}%",
            f"- (based on {ss['residues_assigned']} DSSP-assigned residues)",
        ]
    else:
        lines.append(
            f"- Not available ({ss.get('reason', 'DSSP not installed')}). "
            "Install `mkdssp` (e.g. `apt install dssp` or `conda install -c salilab dssp`) "
            "to enable helix/sheet/coil assignment."
        )
    lines.append("")

    lines.append("## Disulfide Bond Candidates")
    bonds = results["disulfide_candidates"]
    if bonds:
        for b in bonds:
            lines.append(
                f"- {b['residue_1']} — {b['residue_2']} ({b['distance_angstrom']} Å)"
            )
    else:
        lines.append("- None found within 2.5 Å")
    lines.append("")

    domains = results.get("domains")
    if domains is not None:
        lines.append("## Pfam Domains")
        if domains.get("available") and domains.get("chains"):
            for chain_id, hits in domains["chains"].items():
                lines.append(f"### Chain {chain_id}")
                for d in hits:
                    lines.append(
                        f"- {d['pfam_id']} — {d['name']} "
                        f"(residues {d['pdb_start']}–{d['pdb_end']}, "
                        f"via UniProt {d['uniprot_accession']})"
                    )
        elif domains.get("available"):
            lines.append("- No Pfam domains found for this entry")
        else:
            lines.append(f"- Not available ({domains.get('reason', 'lookup failed')})")

    return "\n".join(lines)
