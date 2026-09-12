"""
domains.py — Domain annotation via PDB -> UniProt -> Pfam.

Two public APIs are chained together:

1. EBI PDBe SIFTS mapping API
   https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}
   Gives, for each chain, which UniProt accession(s) it corresponds to and
   the residue-range correspondence between PDB numbering and UniProt
   numbering (a structure often covers only part of the full UniProt
   sequence, and the numbering rarely lines up 1:1 from residue 1).

2. EBI InterPro API (Pfam member database)
   https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{accession}
   Gives Pfam domain hits (accession, name, start/end) in UniProt sequence
   coordinates for that accession.

We combine these two to report Pfam domains in PDB chain/residue-number
coordinates, which is what's actually useful when looking at a structure.

Both calls are network-dependent and best-effort: if either API is
unreachable, unavailable for this entry, or returns an unexpected shape,
functions return an empty result with an "available"/"reason" field rather
than raising, so the rest of the report can still be generated.
"""

import requests

PDBE_MAPPING_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
INTERPRO_PFAM_URL = "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{accession}"
REQUEST_TIMEOUT = 15


def get_uniprot_mappings(pdb_id: str) -> dict:
    """
    Return {chain_id: [mapping_segment, ...]} where each mapping_segment is:
        {
            "uniprot_accession": str,
            "uniprot_name": str,
            "unp_start": int, "unp_end": int,     # UniProt sequence coords
            "pdb_start": int, "pdb_end": int,      # PDB author residue numbers
        }
    A chain can map to multiple UniProt accessions (rare, e.g. fusion
    constructs) or multiple discontiguous segments of the same accession.
    """
    pdb_id = pdb_id.strip().lower()
    url = PDBE_MAPPING_URL.format(pdb_id=pdb_id)

    try:
        resp = requests.get(url, timeout=REQUEST_TIMEOUT)
        resp.raise_for_status()
        data = resp.json()
    except Exception as exc:
        return {"_error": f"PDBe mapping request failed: {exc}"}

    entry = data.get(pdb_id)
    if not entry:
        return {"_error": f"No SIFTS mapping entry found for {pdb_id}"}

    uniprot_block = entry.get("UniProt", {})
    if not uniprot_block:
        return {"_error": f"No UniProt cross-references found for {pdb_id}"}

    by_chain: dict = {}
    for accession, info in uniprot_block.items():
        name = info.get("name", "")
        for mapping in info.get("mappings", []):
            chain_id = mapping.get("chain_id") or mapping.get("struct_asym_id")
            start = mapping.get("start", {})
            end = mapping.get("end", {})
            segment = {
                "uniprot_accession": accession,
                "uniprot_name": name,
                "unp_start": mapping.get("unp_start"),
                "unp_end": mapping.get("unp_end"),
                "pdb_start": start.get("author_residue_number"),
                "pdb_end": end.get("author_residue_number"),
            }
            if None in segment.values():
                continue
            by_chain.setdefault(chain_id, []).append(segment)

    return by_chain


def get_pfam_domains(uniprot_accession: str) -> list:
    """
    Return a list of Pfam domain hits for a UniProt accession, in UniProt
    sequence coordinates:
        [{"pfam_id": "PF00069", "name": "Protein kinase domain",
          "unp_start": 12, "unp_end": 270}, ...]
    """
    url = INTERPRO_PFAM_URL.format(accession=uniprot_accession)

    try:
        resp = requests.get(url, params={"format": "json"}, timeout=REQUEST_TIMEOUT)
        resp.raise_for_status()
        data = resp.json()
    except Exception:
        return []

    domains = []
    for result in data.get("results", []):
        metadata = result.get("metadata", {})
        pfam_id = metadata.get("accession")
        name = metadata.get("name")
        for protein in result.get("proteins", []):
            for location in protein.get("entry_protein_locations", []) or []:
                for fragment in location.get("fragments", []) or []:
                    start, end = fragment.get("start"), fragment.get("end")
                    if start is None or end is None:
                        continue
                    domains.append({
                        "pfam_id": pfam_id,
                        "name": name,
                        "unp_start": start,
                        "unp_end": end,
                    })
    return domains


def _map_unp_range_to_pdb(unp_start: int, unp_end: int, segment: dict):
    """
    Translate a [unp_start, unp_end] UniProt-coordinate range into PDB
    residue numbering, using a single mapping segment's linear correspondence.
    Returns None if the domain doesn't overlap this segment's covered range.
    """
    offset = segment["pdb_start"] - segment["unp_start"]
    overlap_start = max(unp_start, segment["unp_start"])
    overlap_end = min(unp_end, segment["unp_end"])
    if overlap_start > overlap_end:
        return None
    return overlap_start + offset, overlap_end + offset


def annotate_domains(pdb_id: str) -> dict:
    """
    Full pipeline: PDB ID -> per-chain UniProt mapping -> Pfam domains ->
    domains reported in PDB chain/residue coordinates.

    Returns:
        {
            "available": bool,
            "reason": str,               # present only if available is False
            "chains": {
                "A": [
                    {"pfam_id": "PF00069", "name": "Protein kinase domain",
                     "pdb_start": 15, "pdb_end": 270,
                     "uniprot_accession": "P00519"},
                    ...
                ]
            }
        }
    """
    mappings = get_uniprot_mappings(pdb_id)
    if "_error" in mappings:
        return {"available": False, "reason": mappings["_error"], "chains": {}}

    chains_out: dict = {}
    accession_cache: dict = {}

    for chain_id, segments in mappings.items():
        chain_domains = []
        for segment in segments:
            accession = segment["uniprot_accession"]
            if accession not in accession_cache:
                accession_cache[accession] = get_pfam_domains(accession)
            pfam_hits = accession_cache[accession]

            for hit in pfam_hits:
                mapped = _map_unp_range_to_pdb(hit["unp_start"], hit["unp_end"], segment)
                if mapped is None:
                    continue
                pdb_start, pdb_end = mapped
                chain_domains.append({
                    "pfam_id": hit["pfam_id"],
                    "name": hit["name"],
                    "pdb_start": pdb_start,
                    "pdb_end": pdb_end,
                    "uniprot_accession": accession,
                })

        if chain_domains:
            # de-duplicate identical hits from overlapping segments, keep order
            seen = set()
            deduped = []
            for d in chain_domains:
                key = (d["pfam_id"], d["pdb_start"], d["pdb_end"])
                if key not in seen:
                    seen.add(key)
                    deduped.append(d)
            chains_out[chain_id] = sorted(deduped, key=lambda d: d["pdb_start"])

    return {"available": True, "chains": chains_out}
