#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Populate demo chain JSON files with small synthetic feature and annotation tracks.

This script is only intended for the public demo dataset. It makes the included
example JSON files look like a completed pipeline output without requiring
external tools such as mkdssp, IUPred2A, ChimeraX, or online UniProt/SIFTS calls.
It should not be used for real scientific analysis.
"""

import argparse
import csv
import json
from pathlib import Path


DEMO_MAPPING = {
    ("demo_1abc", "A"): ("1ABC", "P00001", "DEMO_PROT_A"),
    ("demo_1abc", "B"): ("1ABC", "P00002", "DEMO_PROT_B"),
    ("demo_2xyz", "A"): ("2XYZ", "P00003", "DEMO_PROT_C"),
    ("demo_2xyz", "C"): ("2XYZ", "P00004", "DEMO_PROT_D"),
}


def read_json(path: Path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def write_json(path: Path, data):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def write_csv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def populate_chain_json(path: Path):
    data = read_json(path)
    entry_id = data.get("entry", {}).get("entry_id", path.parent.parent.name)
    chain_id = data.get("chain", {}).get("chain_id", path.stem)
    sequence = data.get("chain", {}).get("sequence", "")
    n = len(sequence)

    pdb_id, accession, uniprot_id = DEMO_MAPPING.get(
        (entry_id, chain_id),
        (data.get("entry", {}).get("pdb_id", "DEMO"), "P00000", "DEMO_PROTEIN"),
    )

    data["uniprot_mapping"] = {
        "status": "mapped",
        "quality": "HIGH",
        "method": "demo_sifts_exact_match",
        "final_mapping_method": "demo_sifts_exact_match",
        "uniprot_accession": accession,
        "uniprot_id": uniprot_id,
        "coverage": 1.0,
        "identity": 1.0,
        "chain_index_to_uniprot": {str(i): i for i in range(1, n + 1)},
    }
    data["external_annotations"] = {
        "status": "done",
        "source": "synthetic_demo",
        "uniprot_accession": accession,
    }

    q_values = [
        {"chain_index": i, "qscore": round(0.65 + 0.01 * i, 3)}
        for i in range(1, n + 1)
    ]
    secondary_structure = [{"start": 1, "end": min(n, 2), "type": "H"}]
    if n > 3:
        secondary_structure.append({"start": 3, "end": n, "type": "E"})

    tracks = data.setdefault("tracks", {})
    tracks.update(
        {
            "secondary_structure": {
                "status": "done",
                "items": secondary_structure,
                "method": "demo_assigned",
            },
            "buried_residues": {
                "status": "done",
                "items": [{"chain_index": 2, "rsa": 0.12}] if n >= 2 else [],
                "threshold": 0.20,
            },
            "rsa": {
                "status": "done",
                "values": [
                    {"chain_index": i, "rsa": round(0.25 + 0.05 * (i % 3), 3)}
                    for i in range(1, n + 1)
                ],
            },
            "hydropathy": {
                "status": "done",
                "values": [
                    {"chain_index": i, "value": round(0.1 + 0.02 * i, 3)}
                    for i in range(1, n + 1)
                ],
            },
            "disorder": {
                "status": "done",
                "values": [
                    {"chain_index": i, "score": round(0.15 + 0.03 * (i % 2), 3)}
                    for i in range(1, n + 1)
                ],
                "segments": [],
            },
            "atomic_clashes": {"status": "done", "items": []},
            "bond_outliers": {"status": "done", "items": []},
            "cis_peptide": {"status": "done", "items": []},
            "qscore": {
                "status": "done",
                "values": q_values,
                "mean_qscore": round(sum(v["qscore"] for v in q_values) / n, 3) if n else 0,
            },
            "uniprot_features": {
                "status": "done",
                "items": [
                    {
                        "type": "chain",
                        "start": 1,
                        "end": n,
                        "description": "Demo UniProt feature mapped to chain coordinates",
                    }
                ],
            },
            "pfam": {
                "status": "done",
                "items": [{"accession": "PF00000", "name": "Demo_domain", "start": 1, "end": n}],
            },
            "interpro": {
                "status": "done",
                "items": [{"accession": "IPR000000", "name": "Demo InterPro family", "start": 1, "end": n}],
            },
            "domain_family": {"status": "empty", "items": []},
            "active_site_site": {"status": "empty", "items": []},
            "binding_site": {"status": "empty", "items": []},
            "molecule_processing": {"status": "empty", "items": []},
            "ptm_site": {"status": "empty", "items": []},
        }
    )

    data.setdefault("status", {}).update(
        {
            "local_tracks": "done",
            "uniprot_mapping": "done",
            "external_annotations": "done",
            "frontend_tracks": "done",
        }
    )

    write_json(path, data)

    local_row = {
        "entry_id": entry_id,
        "pdb_id": pdb_id,
        "chain_id": chain_id,
        "chain_length": n,
        "secondary_structure_status": "done",
        "rsa_status": "done",
        "disorder_status": "done",
        "qscore_status": "done",
        "mean_qscore": tracks["qscore"]["mean_qscore"],
    }
    annotation_row = {
        "entry_id": entry_id,
        "pdb_id": pdb_id,
        "chain_id": chain_id,
        "uniprot_accession": accession,
        "uniprot_id": uniprot_id,
        "mapping_coverage": 1.0,
        "mapping_identity": 1.0,
        "mapping_quality": "HIGH",
        "mapping_method": "demo_sifts_exact_match",
    }
    group_row = {
        "group_id": f"uniprot_{accession}",
        "group_type": "uniprot_group",
        "entry_id": entry_id,
        "pdb_id": pdb_id,
        "chain_id": chain_id,
        "representative_chain": "yes",
        "group_size": 1,
        "group_rule": "same UniProt accession",
    }
    return local_row, annotation_row, group_row


def main():
    parser = argparse.ArgumentParser(
        description="Populate synthetic feature and annotation tracks for the included demo JSON files."
    )
    parser.add_argument("--json-root", default="example_outputs/json")
    parser.add_argument("--feature-out", default="example_outputs/features/local_tracks_summary.csv")
    parser.add_argument("--annotation-out", default="example_outputs/annotation/uniprot_sifts_mapping.csv")
    parser.add_argument("--group-out", default="example_outputs/annotation/annotation_groups.csv")
    args = parser.parse_args()

    json_root = Path(args.json_root)
    local_rows = []
    annotation_rows = []
    group_rows = []

    for path in sorted(json_root.glob("*/chains/*.json")):
        local_row, annotation_row, group_row = populate_chain_json(path)
        local_rows.append(local_row)
        annotation_rows.append(annotation_row)
        group_rows.append(group_row)

    write_csv(Path(args.feature_out), local_rows)
    write_csv(Path(args.annotation_out), annotation_rows)
    write_csv(Path(args.group_out), group_rows)

    print(f"Updated demo chain JSON files: {len(local_rows)}")
    print(f"Feature summary: {args.feature_out}")
    print(f"Annotation summary: {args.annotation_out}")
    print(f"Annotation groups: {args.group_out}")


if __name__ == "__main__":
    main()
