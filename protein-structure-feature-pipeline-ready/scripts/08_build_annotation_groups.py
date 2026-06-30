#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import json
import os
from pathlib import Path
from datetime import datetime
from collections import Counter


DOMAIN_FAMILY_TYPES = {
    "Domain",
    "Repeat",
    "Region",
    "Motif",
    "Coiled coil",
    "Compositional bias",
    "Zinc finger",
    "Topological domain",
    "Transmembrane",
    "Intramembrane",
}

ACTIVE_SITE_TYPES = {
    "Active site",
    "Site",
}

BINDING_SITE_TYPES = {
    "Binding site",
    "Nucleotide binding",
    "DNA binding",
    "Calcium binding",
    "Metal binding",
}

MOLECULE_PROCESSING_TYPES = {
    "Initiator methionine",
    "Signal peptide",
    "Transit peptide",
    "Propeptide",
    "Chain",
    "Peptide",
}

PTM_SITE_TYPES = {
    "Modified residue",
    "Lipidation",
    "Glycosylation",
    "Disulfide bond",
    "Cross-link",
}


def now_iso():
    return datetime.now().isoformat(timespec="seconds")


def read_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def to_builtin_for_json(obj):
    if obj is None or isinstance(obj, (str, int, float, bool)):
        return obj
    if isinstance(obj, dict):
        return {str(to_builtin_for_json(k)): to_builtin_for_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_builtin_for_json(x) for x in obj]
    if isinstance(obj, set):
        return [to_builtin_for_json(x) for x in sorted(obj, key=lambda x: str(x))]
    if isinstance(obj, Path):
        return str(obj)
    if hasattr(obj, "item"):
        try:
            return to_builtin_for_json(obj.item())
        except Exception:
            pass
    if hasattr(obj, "tolist"):
        try:
            return to_builtin_for_json(obj.tolist())
        except Exception:
            pass
    return str(obj)


def write_json_atomic(path, data):
    tmp = str(path) + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(to_builtin_for_json(data), f, ensure_ascii=False, indent=2)
    os.replace(tmp, path)


def safe_items(track_obj):
    items = track_obj.get("items", [])
    return items if isinstance(items, list) else []


def clone_item(item, group_name, source_section):
    x = dict(item)
    x["annotation_group"] = group_name
    x["source_section"] = source_section
    return x


def make_track(items, source, external_status):
    if external_status == "unreliable_mapping":
        return {
            "status": "unreliable_mapping",
            "type": "features",
            "source": "not generated because final UniProt mapping is unreliable",
            "coordinate_system": "chain_index_1_based",
            "items": [],
            "item_count": 0,
        }

    status = "done" if items else "empty"
    return {
        "status": status,
        "type": "features",
        "source": source,
        "coordinate_system": "chain_index_1_based",
        "items": items,
        "item_count": len(items),
    }


def build_groups_for_chain_json(d):
    tracks = d.setdefault("tracks", {})
    ext = d.get("external_annotations", {})
    external_status = ext.get("status", "")

    uniprot_features = safe_items(tracks.get("uniprot_features", {}))
    pfam_items = safe_items(tracks.get("pfam", {}))
    interpro_items = safe_items(tracks.get("interpro", {}))

    domain_family = []
    active_site_site = []
    binding_site = []
    molecule_processing = []
    ptm_site = []

    for item in uniprot_features:
        ftype = item.get("feature_type", "")

        if ftype in DOMAIN_FAMILY_TYPES:
            domain_family.append(clone_item(item, "domain_family", "uniprot_features"))

        if ftype in ACTIVE_SITE_TYPES:
            active_site_site.append(clone_item(item, "active_site_site", "uniprot_features"))

        if ftype in BINDING_SITE_TYPES:
            binding_site.append(clone_item(item, "binding_site", "uniprot_features"))

        if ftype in MOLECULE_PROCESSING_TYPES:
            molecule_processing.append(clone_item(item, "molecule_processing", "uniprot_features"))

        if ftype in PTM_SITE_TYPES:
            ptm_site.append(clone_item(item, "ptm_site", "uniprot_features"))

    for item in pfam_items:
        domain_family.append(clone_item(item, "domain_family", "pfam"))

    for item in interpro_items:
        domain_family.append(clone_item(item, "domain_family", "interpro"))

    tracks["domain_family"] = make_track(
        domain_family,
        "UniProt domain-like features + SIFTS Pfam + SIFTS InterPro",
        external_status,
    )

    tracks["active_site_site"] = make_track(
        active_site_site,
        "UniProt Active site / Site features",
        external_status,
    )

    tracks["binding_site"] = make_track(
        binding_site,
        "UniProt binding-related features",
        external_status,
    )

    tracks["molecule_processing"] = make_track(
        molecule_processing,
        "UniProt molecule processing features",
        external_status,
    )

    tracks["ptm_site"] = make_track(
        ptm_site,
        "UniProt PTM / modification features",
        external_status,
    )

    d["tracks"] = tracks

    d.setdefault("external_annotations", {})
    d["external_annotations"]["derived_track_counts"] = {
        "domain_family": tracks["domain_family"]["item_count"],
        "active_site_site": tracks["active_site_site"]["item_count"],
        "binding_site": tracks["binding_site"]["item_count"],
        "molecule_processing": tracks["molecule_processing"]["item_count"],
        "ptm_site": tracks["ptm_site"]["item_count"],
    }

    d.setdefault("status", {})
    d["status"]["domain_family"] = tracks["domain_family"]["status"]
    d["status"]["active_site_site"] = tracks["active_site_site"]["status"]
    d["status"]["binding_site"] = tracks["binding_site"]["status"]
    d["status"]["molecule_processing"] = tracks["molecule_processing"]["status"]
    d["status"]["ptm_site"] = tracks["ptm_site"]["status"]

    d.setdefault("logs", [])
    d["logs"].append({
        "time": now_iso(),
        "step": "build_annotation_groups",
        "status": "done",
        "message": "Derived grouped annotation tracks from existing uniprot_features / pfam / interpro tracks.",
    })

    return d


def iter_chain_jsons(json_root, manifest=None, limit_entries=0):
    json_root = Path(json_root)

    if manifest:
        seen = set()
        count_entries = 0
        with open(manifest, newline="", encoding="utf-8", errors="ignore") as f:
            reader = csv.DictReader(f)
            for r in reader:
                entry_id = r.get("entry_id", "").strip()
                if not entry_id or entry_id in seen:
                    continue
                seen.add(entry_id)

                has_structure = r.get("has_structure", "").lower()
                skeleton_status = r.get("skeleton_status", "").lower()

                if has_structure and has_structure not in {"yes", "true", "1"}:
                    continue
                if skeleton_status and skeleton_status != "done":
                    continue

                count_entries += 1
                if limit_entries and count_entries > limit_entries:
                    break

                chain_dir = json_root / entry_id / "chains"
                if chain_dir.is_dir():
                    for jf in sorted(chain_dir.glob("*.json")):
                        yield jf
    else:
        entry_dirs = sorted([p for p in json_root.iterdir() if p.is_dir()])
        if limit_entries:
            entry_dirs = entry_dirs[:limit_entries]
        for entry_dir in entry_dirs:
            chain_dir = entry_dir / "chains"
            if chain_dir.is_dir():
                for jf in sorted(chain_dir.glob("*.json")):
                    yield jf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json-root", default="example_outputs/json")
    ap.add_argument("--manifest", default="example_outputs/manifest/manifest_skeleton.csv")
    ap.add_argument("--summary-out", default="example_outputs/qc/annotation_groups_summary.csv")
    ap.add_argument("--limit-entries", type=int, default=0)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    summary_out = Path(args.summary_out)
    summary_out.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    status_counters = {
        "domain_family": Counter(),
        "active_site_site": Counter(),
        "binding_site": Counter(),
        "molecule_processing": Counter(),
        "ptm_site": Counter(),
    }

    for jf in iter_chain_jsons(args.json_root, manifest=args.manifest, limit_entries=args.limit_entries):
        try:
            d = read_json(jf)
            d = build_groups_for_chain_json(d)

            if not args.dry_run:
                write_json_atomic(jf, d)

            tracks = d.get("tracks", {})
            row = {
                "entry_id": d.get("entry", {}).get("entry_id", jf.parent.parent.name),
                "pdb_id": d.get("entry", {}).get("pdb_id", ""),
                "chain_id": d.get("chain", {}).get("chain_id", jf.stem),
                "json_path": str(jf),
                "mapping_quality": d.get("uniprot_mapping", {}).get("quality", ""),
                "mapping_method": d.get("uniprot_mapping", {}).get("final_mapping_method", d.get("uniprot_mapping", {}).get("method", "")),
                "external_annotation_status": d.get("external_annotations", {}).get("status", ""),

                "domain_family_status": tracks.get("domain_family", {}).get("status", ""),
                "domain_family_count": tracks.get("domain_family", {}).get("item_count", 0),

                "active_site_site_status": tracks.get("active_site_site", {}).get("status", ""),
                "active_site_site_count": tracks.get("active_site_site", {}).get("item_count", 0),

                "binding_site_status": tracks.get("binding_site", {}).get("status", ""),
                "binding_site_count": tracks.get("binding_site", {}).get("item_count", 0),

                "molecule_processing_status": tracks.get("molecule_processing", {}).get("status", ""),
                "molecule_processing_count": tracks.get("molecule_processing", {}).get("item_count", 0),

                "ptm_site_status": tracks.get("ptm_site", {}).get("status", ""),
                "ptm_site_count": tracks.get("ptm_site", {}).get("item_count", 0),
                "error": "",
            }
            rows.append(row)

            for k in status_counters:
                status_counters[k][row[f"{k}_status"]] += 1

        except Exception as e:
            rows.append({
                "entry_id": jf.parent.parent.name,
                "pdb_id": "",
                "chain_id": jf.stem,
                "json_path": str(jf),
                "mapping_quality": "",
                "mapping_method": "",
                "external_annotation_status": "failed",
                "domain_family_status": "failed",
                "domain_family_count": 0,
                "active_site_site_status": "failed",
                "active_site_site_count": 0,
                "binding_site_status": "failed",
                "binding_site_count": 0,
                "molecule_processing_status": "failed",
                "molecule_processing_count": 0,
                "ptm_site_status": "failed",
                "ptm_site_count": 0,
                "error": repr(e),
            })
            print("[ERROR]", jf, repr(e))

    fields = [
        "entry_id", "pdb_id", "chain_id", "json_path",
        "mapping_quality", "mapping_method", "external_annotation_status",
        "domain_family_status", "domain_family_count",
        "active_site_site_status", "active_site_site_count",
        "binding_site_status", "binding_site_count",
        "molecule_processing_status", "molecule_processing_count",
        "ptm_site_status", "ptm_site_count",
        "error",
    ]

    with open(summary_out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    print("chains processed:", len(rows))
    print("dry_run:", args.dry_run)
    print("summary:", summary_out)

    for k, c in status_counters.items():
        print(k, dict(c))


if __name__ == "__main__":
    main()