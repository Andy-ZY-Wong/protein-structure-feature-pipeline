#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import json
from pathlib import Path
from collections import Counter


LOCAL_TRACKS = [
    "secondary_structure",
    "buried_residues",
    "rsa",
    "hydropathy",
    "disorder",
    "atomic_clashes",
    "bond_outliers",
    "cis_peptide",
    "qscore",
]

ANNOTATION_TRACKS = [
    "uniprot_features",
    "pfam",
    "interpro",
    "domain_family",
    "active_site_site",
    "binding_site",
    "molecule_processing",
    "ptm_site",
]


def read_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def safe_len_items(obj):
    """Return the number of records stored in a track object."""
    if "items" in obj and isinstance(obj.get("items"), list):
        return len(obj.get("items", []))
    if "values" in obj and isinstance(obj.get("values"), list):
        return len(obj.get("values", []))
    return obj.get("item_count", obj.get("mapped_count", ""))


def get_track_obj(d, name):
    tracks = d.get("tracks", {})
    if name in tracks:
        return tracks.get(name, {})

    # Compatibility with older JSON files where some fields may be at the top level
    if name in d:
        return d.get(name, {})

    # Compatibility with status fields that only contain a status string
    status = d.get("status", {}).get(name, "")
    if status:
        return {"status": status}

    return {}


def get_track_status(d, name):
    obj = get_track_obj(d, name)
    status = obj.get("status", "")
    if not status:
        status = d.get("status", {}).get(name, "")
    return status


def get_track_count(d, name):
    obj = get_track_obj(d, name)
    return safe_len_items(obj)


def is_bad_local_track(track_name, status, count):
    """
    Local computed parameters should generally be done.
    skipped / failed / partial / missing should all be recorded.
    If qscore is done but count=0, record it as well.
    """
    if not status:
        return True

    if status in {"failed", "error", "missing", "skipped"}:
        return True

    if track_name == "qscore":
        try:
            if status == "done" and int(float(count or 0)) == 0:
                return True
        except Exception:
            pass

    if status not in {"done", "partial"}:
        return True

    return False


def is_bad_annotation_track(track_name, status, mapping_quality, external_status):
    """
    Definition of annotation issues:
    1. When mapping is unreliable, annotation=unreliable_mapping is reasonable, but it is still exported to annotation_problem CSV,
       because such chains cannot reliably display external annotations.
    2. When mapping is reliable, annotation tracks should not be failed/missing/unreliable_mapping.
    3. empty is not necessarily an error; it means no annotation of that type is present.
    """
    reliable = mapping_quality in {"HIGH", "ACCEPTABLE", "PARTIAL_USABLE"}

    if not status:
        return True

    if reliable:
        if status in {"failed", "error", "missing", "unreliable_mapping"}:
            return True
        return False

    # Unreliable mapping is recorded as an annotation issue, but not as a script error
    if mapping_quality in {"PROBLEM", "FAILED", ""}:
        return True

    return False


def iter_chain_jsons(json_root, manifest=None):
    json_root = Path(json_root)

    if manifest:
        seen = set()
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

                chain_dir = json_root / entry_id / "chains"
                if chain_dir.is_dir():
                    for jf in sorted(chain_dir.glob("*.json")):
                        yield jf
    else:
        for entry_dir in sorted(json_root.iterdir()):
            chain_dir = entry_dir / "chains"
            if chain_dir.is_dir():
                for jf in sorted(chain_dir.glob("*.json")):
                    yield jf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json-root", default="example_outputs/json")
    ap.add_argument("--manifest", default="example_outputs/manifest/manifest_skeleton.csv")
    ap.add_argument("--summary-out", default="example_outputs/qc/all_chain_qc_summary.csv")
    ap.add_argument("--problem-out", default="example_outputs/qc/all_chain_qc_problem_chains.csv")
    ap.add_argument("--annotation-problem-out", default="example_outputs/qc/annotation_problem_chains.csv")
    ap.add_argument("--local-problem-out", default="example_outputs/qc/local_parameter_problem_chains.csv")
    args = ap.parse_args()

    summary_rows = []
    problem_rows = []
    annotation_problem_rows = []
    local_problem_rows = []

    counters = Counter()
    track_status_counter = Counter()
    annotation_status_counter = Counter()

    for jf in iter_chain_jsons(args.json_root, args.manifest):
        try:
            d = read_json(jf)
        except Exception as e:
            row = {
                "entry_id": jf.parent.parent.name,
                "pdb_id": "",
                "chain_id": jf.stem,
                "json_path": str(jf),
                "problem_type": "json_read_error",
                "problem_detail": repr(e),
            }
            problem_rows.append(row)
            continue

        entry_id = d.get("entry", {}).get("entry_id", jf.parent.parent.name)
        pdb_id = d.get("entry", {}).get("pdb_id", "")
        chain_id = d.get("chain", {}).get("chain_id", jf.stem)
        chain_length = d.get("chain", {}).get("length", "")

        m = d.get("uniprot_mapping", {})
        external = d.get("external_annotations", {})

        mapping_status = m.get("status", "")
        mapping_quality = m.get("quality", "")
        mapping_method = m.get("final_mapping_method", m.get("method", ""))
        mapping_coverage = m.get("coverage", "")
        mapping_identity = m.get("identity", "")
        external_status = external.get("status", "")

        row = {
            "entry_id": entry_id,
            "pdb_id": pdb_id,
            "chain_id": chain_id,
            "chain_length": chain_length,
            "json_path": str(jf),

            "mapping_status": mapping_status,
            "mapping_quality": mapping_quality,
            "mapping_method": mapping_method,
            "mapping_coverage": mapping_coverage,
            "mapping_identity": mapping_identity,
            "external_annotation_status": external_status,
        }

        local_flags = []
        annotation_flags = []

        for t in LOCAL_TRACKS:
            st = get_track_status(d, t)
            cnt = get_track_count(d, t)
            row[f"{t}_status"] = st
            row[f"{t}_count"] = cnt

            track_status_counter[(t, st)] += 1

            if is_bad_local_track(t, st, cnt):
                local_flags.append(f"{t}:{st or 'missing'}")

        for t in ANNOTATION_TRACKS:
            st = get_track_status(d, t)
            cnt = get_track_count(d, t)
            row[f"{t}_status"] = st
            row[f"{t}_count"] = cnt

            annotation_status_counter[(t, st)] += 1

            if is_bad_annotation_track(t, st, mapping_quality, external_status):
                annotation_flags.append(f"{t}:{st or 'missing'}")

        row["local_parameter_problem"] = "yes" if local_flags else "no"
        row["local_parameter_problem_flags"] = ";".join(local_flags)

        row["annotation_problem"] = "yes" if annotation_flags else "no"
        row["annotation_problem_flags"] = ";".join(annotation_flags)

        row["any_problem"] = "yes" if local_flags or annotation_flags else "no"

        summary_rows.append(row)

        counters["chains"] += 1
        counters[f"mapping_quality::{mapping_quality}"] += 1
        counters[f"external::{external_status}"] += 1

        if local_flags:
            local_problem_rows.append(row)

        if annotation_flags:
            annotation_problem_rows.append(row)

        if local_flags or annotation_flags:
            problem_rows.append(row)

    fields = [
        "entry_id", "pdb_id", "chain_id", "chain_length", "json_path",
        "mapping_status", "mapping_quality", "mapping_method",
        "mapping_coverage", "mapping_identity",
        "external_annotation_status",
    ]

    for t in LOCAL_TRACKS:
        fields.extend([f"{t}_status", f"{t}_count"])

    for t in ANNOTATION_TRACKS:
        fields.extend([f"{t}_status", f"{t}_count"])

    fields.extend([
        "local_parameter_problem",
        "local_parameter_problem_flags",
        "annotation_problem",
        "annotation_problem_flags",
        "any_problem",
    ])

    for path, rows in [
        (args.summary_out, summary_rows),
        (args.problem_out, problem_rows),
        (args.annotation_problem_out, annotation_problem_rows),
        (args.local_problem_out, local_problem_rows),
    ]:
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)

    print("chains:", len(summary_rows))
    print("summary:", args.summary_out)
    print("all problems:", len(problem_rows), args.problem_out)
    print("annotation problems:", len(annotation_problem_rows), args.annotation_problem_out)
    print("local parameter problems:", len(local_problem_rows), args.local_problem_out)

    print("\nMapping / external counters:")
    for k, v in counters.most_common():
        print(k, v)

    print("\nLocal track status:")
    for (track, status), v in sorted(track_status_counter.items()):
        print(track, status, v)

    print("\nAnnotation track status:")
    for (track, status), v in sorted(annotation_status_counter.items()):
        print(track, status, v)


if __name__ == "__main__":
    main()