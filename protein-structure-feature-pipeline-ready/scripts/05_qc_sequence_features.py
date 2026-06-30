#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
05_qc_sequence_features.py

功能：
对已经计算完成的 sequence features 做全库 QC。

检查对象：
1. Secondary Structure
2. RSA
3. Buried residues
4. Hydropathy
5. Disorder
6. Atomic clashes
7. Bond outliers
8. Q-score

输出：
1. entry-level QC 汇总表
2. chain-level QC 汇总表
3. problem-level 详细问题表

不会重新计算任何 feature。
只读取已有 JSON 并检查数据是否合理。
"""

import argparse
import csv
import json
import math
from pathlib import Path
from collections import Counter, defaultdict


REQUIRED_FEATURES = [
    "secondary_structure",
    "rsa",
    "buried_residues",
    "hydropathy",
    "disorder",
    "atomic_clashes",
    "bond_outliers",
    "qscore",
]


def read_json(path: Path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def is_number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(float(x))


def add_problem(problems, entry_id, chain_id, feature, severity, problem_type, message, json_path):
    problems.append({
        "entry_id": entry_id,
        "chain_id": chain_id,
        "feature": feature,
        "severity": severity,
        "problem_type": problem_type,
        "message": message,
        "json_path": str(json_path),
    })


def get_feature(d, name):
    return d.get("computed_features", {}).get(name, {})


def check_positions(
    positions,
    L,
    problems,
    entry_id,
    chain_id,
    feature,
    json_path
):
    if not isinstance(positions, list):
        add_problem(
            problems, entry_id, chain_id, feature, "ERROR",
            "positions_not_list",
            "positions is not a list",
            json_path
        )
        return

    bad = [
        p for p in positions
        if not isinstance(p, int) or p < 1 or p > L
    ]

    if bad:
        add_problem(
            problems, entry_id, chain_id, feature, "ERROR",
            "position_out_of_range",
            f"{len(bad)} positions out of range 1..{L}; examples={bad[:10]}",
            json_path
        )


def rebuild_segments_from_labels(labels):
    segs = []
    cur = None
    start = None

    for i, label in enumerate(labels, start=1):
        if label not in {"H", "E"}:
            if cur is not None:
                segs.append({
                    "start": start,
                    "end": i - 1,
                    "label": cur
                })
            cur = None
            start = None
            continue

        if cur is None:
            cur = label
            start = i
        elif label != cur:
            segs.append({
                "start": start,
                "end": i - 1,
                "label": cur
            })
            cur = label
            start = i

    if cur is not None:
        segs.append({
            "start": start,
            "end": len(labels),
            "label": cur
        })

    return segs


def normalize_segments(items):
    out = []

    if not isinstance(items, list):
        return out

    for x in items:
        if not isinstance(x, dict):
            continue

        try:
            s = int(x.get("start"))
            e = int(x.get("end"))
            label = str(x.get("label"))
        except Exception:
            continue

        out.append({
            "start": s,
            "end": e,
            "label": label
        })

    return out


def rebuild_disorder_segments(values, threshold):
    segs = []
    start = None

    for i, v in enumerate(values, start=1):
        is_dis = is_number(v) and float(v) >= threshold

        if is_dis and start is None:
            start = i

        if not is_dis and start is not None:
            segs.append({
                "start": start,
                "end": i - 1,
                "label": "disordered"
            })
            start = None

    if start is not None:
        segs.append({
            "start": start,
            "end": len(values),
            "label": "disordered"
        })

    return segs


def check_secondary_structure(d, L, problems, entry_id, chain_id, json_path):
    obj = get_feature(d, "secondary_structure")
    status = obj.get("status", "missing")

    if status != "done":
        add_problem(
            problems, entry_id, chain_id, "secondary_structure", "WARNING",
            "feature_not_done",
            f"secondary_structure status={status}",
            json_path
        )
        return

    pp = obj.get("per_position", [])

    if len(pp) != L:
        add_problem(
            problems, entry_id, chain_id, "secondary_structure", "ERROR",
            "length_mismatch",
            f"per_position length {len(pp)} != chain length {L}",
            json_path
        )
        return

    bad_labels = sorted(set(x for x in pp if x not in {"H", "E", "C"}))
    if bad_labels:
        add_problem(
            problems, entry_id, chain_id, "secondary_structure", "ERROR",
            "bad_labels",
            f"invalid secondary structure labels: {bad_labels}",
            json_path
        )

    expected = rebuild_segments_from_labels(pp)
    observed = normalize_segments(obj.get("items", []))

    if expected != observed:
        add_problem(
            problems, entry_id, chain_id, "secondary_structure", "WARNING",
            "segments_inconsistent",
            f"items segments inconsistent with per_position; expected={len(expected)}, observed={len(observed)}",
            json_path
        )


def check_rsa(d, L, problems, entry_id, chain_id, json_path):
    obj = get_feature(d, "rsa")
    status = obj.get("status", "missing")

    if status != "done":
        add_problem(
            problems, entry_id, chain_id, "rsa", "ERROR",
            "feature_not_done",
            f"rsa status={status}",
            json_path
        )
        return

    values = obj.get("values", [])

    if len(values) != L:
        add_problem(
            problems, entry_id, chain_id, "rsa", "ERROR",
            "length_mismatch",
            f"rsa values length {len(values)} != chain length {L}",
            json_path
        )
        return

    bad_negative = []
    suspicious_high = []
    null_count = 0

    for i, v in enumerate(values, start=1):
        if v is None:
            null_count += 1
            continue

        if not is_number(v):
            add_problem(
                problems, entry_id, chain_id, "rsa", "ERROR",
                "bad_value",
                f"rsa value at position {i} is not numeric: {v}",
                json_path
            )
            continue

        v = float(v)

        if v < 0:
            bad_negative.append((i, v))

        if v > 2.0:
            suspicious_high.append((i, v))

    if bad_negative:
        add_problem(
            problems, entry_id, chain_id, "rsa", "ERROR",
            "negative_rsa",
            f"{len(bad_negative)} RSA values < 0; examples={bad_negative[:10]}",
            json_path
        )

    if suspicious_high:
        add_problem(
            problems, entry_id, chain_id, "rsa", "WARNING",
            "very_high_rsa",
            f"{len(suspicious_high)} RSA values > 2.0; examples={suspicious_high[:10]}",
            json_path
        )

    if null_count > 0:
        add_problem(
            problems, entry_id, chain_id, "rsa", "WARNING",
            "rsa_null_values",
            f"{null_count} RSA values are null",
            json_path
        )


def check_buried(d, L, problems, entry_id, chain_id, json_path):
    obj = get_feature(d, "buried_residues")
    status = obj.get("status", "missing")

    if status != "done":
        add_problem(
            problems, entry_id, chain_id, "buried_residues", "ERROR",
            "feature_not_done",
            f"buried_residues status={status}",
            json_path
        )
        return

    positions = obj.get("positions", [])
    check_positions(positions, L, problems, entry_id, chain_id, "buried_residues", json_path)

    rsa = get_feature(d, "rsa")
    if rsa.get("status") == "done":
        values = rsa.get("values", [])

        if len(values) == L:
            threshold = obj.get("source", {}).get("threshold", 0.2)
            try:
                threshold = float(threshold)
            except Exception:
                threshold = 0.2

            expected = sorted([
                i for i, v in enumerate(values, start=1)
                if is_number(v) and float(v) < threshold
            ])
            observed = sorted(positions)

            if expected != observed:
                add_problem(
                    problems, entry_id, chain_id, "buried_residues", "WARNING",
                    "buried_rsa_inconsistent",
                    f"buried positions inconsistent with rsa<threshold; expected={len(expected)}, observed={len(observed)}",
                    json_path
                )


def check_hydropathy(d, L, problems, entry_id, chain_id, json_path):
    obj = get_feature(d, "hydropathy")
    status = obj.get("status", "missing")

    if status != "done":
        add_problem(
            problems, entry_id, chain_id, "hydropathy", "ERROR",
            "feature_not_done",
            f"hydropathy status={status}",
            json_path
        )
        return

    values = obj.get("values", [])

    if len(values) != L:
        add_problem(
            problems, entry_id, chain_id, "hydropathy", "ERROR",
            "length_mismatch",
            f"hydropathy values length {len(values)} != chain length {L}",
            json_path
        )
        return

    bad = []
    out_range = []

    for i, v in enumerate(values, start=1):
        if not is_number(v):
            bad.append((i, v))
            continue

        v = float(v)
        if v < -4.6 or v > 4.6:
            out_range.append((i, v))

    if bad:
        add_problem(
            problems, entry_id, chain_id, "hydropathy", "ERROR",
            "bad_value",
            f"{len(bad)} non-numeric hydropathy values; examples={bad[:10]}",
            json_path
        )

    if out_range:
        add_problem(
            problems, entry_id, chain_id, "hydropathy", "ERROR",
            "out_of_range",
            f"{len(out_range)} hydropathy values outside [-4.6, 4.6]; examples={out_range[:10]}",
            json_path
        )


def check_disorder(d, L, problems, entry_id, chain_id, json_path):
    obj = get_feature(d, "disorder")
    status = obj.get("status", "missing")

    if status != "done":
        add_problem(
            problems, entry_id, chain_id, "disorder", "ERROR",
            "feature_not_done",
            f"disorder status={status}",
            json_path
        )
        return

    values = obj.get("values", [])

    if len(values) != L:
        add_problem(
            problems, entry_id, chain_id, "disorder", "ERROR",
            "length_mismatch",
            f"disorder values length {len(values)} != chain length {L}",
            json_path
        )
        return

    bad = []
    out_range = []

    for i, v in enumerate(values, start=1):
        if not is_number(v):
            bad.append((i, v))
            continue

        v = float(v)
        if v < 0 or v > 1:
            out_range.append((i, v))

    if bad:
        add_problem(
            problems, entry_id, chain_id, "disorder", "ERROR",
            "bad_value",
            f"{len(bad)} non-numeric disorder values; examples={bad[:10]}",
            json_path
        )

    if out_range:
        add_problem(
            problems, entry_id, chain_id, "disorder", "ERROR",
            "out_of_range",
            f"{len(out_range)} disorder values outside [0, 1]; examples={out_range[:10]}",
            json_path
        )

    threshold = obj.get("source", {}).get("threshold", 0.5)
    try:
        threshold = float(threshold)
    except Exception:
        threshold = 0.5

    expected = rebuild_disorder_segments(values, threshold)
    observed = normalize_segments(obj.get("segments", []))

    if expected != observed:
        add_problem(
            problems, entry_id, chain_id, "disorder", "WARNING",
            "segments_inconsistent",
            f"disorder segments inconsistent with values>=threshold; expected={len(expected)}, observed={len(observed)}",
            json_path
        )


def check_atomic_clashes(d, L, problems, entry_id, chain_id, json_path, clash_per_res_warn):
    obj = get_feature(d, "atomic_clashes")
    status = obj.get("status", "missing")

    if status != "done":
        add_problem(
            problems, entry_id, chain_id, "atomic_clashes", "ERROR",
            "feature_not_done",
            f"atomic_clashes status={status}",
            json_path
        )
        return

    positions = obj.get("positions", [])
    items = obj.get("items", [])

    check_positions(positions, L, problems, entry_id, chain_id, "atomic_clashes", json_path)

    if not isinstance(items, list):
        add_problem(
            problems, entry_id, chain_id, "atomic_clashes", "ERROR",
            "items_not_list",
            "atomic_clashes items is not a list",
            json_path
        )
        return

    item_pos = set()
    bad_items = 0

    for item in items:
        if not isinstance(item, dict):
            bad_items += 1
            continue

        pos = item.get("pos")
        if isinstance(pos, int):
            item_pos.add(pos)

        overlap = item.get("overlap")
        distance = item.get("distance")

        if overlap is not None and is_number(overlap) and float(overlap) < 0:
            bad_items += 1

        if distance is not None and is_number(distance) and float(distance) < 0:
            bad_items += 1

    if bad_items:
        add_problem(
            problems, entry_id, chain_id, "atomic_clashes", "WARNING",
            "suspicious_items",
            f"{bad_items} suspicious atomic clash items",
            json_path
        )

    pos_set = set(positions)
    if not pos_set.issubset(item_pos) and items:
        missing = sorted(pos_set - item_pos)
        add_problem(
            problems, entry_id, chain_id, "atomic_clashes", "WARNING",
            "positions_not_in_items",
            f"{len(missing)} positions not found in items; examples={missing[:10]}",
            json_path
        )

    clash_count = obj.get("clash_pair_count", len(items))
    if is_number(clash_count) and L > 0:
        ratio = float(clash_count) / L
        if ratio > clash_per_res_warn:
            add_problem(
                problems, entry_id, chain_id, "atomic_clashes", "WARNING",
                "high_clash_density",
                f"clash_pair_count/length={ratio:.3f} > {clash_per_res_warn}",
                json_path
            )


def check_bond_outliers(d, L, problems, entry_id, chain_id, json_path, bond_per_res_warn):
    obj = get_feature(d, "bond_outliers")
    status = obj.get("status", "missing")

    if status != "done":
        add_problem(
            problems, entry_id, chain_id, "bond_outliers", "ERROR",
            "feature_not_done",
            f"bond_outliers status={status}",
            json_path
        )
        return

    positions = obj.get("positions", [])
    items = obj.get("items", [])

    check_positions(positions, L, problems, entry_id, chain_id, "bond_outliers", json_path)

    if not isinstance(items, list):
        add_problem(
            problems, entry_id, chain_id, "bond_outliers", "ERROR",
            "items_not_list",
            "bond_outliers items is not a list",
            json_path
        )
        return

    bad_items = 0

    for item in items:
        if not isinstance(item, dict):
            bad_items += 1
            continue

        pos = item.get("pos")
        if pos is not None and (not isinstance(pos, int) or pos < 1 or pos > L):
            bad_items += 1

        observed = item.get("observed")
        if observed is not None and is_number(observed) and float(observed) <= 0:
            bad_items += 1

    if bad_items:
        add_problem(
            problems, entry_id, chain_id, "bond_outliers", "WARNING",
            "suspicious_items",
            f"{bad_items} suspicious bond outlier items",
            json_path
        )

    outlier_count = obj.get("outlier_count", len(items))
    if is_number(outlier_count) and L > 0:
        ratio = float(outlier_count) / L
        if ratio > bond_per_res_warn:
            add_problem(
                problems, entry_id, chain_id, "bond_outliers", "WARNING",
                "high_bond_outlier_density",
                f"outlier_count/length={ratio:.3f} > {bond_per_res_warn}",
                json_path
            )


def check_qscore(d, L, problems, entry_id, chain_id, json_path, qscore_min_coverage):
    obj = get_feature(d, "qscore")
    status = obj.get("status", "missing")

    if status != "done":
        add_problem(
            problems, entry_id, chain_id, "qscore", "ERROR",
            "feature_not_done",
            f"qscore status={status}",
            json_path
        )
        return

    values = obj.get("values", [])

    if len(values) != L:
        add_problem(
            problems, entry_id, chain_id, "qscore", "ERROR",
            "length_mismatch",
            f"qscore values length {len(values)} != chain length {L}",
            json_path
        )
        return

    valid = []
    suspicious = []

    for i, v in enumerate(values, start=1):
        if v is None:
            continue

        if not is_number(v):
            add_problem(
                problems, entry_id, chain_id, "qscore", "ERROR",
                "bad_value",
                f"qscore value at position {i} is not numeric: {v}",
                json_path
            )
            continue

        v = float(v)
        valid.append(v)

        if v < -1.0 or v > 1.5:
            suspicious.append((i, v))

    mapped_count = obj.get("mapped_count", None)
    if mapped_count is not None and is_number(mapped_count):
        if int(mapped_count) != len(valid):
            add_problem(
                problems, entry_id, chain_id, "qscore", "WARNING",
                "mapped_count_inconsistent",
                f"mapped_count={mapped_count}, non-null values={len(valid)}",
                json_path
            )

    coverage = obj.get("coverage", None)
    actual_coverage = len(valid) / L if L > 0 else 0.0

    if coverage is not None and is_number(coverage):
        if abs(float(coverage) - actual_coverage) > 0.01:
            add_problem(
                problems, entry_id, chain_id, "qscore", "WARNING",
                "coverage_inconsistent",
                f"coverage={coverage}, actual={actual_coverage:.4f}",
                json_path
            )

    if actual_coverage < qscore_min_coverage:
        add_problem(
            problems, entry_id, chain_id, "qscore", "WARNING",
            "low_coverage",
            f"qscore coverage={actual_coverage:.4f} < {qscore_min_coverage}",
            json_path
        )

    if suspicious:
        add_problem(
            problems, entry_id, chain_id, "qscore", "WARNING",
            "suspicious_value_range",
            f"{len(suspicious)} qscore values outside [-1.0, 1.5]; examples={suspicious[:10]}",
            json_path
        )

    mean = obj.get("mean", None)
    if mean is not None and is_number(mean) and valid:
        actual_mean = sum(valid) / len(valid)
        if abs(float(mean) - actual_mean) > 0.01:
            add_problem(
                problems, entry_id, chain_id, "qscore", "WARNING",
                "mean_inconsistent",
                f"mean={mean}, actual={actual_mean:.4f}",
                json_path
            )


def qc_one_chain(
    chain_json_path,
    entry_id,
    args
):
    problems = []

    try:
        d = read_json(chain_json_path)
    except Exception as e:
        add_problem(
            problems, entry_id, chain_json_path.stem, "json", "ERROR",
            "json_read_error",
            repr(e),
            chain_json_path
        )
        return None, problems

    chain = d.get("chain", {})
    chain_id = str(chain.get("chain_id", chain_json_path.stem))

    try:
        L = int(chain.get("length"))
    except Exception:
        L = 0

    if L <= 0:
        add_problem(
            problems, entry_id, chain_id, "chain", "ERROR",
            "bad_chain_length",
            f"invalid chain length: {chain.get('length')}",
            chain_json_path
        )
        return d, problems

    seq = chain.get("sequence", "")
    if isinstance(seq, str) and seq and len(seq) != L:
        add_problem(
            problems, entry_id, chain_id, "chain", "WARNING",
            "sequence_length_mismatch",
            f"sequence length {len(seq)} != chain length {L}",
            chain_json_path
        )

    residue_mapping = d.get("residue_mapping", [])
    if isinstance(residue_mapping, list) and residue_mapping:
        if len(residue_mapping) != L:
            add_problem(
                problems, entry_id, chain_id, "chain", "WARNING",
                "residue_mapping_length_mismatch",
                f"residue_mapping length {len(residue_mapping)} != chain length {L}",
                chain_json_path
            )
    else:
        add_problem(
            problems, entry_id, chain_id, "chain", "WARNING",
            "missing_residue_mapping",
            "residue_mapping missing or empty",
            chain_json_path
        )

    # Required features
    computed = d.get("computed_features", {})
    if not isinstance(computed, dict):
        add_problem(
            problems, entry_id, chain_id, "computed_features", "ERROR",
            "computed_features_not_dict",
            "computed_features is not a dict",
            chain_json_path
        )
        return d, problems

    for feature in REQUIRED_FEATURES:
        if feature not in computed:
            add_problem(
                problems, entry_id, chain_id, feature, "ERROR",
                "feature_missing",
                f"{feature} missing from computed_features",
                chain_json_path
            )

    if "secondary_structure" in computed:
        check_secondary_structure(d, L, problems, entry_id, chain_id, chain_json_path)

    if "rsa" in computed:
        check_rsa(d, L, problems, entry_id, chain_id, chain_json_path)

    if "buried_residues" in computed:
        check_buried(d, L, problems, entry_id, chain_id, chain_json_path)

    if "hydropathy" in computed:
        check_hydropathy(d, L, problems, entry_id, chain_id, chain_json_path)

    if "disorder" in computed:
        check_disorder(d, L, problems, entry_id, chain_id, chain_json_path)

    if "atomic_clashes" in computed:
        check_atomic_clashes(
            d, L, problems, entry_id, chain_id, chain_json_path,
            clash_per_res_warn=args.clash_per_res_warn
        )

    if "bond_outliers" in computed:
        check_bond_outliers(
            d, L, problems, entry_id, chain_id, chain_json_path,
            bond_per_res_warn=args.bond_per_res_warn
        )

    if "qscore" in computed:
        check_qscore(
            d, L, problems, entry_id, chain_id, chain_json_path,
            qscore_min_coverage=args.qscore_min_coverage
        )

    return d, problems


def read_manifest_entries(path):
    rows = []

    with open(path, newline="", encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f)
        for row in reader:
            entry_id = row.get("entry_id", "").strip()
            if not entry_id:
                continue

            has_structure = row.get("has_structure", "").lower()
            skeleton_status = row.get("skeleton_status", "").lower()

            if has_structure and has_structure not in {"yes", "true", "1"}:
                continue

            if skeleton_status and skeleton_status != "done":
                continue

            rows.append(row)

    return rows


def write_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main():
    p = argparse.ArgumentParser(description="QC computed sequence features.")

    p.add_argument(
        "--manifest",
        default="example_outputs/manifest/manifest_skeleton.csv"
    )
    p.add_argument(
        "--json-root",
        default="example_outputs/json"
    )
    p.add_argument(
        "--qc-dir",
        default="example_outputs/qc"
    )
    p.add_argument(
        "--qscore-min-coverage",
        type=float,
        default=0.8,
        help="Q-score coverage below this threshold will be marked as WARNING."
    )
    p.add_argument(
        "--clash-per-res-warn",
        type=float,
        default=2.0,
        help="atomic clash_pair_count / chain_length above this threshold will be WARNING."
    )
    p.add_argument(
        "--bond-per-res-warn",
        type=float,
        default=0.2,
        help="bond_outlier_count / chain_length above this threshold will be WARNING."
    )
    p.add_argument(
        "--limit",
        type=int,
        default=0
    )

    args = p.parse_args()

    manifest = Path(args.manifest)
    json_root = Path(args.json_root)
    qc_dir = Path(args.qc_dir)

    entries = read_manifest_entries(manifest)

    if args.limit and args.limit > 0:
        entries = entries[:args.limit]

    print("========== Sequence Feature QC ==========")
    print(f"Manifest       : {manifest}")
    print(f"JSON root      : {json_root}")
    print(f"QC dir         : {qc_dir}")
    print(f"Entries        : {len(entries)}")
    print(f"Q-score min coverage: {args.qscore_min_coverage}")
    print("=========================================")

    all_problem_rows = []
    chain_rows = []
    entry_agg = defaultdict(lambda: {
        "entry_id": "",
        "cfdb_id": "",
        "emdb_id": "",
        "pdb_id": "",
        "structure_path": "",
        "qc_status": "PASS",
        "chain_count": 0,
        "error_count": 0,
        "warning_count": 0,
        "problem_chain_count": 0,
        "problem_chains": set(),
        "problem_features": Counter(),
        "problem_types": Counter(),
    })

    global_status = Counter()
    feature_status = defaultdict(Counter)

    for i, row in enumerate(entries, start=1):
        entry_id = row.get("entry_id", "").strip()
        entry_dir = json_root / entry_id
        chain_dir = entry_dir / "chains"

        print(f"[INFO] ({i}/{len(entries)}) QC entry: {entry_id}")

        agg = entry_agg[entry_id]
        agg["entry_id"] = entry_id
        agg["cfdb_id"] = row.get("cfdb_id", "")
        agg["emdb_id"] = row.get("emdb_id", "")
        agg["pdb_id"] = row.get("pdb_id", "")
        agg["structure_path"] = row.get("structure_path", "")

        if not chain_dir.is_dir():
            problem = {
                "entry_id": entry_id,
                "chain_id": "",
                "feature": "entry",
                "severity": "ERROR",
                "problem_type": "missing_chain_dir",
                "message": f"chain directory not found: {chain_dir}",
                "json_path": str(chain_dir),
            }
            all_problem_rows.append(problem)
            agg["error_count"] += 1
            agg["problem_features"]["entry"] += 1
            agg["problem_types"]["missing_chain_dir"] += 1
            continue

        chain_files = sorted(chain_dir.glob("*.json"))

        if not chain_files:
            problem = {
                "entry_id": entry_id,
                "chain_id": "",
                "feature": "entry",
                "severity": "ERROR",
                "problem_type": "no_chain_json",
                "message": f"no chain json found under {chain_dir}",
                "json_path": str(chain_dir),
            }
            all_problem_rows.append(problem)
            agg["error_count"] += 1
            agg["problem_features"]["entry"] += 1
            agg["problem_types"]["no_chain_json"] += 1
            continue

        for jf in chain_files:
            d, problems = qc_one_chain(jf, entry_id, args)

            chain_id = jf.stem
            L = ""

            if d:
                chain_id = str(d.get("chain", {}).get("chain_id", jf.stem))
                L = d.get("chain", {}).get("length", "")

                feats = d.get("computed_features", {})
                for feat in REQUIRED_FEATURES:
                    feature_status[feat][feats.get(feat, {}).get("status", "missing")] += 1

            error_count = sum(1 for x in problems if x["severity"] == "ERROR")
            warning_count = sum(1 for x in problems if x["severity"] == "WARNING")

            if error_count > 0:
                qc_status = "FAIL"
            elif warning_count > 0:
                qc_status = "WARNING"
            else:
                qc_status = "PASS"

            global_status[qc_status] += 1

            chain_rows.append({
                "entry_id": entry_id,
                "chain_id": chain_id,
                "chain_length": L,
                "qc_status": qc_status,
                "error_count": error_count,
                "warning_count": warning_count,
                "json_path": str(jf),
            })

            agg["chain_count"] += 1
            agg["error_count"] += error_count
            agg["warning_count"] += warning_count

            if problems:
                agg["problem_chains"].add(chain_id)

            for problem in problems:
                all_problem_rows.append(problem)
                agg["problem_features"][problem["feature"]] += 1
                agg["problem_types"][problem["problem_type"]] += 1

    entry_rows = []

    for entry_id, agg in sorted(entry_agg.items()):
        if agg["error_count"] > 0:
            qc_status = "FAIL"
        elif agg["warning_count"] > 0:
            qc_status = "WARNING"
        else:
            qc_status = "PASS"

        entry_rows.append({
            "entry_id": entry_id,
            "cfdb_id": agg["cfdb_id"],
            "emdb_id": agg["emdb_id"],
            "pdb_id": agg["pdb_id"],
            "structure_path": agg["structure_path"],
            "qc_status": qc_status,
            "chain_count": agg["chain_count"],
            "problem_chain_count": len(agg["problem_chains"]),
            "problem_chains": ",".join(sorted(agg["problem_chains"])),
            "error_count": agg["error_count"],
            "warning_count": agg["warning_count"],
            "problem_features": ";".join(f"{k}:{v}" for k, v in sorted(agg["problem_features"].items())),
            "problem_types": ";".join(f"{k}:{v}" for k, v in sorted(agg["problem_types"].items())),
        })

    entry_out = qc_dir / "sequence_tracks_qc_entries.csv"
    chain_out = qc_dir / "sequence_tracks_qc_chains.csv"
    problem_out = qc_dir / "sequence_tracks_qc_problems.csv"
    summary_out = qc_dir / "sequence_tracks_qc_summary.txt"

    write_csv(
        entry_out,
        entry_rows,
        [
            "entry_id",
            "cfdb_id",
            "emdb_id",
            "pdb_id",
            "structure_path",
            "qc_status",
            "chain_count",
            "problem_chain_count",
            "problem_chains",
            "error_count",
            "warning_count",
            "problem_features",
            "problem_types",
        ]
    )

    write_csv(
        chain_out,
        chain_rows,
        [
            "entry_id",
            "chain_id",
            "chain_length",
            "qc_status",
            "error_count",
            "warning_count",
            "json_path",
        ]
    )

    write_csv(
        problem_out,
        all_problem_rows,
        [
            "entry_id",
            "chain_id",
            "feature",
            "severity",
            "problem_type",
            "message",
            "json_path",
        ]
    )

    with open(summary_out, "w", encoding="utf-8") as f:
        f.write("========== Sequence Feature QC Summary ==========\n")
        f.write(f"entries checked: {len(entry_rows)}\n")
        f.write(f"chains checked: {len(chain_rows)}\n")
        f.write(f"problems total: {len(all_problem_rows)}\n\n")

        f.write("chain qc status:\n")
        for k, v in global_status.items():
            f.write(f"  {k}: {v}\n")

        f.write("\nfeature status:\n")
        for feat in REQUIRED_FEATURES:
            f.write(f"  {feat}: {dict(feature_status[feat])}\n")

        f.write("\nproblem severity:\n")
        sev = Counter(x["severity"] for x in all_problem_rows)
        for k, v in sev.items():
            f.write(f"  {k}: {v}\n")

        f.write("\nproblem types:\n")
        pt = Counter(x["problem_type"] for x in all_problem_rows)
        for k, v in pt.most_common():
            f.write(f"  {k}: {v}\n")

    print("========== QC Finished ==========")
    print(f"Entry QC table   : {entry_out}")
    print(f"Chain QC table   : {chain_out}")
    print(f"Problem table    : {problem_out}")
    print(f"Summary          : {summary_out}")
    print("=================================")


if __name__ == "__main__":
    main()