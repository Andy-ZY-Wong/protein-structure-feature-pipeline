#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Run emdbMatch on selected problem entries and classify whether the generated
toRCSB mapping is usable for downstream UniProt mapping repair.

Key changes compared with the previous runner:
1. Default emdbMatch script is changed to emdbMatch_relaxed.py.
2. The result status is no longer judged only by returncode.
3. The main QC target is toRCSB coverage against existing JSON chains.
4. A non-zero returncode can still be retained as usable if toRCSB records exist.
"""

import argparse
import csv
import shutil
import subprocess
from pathlib import Path
from collections import defaultdict, Counter


def read_manifest(path):
    out = {}
    with open(path, newline="", encoding="utf-8", errors="ignore") as f:
        for r in csv.DictReader(f):
            eid = r.get("entry_id", "").strip()
            if eid:
                out[eid] = r
    return out


def read_entries(path, entry_col="entry_id"):
    entries = []
    with open(path, newline="", encoding="utf-8-sig", errors="ignore") as f:
        reader = csv.DictReader(f)
        cols = reader.fieldnames or []

        if entry_col not in cols:
            if "entry_id" in cols:
                entry_col = "entry_id"
            else:
                raise ValueError(f"Cannot find entry column in {path}. Columns={cols}")

        for r in reader:
            eid = r.get(entry_col, "").strip()
            if eid:
                entries.append(eid)

    # Deduplicate while preserving order.
    return list(dict.fromkeys(entries))


def load_aux_index(path):
    idx = defaultdict(list)
    with open(path, newline="", encoding="utf-8", errors="ignore") as f:
        for r in csv.DictReader(f):
            pdb = r.get("pdb_id", "").strip().lower()
            if pdb:
                idx[pdb].append(r)
    return idx


def add_same_stem_siblings(path):
    """
    If a PDB id points to 250523/9397.cif, also collect sibling files such as:
    250523/9397.faa, 250523/9397.manAlign, 250523/9397.fullLen, etc.
    """
    path = Path(path)
    parent = path.parent
    stem = path.stem

    wanted = [
        "pdb", "cif", "fasta", "fa", "faa",
        "fullLen", "manAlign", "autoAlign",
        "toRCSB", "autoDel", "manDel",
    ]

    files = {}
    for suf in wanted:
        q = parent / f"{stem}.{suf}"
        if q.exists():
            files[suf] = q

    return files


def candidate_groups_for_pdb(pdb_id, aux_index):
    """
    Return candidate batch groups for one PDB id.

    Each group is keyed by the same stem path without suffix, e.g.
    example_data/raw/batch/example_case
    """
    groups = {}

    for r in aux_index.get(pdb_id.lower(), []):
        p = Path(r.get("path", ""))
        if not p.exists():
            continue

        stem_key = str(p.with_suffix(""))
        group = groups.setdefault(stem_key, {})
        group[r.get("suffix", p.suffix.lstrip("."))] = p

        sibs = add_same_stem_siblings(p)
        for suf, q in sibs.items():
            group[suf] = q

    return groups


def find_web_file(entry_id, web_roots, suffix):
    for root in web_roots:
        p = Path(root) / f"{entry_id}.{suffix}"
        if p.exists():
            return p
    return None


def copy_to_work(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return dst


def count_fasta_records(path):
    path = Path(path)
    if not path.exists():
        return 0

    n = 0
    try:
        with open(path, encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith(">"):
                    n += 1
    except Exception:
        return 0

    return n


def parse_torcsb_local_chains(path):
    """
    Parse .toRCSB produced by emdbMatch.py.

    Header format from emdbMatch.py:
        >{t2c[t]}_{n}

    Example:
        >6nm5_1_3D_0

    Interpreted as:
        local chain = 0
        RCSB/template chain = 3D
    """
    path = Path(path)
    out = {}

    if not path.exists():
        return out

    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                continue

            h = line[1:].split()[0]
            parts = h.split("_")

            if len(parts) < 4:
                continue

            local_chain = parts[-1]
            rcsb_chain = parts[-2]

            out[local_chain] = {
                "rcsb_chain": rcsb_chain,
                "header": h,
            }

    return out


def get_json_chains(json_root, entry_id):
    chain_dir = Path(json_root) / entry_id / "chains"

    if not chain_dir.is_dir():
        return set()

    return {p.stem for p in chain_dir.glob("*.json")}


def sort_chain_ids(chains):
    return sorted(chains, key=lambda x: (len(x), x))


def classify_emdbmatch_result(returncode, out_toRCSB, out_fasta, out_autoDel, json_root, entry_id):
    """
    Classify the emdbMatch result by toRCSB coverage, not merely by returncode.

    This is safer for our mapping-repair goal because .toRCSB is the useful
    local_chain -> RCSB_chain bridge. fasta/autoDel generation can fail in
    the original emdbMatch for descD[seq] issues, but .toRCSB may still be usable.
    """
    torcsb_map = parse_torcsb_local_chains(out_toRCSB)
    torcsb_chains = set(torcsb_map)
    json_chains = get_json_chains(json_root, entry_id)

    both = json_chains & torcsb_chains
    only_json = json_chains - torcsb_chains
    only_torcsb = torcsb_chains - json_chains

    toRCSB_records = len(torcsb_chains)
    fasta_records = count_fasta_records(out_fasta)
    autoDel_records = count_fasta_records(out_autoDel)

    json_chain_count = len(json_chains)
    torcsb_chain_count = len(torcsb_chains)
    torcsb_json_overlap = len(both)
    torcsb_json_coverage = (
        torcsb_json_overlap / json_chain_count
        if json_chain_count else 0.0
    )

    if toRCSB_records == 0:
        status = "failed_no_toRCSB"
        reason = "no usable toRCSB records"
    elif json_chain_count > 0 and torcsb_json_overlap == json_chain_count:
        if returncode == 0:
            status = "done_full_toRCSB"
            reason = ""
        else:
            status = "partial_toRCSB_usable_full_overlap_returncode_nonzero"
            reason = "returncode non-zero but all JSON chains are covered by toRCSB"
    elif json_chain_count > 0 and torcsb_json_overlap > 0:
        status = "partial_toRCSB_usable"
        reason = "some JSON chains are covered by toRCSB; missing chains should not be repaired by emdbMatch"
    elif json_chain_count == 0 and toRCSB_records > 0:
        if returncode == 0:
            status = "done_toRCSB_no_json_check"
            reason = "toRCSB exists but JSON chain directory was not found or empty"
        else:
            status = "partial_toRCSB_usable_no_json_check_returncode_nonzero"
            reason = "returncode non-zero but toRCSB exists; JSON chain directory was not found or empty"
    else:
        status = "failed_no_json_overlap"
        reason = "toRCSB exists but none of the local chains overlap JSON chains"

    return {
        "status": status,
        "returncode": returncode,
        "toRCSB_records": toRCSB_records,
        "fasta_records": fasta_records,
        "autoDel_records": autoDel_records,
        "json_chain_count": json_chain_count,
        "torcsb_chain_count": torcsb_chain_count,
        "torcsb_json_overlap": torcsb_json_overlap,
        "torcsb_json_coverage": f"{torcsb_json_coverage:.6f}",
        "missing_json_chain_count": len(only_json),
        "extra_torcsb_chain_count": len(only_torcsb),
        "missing_json_chains": ",".join(sort_chain_ids(only_json)),
        "extra_torcsb_chains": ",".join(sort_chain_ids(only_torcsb)),
        "reason": reason,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--entries-csv", required=True)
    ap.add_argument("--manifest", default="example_outputs/manifest/manifest_skeleton.csv")
    ap.add_argument("--aux-index", default="example_outputs/qc/batch_aux_file_index.csv")
    ap.add_argument("--web-roots", default="example_data/raw/web")
    ap.add_argument("--emdbmatch", default="scripts/emdbMatch_relaxed.py")
    ap.add_argument("--work-root", default="example_outputs/emdbmatch_work")
    ap.add_argument("--json-root", default="example_outputs/json")
    ap.add_argument("--summary-out", default="example_outputs/qc/emdbmatch_nosource_modify_summary.csv")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    manifest = read_manifest(args.manifest)
    aux_index = load_aux_index(args.aux_index)
    web_roots = [x.strip() for x in args.web_roots.split(",") if x.strip()]
    entries = read_entries(args.entries_csv)

    if args.limit:
        entries = entries[:args.limit]

    work_root = Path(args.work_root)
    work_root.mkdir(parents=True, exist_ok=True)

    rows = []
    status_counter = Counter()

    for i, eid in enumerate(entries, 1):
        m = manifest.get(eid, {})
        pdb_id = m.get("pdb_id", "").strip().lower()

        row_base = {
            "entry_id": eid,
            "pdb_id": pdb_id,
        }

        if not pdb_id:
            row = {
                **row_base,
                "status": "failed_missing_pdb_id",
                "reason": "missing pdb_id",
            }
            rows.append(row)
            status_counter[row["status"]] += 1
            continue

        web_pdb = find_web_file(eid, web_roots, "pdb")
        if not web_pdb:
            row = {
                **row_base,
                "status": "failed_missing_web_pdb",
                "reason": "missing web pdb",
            }
            rows.append(row)
            status_counter[row["status"]] += 1
            continue

        groups = candidate_groups_for_pdb(pdb_id, aux_index)

        # Select groups with all three files required by emdbMatch.py:
        # .manAlign, .faa, and .cif.
        valid_groups = []
        for stem, g in groups.items():
            if "manAlign" in g and "faa" in g and "cif" in g:
                valid_groups.append((stem, g))

        if not valid_groups:
            row = {
                **row_base,
                "status": "failed_no_valid_batch_group",
                "reason": "no valid batch group with manAlign+faa+cif",
                "web_pdb": str(web_pdb),
                "candidate_group_count": len(groups),
                "valid_group_count": 0,
            }
            rows.append(row)
            status_counter[row["status"]] += 1
            continue

        # Default: use the first complete group.
        # If one PDB maps to multiple batch groups, a sequence-similarity based
        # ranking can be added later.
        stem, g = valid_groups[0]

        out_dir = work_root / eid
        out_dir.mkdir(parents=True, exist_ok=True)

        safe_stem = f"{eid}_{pdb_id}_{Path(stem).name}"

        local_pdb = copy_to_work(web_pdb, out_dir / f"{safe_stem}.pdb")
        local_manAlign = copy_to_work(g["manAlign"], out_dir / f"{safe_stem}.manAlign")
        local_faa = copy_to_work(g["faa"], out_dir / f"{safe_stem}.faa")
        local_cif = copy_to_work(g["cif"], out_dir / f"{safe_stem}.cif")

        out_autoDel = out_dir / f"{safe_stem}.autoDel"
        out_log = out_dir / f"{safe_stem}.emdbMatch.log"
        out_err = out_dir / f"{safe_stem}.emdbMatch.err"

        cmd = [
            "python",
            args.emdbmatch,
            "-i", str(local_pdb),
            "-s", str(local_manAlign),
            "-t", str(local_cif),
            "-o", str(out_autoDel),
        ]

        print(f"[{i}/{len(entries)}] {eid} pdb={pdb_id}")
        print("CMD:", " ".join(cmd))

        if args.dry_run:
            json_chains = get_json_chains(args.json_root, eid)
            row = {
                **row_base,
                "status": "dry_run",
                "web_pdb": str(web_pdb),
                "batch_stem": stem,
                "batch_manAlign": str(g["manAlign"]),
                "batch_faa": str(g["faa"]),
                "batch_cif": str(g["cif"]),
                "candidate_group_count": len(groups),
                "valid_group_count": len(valid_groups),
                "work_dir": str(out_dir),
                "cmd": " ".join(cmd),
                "json_chain_count": len(json_chains),
            }
            rows.append(row)
            status_counter[row["status"]] += 1
            continue

        with open(out_log, "w", encoding="utf-8") as fo, open(out_err, "w", encoding="utf-8") as fe:
            ret = subprocess.run(cmd, stdout=fo, stderr=fe)

        out_toRCSB = Path(str(local_manAlign).replace("manAlign", "toRCSB"))
        out_fasta = Path(str(local_manAlign).replace("manAlign", "fasta"))

        result_info = classify_emdbmatch_result(
            returncode=ret.returncode,
            out_toRCSB=out_toRCSB,
            out_fasta=out_fasta,
            out_autoDel=out_autoDel,
            json_root=args.json_root,
            entry_id=eid,
        )

        status_counter[result_info["status"]] += 1

        rows.append({
            **row_base,
            **result_info,
            "web_pdb": str(web_pdb),
            "batch_stem": stem,
            "batch_manAlign": str(g["manAlign"]),
            "batch_faa": str(g["faa"]),
            "batch_cif": str(g["cif"]),
            "batch_pdb": str(g["pdb"]) if "pdb" in g else "",
            "batch_toRCSB": str(g["toRCSB"]) if "toRCSB" in g else "",
            "batch_autoDel": str(g["autoDel"]) if "autoDel" in g else "",
            "candidate_group_count": len(groups),
            "valid_group_count": len(valid_groups),
            "work_dir": str(out_dir),
            "local_pdb": str(local_pdb),
            "local_manAlign": str(local_manAlign),
            "local_faa": str(local_faa),
            "local_cif": str(local_cif),
            "out_toRCSB": str(out_toRCSB),
            "out_fasta": str(out_fasta),
            "out_autoDel": str(out_autoDel),
            "stdout_log": str(out_log),
            "stderr_log": str(out_err),
            "cmd": " ".join(cmd),
        })

    out = Path(args.summary_out)
    out.parent.mkdir(parents=True, exist_ok=True)

    fields = []
    for r in rows:
        for k in r.keys():
            if k not in fields:
                fields.append(k)

    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    print("summary:", out)
    print("status:", dict(status_counter))


if __name__ == "__main__":
    main()