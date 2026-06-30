#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import re
from pathlib import Path
from collections import defaultdict, Counter


THREE2ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    "MSE": "M", "SEC": "U", "PYL": "O", "UNK": "X",
}


def norm_seq(s):
    return re.sub(r"[^A-Za-z]", "", str(s or "")).upper().replace("*", "")


def read_manifest(manifest):
    out = {}
    p = Path(manifest)
    if not p.exists():
        return out

    with open(p, newline="", encoding="utf-8", errors="ignore") as f:
        for r in csv.DictReader(f):
            eid = r.get("entry_id", "").strip()
            if not eid:
                continue
            out[eid] = r

    return out


def read_pdb_ca_seq(path):
    chains = defaultdict(list)
    seen = set()

    try:
        with open(path, encoding="utf-8", errors="ignore") as f:
            for l in f:
                if not l.startswith("ATOM"):
                    continue
                if l[12:16].strip() != "CA":
                    continue

                # Handle ordinary single-character chains and possible two-character chains
                chain = l[20:22].strip() or l[21].strip()
                resname = l[17:20].strip()
                resseq = l[22:27].strip()

                key = (chain, resseq)
                if key in seen:
                    continue
                seen.add(key)

                chains[chain].append(THREE2ONE.get(resname, "X"))
    except Exception:
        return {}

    return {c: "".join(s) for c, s in chains.items()}


def read_fasta_like(path):
    records = []
    header = None
    seq = []

    try:
        with open(path, encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    if header is not None:
                        records.append((header, norm_seq("".join(seq))))
                    header = line[1:].strip()
                    seq = []
                else:
                    seq.append(line)

        if header is not None:
            records.append((header, norm_seq("".join(seq))))
    except Exception:
        return []

    return [(h, s) for h, s in records if s]


def seq_similarity(a, b):
    """
    Return a simple but practical similarity summary:
    - exact: whether sequences are exactly identical
    - containment: whether one sequence contains the other
    - identity_shorter: per-position identity over the shorter length
    - coverage_shorter: shorter length / longer length
    """
    a = norm_seq(a)
    b = norm_seq(b)

    if not a or not b:
        return {
            "exact": False,
            "containment": False,
            "identity_shorter": 0.0,
            "coverage_shorter": 0.0,
            "len_a": len(a),
            "len_b": len(b),
        }

    exact = a == b
    containment = (a in b) or (b in a)

    n = min(len(a), len(b))
    m = sum(1 for x, y in zip(a[:n], b[:n]) if x == y)

    return {
        "exact": exact,
        "containment": containment,
        "identity_shorter": m / n if n else 0.0,
        "coverage_shorter": n / max(len(a), len(b)),
        "len_a": len(a),
        "len_b": len(b),
    }


def best_match(query_records, target_records):
    """
    query_records: [(name, seq)]
    target_records: [(name, seq)]

    Return the best target match for each query and summarize the overall result.
    """
    rows = []
    exact_count = 0
    containment_count = 0

    for qh, qs in query_records:
        best = None

        for th, ts in target_records:
            sim = seq_similarity(qs, ts)

            score = (
                1 if sim["exact"] else 0,
                1 if sim["containment"] else 0,
                sim["identity_shorter"],
                sim["coverage_shorter"],
            )

            item = {
                "query_header": qh,
                "target_header": th,
                **sim,
                "score": score,
            }

            if best is None or item["score"] > best["score"]:
                best = item

        if best:
            rows.append(best)
            if best["exact"]:
                exact_count += 1
            if best["containment"]:
                containment_count += 1

    avg_identity = sum(r["identity_shorter"] for r in rows) / len(rows) if rows else 0.0
    avg_coverage = sum(r["coverage_shorter"] for r in rows) / len(rows) if rows else 0.0

    return {
        "query_count": len(query_records),
        "target_count": len(target_records),
        "matched_count": len(rows),
        "exact_count": exact_count,
        "containment_count": containment_count,
        "avg_identity_shorter": avg_identity,
        "avg_coverage_shorter": avg_coverage,
        "best_rows": rows,
    }


def load_aux_index(path):
    idx = defaultdict(list)
    p = Path(path)
    if not p.exists():
        return idx

    with open(p, newline="", encoding="utf-8", errors="ignore") as f:
        for r in csv.DictReader(f):
            pdb = r.get("pdb_id", "").strip().lower()
            if pdb:
                idx[pdb].append(r)

    return idx


def find_web_files(entry_id, web_roots):
    files = []

    for root in web_roots:
        root = Path(root)
        if not root.exists():
            continue

        for suf in ["pdb", "fasta", "fa"]:
            p = root / f"{entry_id}.{suf}"
            if p.exists():
                files.append((suf, p))

    return files


def candidate_batch_files(entry_id, pdb_id, aux_index, batch_root, max_name_hits=100):
    files = {}

    pdb_id = str(pdb_id or "").lower()

    def add_file(path, suffix=None):
        path = Path(path)
        if path.exists() and path.is_file():
            suf = suffix or path.name.split(".")[-1]
            files[str(path)] = (suf, path)

    def add_same_stem_siblings(path):
        """
        If 250523/9397.cif is found by pdb_id, also add 250523/9397.faa,
        250523/9397.manAlign, 250523/9397.fullLen, etc.
        """
        path = Path(path)
        parent = path.parent
        stem = path.stem

        wanted_suffixes = [
            "pdb", "cif", "fasta", "fa", "faa",
            "fullLen", "manAlign", "autoAlign",
            "toRCSB", "autoDel", "manDel",
        ]

        for suf in wanted_suffixes:
            q = parent / f"{stem}.{suf}"
            if q.exists():
                add_file(q, suf)

        # Handle files without dot suffixes or with special naming
        try:
            for q in parent.iterdir():
                if q.is_file() and q.name.startswith(stem):
                    add_file(q)
        except Exception:
            pass

    # 1. Search batch_aux_file_index.csv by pdb_id
    for r in aux_index.get(pdb_id, []):
        p = Path(r.get("path", ""))
        if p.exists():
            add_file(p, r.get("suffix", ""))
            add_same_stem_siblings(p)

    # 2. Also add files whose names match entry_id / pdb_id
    root = Path(batch_root)
    if root.exists():
        patterns = []
        if entry_id:
            patterns.append(f"*{entry_id}*")
        if pdb_id:
            patterns.append(f"*{pdb_id}*")
            patterns.append(f"*{pdb_id.upper()}*")

        hits = []
        for pat in patterns:
            try:
                hits.extend(root.rglob(pat))
            except Exception:
                pass
            if len(hits) > max_name_hits:
                break

        for p in hits[:max_name_hits]:
            if p.is_file():
                add_file(p)
                add_same_stem_siblings(p)

    return list(files.values())


def file_to_records(suffix, path):
    suffix = str(suffix or "").lower()

    if suffix == "pdb" or str(path).lower().endswith(".pdb"):
        d = read_pdb_ca_seq(path)
        return [(f"{Path(path).name}|chain:{c}", s) for c, s in d.items()]

    # Treat these files as FASTA-like inputs
    if suffix in {
        "fasta", "fa", "faa", "fulllen", "manalign", "autoalign",
        "torcsb", "autodel", "mandel", "fas"
    }:
        return [(f"{Path(path).name}|{h}", s) for h, s in read_fasta_like(path)]

    return []


def web_to_records(web_files):
    out = {}

    for suf, p in web_files:
        if suf == "pdb":
            d = read_pdb_ca_seq(p)
            out["web_pdb"] = {
                "path": str(p),
                "records": [(f"{p.name}|chain:{c}", s) for c, s in d.items()],
            }
        elif suf in {"fasta", "fa"}:
            out["web_fasta"] = {
                "path": str(p),
                "records": [(f"{p.name}|{h}", s) for h, s in read_fasta_like(p)],
            }

    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--entry-id", required=True)
    ap.add_argument("--manifest", default="example_outputs/manifest/manifest_skeleton.csv")
    ap.add_argument("--web-roots", default="example_data/raw/web")
    ap.add_argument("--batch-root", default="example_data/raw/batch")
    ap.add_argument("--aux-index", default="example_outputs/qc/batch_aux_file_index.csv")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    manifest = read_manifest(args.manifest)
    row = manifest.get(args.entry_id, {})
    pdb_id = row.get("pdb_id", "").strip().lower()

    web_roots = [x.strip() for x in args.web_roots.split(",") if x.strip()]
    web_files = find_web_files(args.entry_id, web_roots)
    web_records = web_to_records(web_files)

    aux_index = load_aux_index(args.aux_index)
    batch_files = candidate_batch_files(args.entry_id, pdb_id, aux_index, args.batch_root)

    rows = []

    for web_type, web_obj in web_records.items():
        q_records = web_obj["records"]

        for suffix, bp in batch_files:
            t_records = file_to_records(suffix, bp)
            if not t_records:
                continue

            bm = best_match(q_records, t_records)

            rows.append({
                "entry_id": args.entry_id,
                "pdb_id": pdb_id,
                "web_type": web_type,
                "web_path": web_obj["path"],
                "batch_suffix": suffix,
                "batch_path": str(bp),
                "web_record_count": bm["query_count"],
                "batch_record_count": bm["target_count"],
                "matched_count": bm["matched_count"],
                "exact_count": bm["exact_count"],
                "containment_count": bm["containment_count"],
                "avg_identity_shorter": round(bm["avg_identity_shorter"], 4),
                "avg_coverage_shorter": round(bm["avg_coverage_shorter"], 4),
            })

    rows.sort(
        key=lambda r: (
            r["web_type"],
            int(r["exact_count"]),
            int(r["containment_count"]),
            float(r["avg_identity_shorter"]),
            float(r["avg_coverage_shorter"]),
        ),
        reverse=True,
    )

    out = Path(args.out) if args.out else Path(f"example_outputs/qc/web_batch_similarity_{args.entry_id}.csv")
    out.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "entry_id",
        "pdb_id",
        "web_type",
        "web_path",
        "batch_suffix",
        "batch_path",
        "web_record_count",
        "batch_record_count",
        "matched_count",
        "exact_count",
        "containment_count",
        "avg_identity_shorter",
        "avg_coverage_shorter",
    ]

    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    print("entry_id:", args.entry_id)
    print("pdb_id:", pdb_id)
    print("web files:")
    for suf, p in web_files:
        print(" ", suf, p)

    print("batch candidate files:", len(batch_files))
    print("summary:", out)

    print("\nTop matches:")
    for r in rows[:30]:
        print(
            r["web_type"],
            r["batch_suffix"],
            "exact=", r["exact_count"],
            "contain=", r["containment_count"],
            "id=", r["avg_identity_shorter"],
            "cov=", r["avg_coverage_shorter"],
            r["batch_path"],
        )


if __name__ == "__main__":
    main()