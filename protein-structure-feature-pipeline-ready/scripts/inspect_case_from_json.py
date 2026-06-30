#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
import subprocess
from pathlib import Path


INTEREST_SUFFIX = [
    ".cif",
    ".pdb",
    ".refine",
    ".cryofold",
    ".fa",
    ".faa",
    ".fas",
    ".fasta",
    ".fullLen",
    ".autoAlign",
    ".manAlign",
    ".autoDel",
    ".manDel",
    ".toRCSB",
    ".txt",
    ".log",
]


def read_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def run_file_cmd(path):
    try:
        r = subprocess.run(
            ["file", "-b", str(path)],
            text=True,
            capture_output=True,
            timeout=5,
        )
        return r.stdout.strip()
    except Exception as e:
        return repr(e)


def safe_head(path, n=100):
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            lines = []
            for i, line in enumerate(f):
                if i >= n:
                    break
                lines.append(line.rstrip("\n"))
        return "\n".join(lines)
    except Exception as e:
        return f"[HEAD ERROR] {repr(e)}"


def add_file(files, p):
    if not p:
        return
    try:
        p = Path(p)
        if p.exists():
            files.add(p)
        rp = p.resolve()
        if rp.exists():
            files.add(rp)
    except Exception:
        pass


def same_stem_candidates(path):
    out = []
    p = Path(path)
    parent = p.parent

    if not parent.exists():
        return out

    stem = p.name
    for suf in INTEREST_SUFFIX:
        if stem.endswith(suf):
            stem = stem[: -len(suf)]
            break

    for x in parent.iterdir():
        if not x.exists():
            continue
        if x.name.startswith(stem):
            out.append(x)

    return out


def find_by_patterns(root, patterns, max_hits=5000):
    root = Path(root)
    hits = []

    for pat in patterns:
        try:
            for p in root.rglob(f"*{pat}*"):
                if p.exists():
                    hits.append(p)
                if len(hits) >= max_hits:
                    return hits
        except Exception:
            pass

    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--entry-id", required=True)
    ap.add_argument("--json-root", default="example_outputs/json")
    ap.add_argument("--batch-root", default="example_data/raw/batch")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    entry_id = args.entry_id
    json_root = Path(args.json_root)
    batch_root = Path(args.batch_root)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    chain_dir = json_root / entry_id / "chains"
    chain_jsons = sorted(chain_dir.glob("*.json"))

    if not chain_jsons:
        raise SystemExit(f"No chain JSON found: {chain_dir}")

    first = read_json(chain_jsons[0])
    pdb_id = str(first.get("entry", {}).get("pdb_id", "")).strip()
    emdb_id = str(first.get("entry", {}).get("emdb_id", "")).strip()
    cfdb_id = str(first.get("entry", {}).get("cfdb_id", entry_id)).strip()

    files = set()

    # 1. Read source_files from all chain JSON files
    for jf in chain_jsons:
        d = read_json(jf)
        sf = d.get("source_files", {})
        for key in ["structure_path", "fasta_path"]:
            p = sf.get(key)
            if not p:
                continue
            add_file(files, Path(p))
            add_file(files, Path(p).resolve())

            # Files in the same directory with the same stem
            for cand in same_stem_candidates(Path(p)):
                add_file(files, cand)

            rp = Path(p).resolve()
            for cand in same_stem_candidates(rp):
                add_file(files, cand)

    # 2. Search batch filenames using entry_id / pdb_id / cfdb_id / emdb_id
    patterns = set()
    for x in [entry_id, pdb_id, pdb_id.lower(), pdb_id.upper(), cfdb_id, emdb_id]:
        if x and x.lower() not in {"na", "none", "null"}:
            patterns.add(x)

    for p in find_by_patterns(batch_root, patterns):
        add_file(files, p)

    files = sorted(files, key=lambda x: str(x))

    files_txt = out_dir / "files.txt"
    with open(files_txt, "w", encoding="utf-8") as f:
        for p in files:
            f.write(str(p) + "\n")

    print("entry:", entry_id)
    print("pdb_id:", pdb_id)
    print("emdb_id:", emdb_id)
    print("files found:", len(files))
    print("wrote:", files_txt)

    # 3. Output file heads
    for p in files:
        safe_name = str(p).replace("/", "_").replace(" ", "_")
        out = out_dir / f"{safe_name}.head.txt"

        with open(out, "w", encoding="utf-8") as f:
            f.write("=" * 100 + "\n")
            f.write(f"FILE: {p}\n")
            try:
                f.write(f"REALPATH: {p.resolve()}\n")
            except Exception:
                pass
            f.write(f"TYPE: {run_file_cmd(p)}\n")
            try:
                f.write(f"SIZE: {p.stat().st_size}\n")
            except Exception:
                pass
            f.write("=" * 100 + "\n")
            f.write(safe_head(p, n=120))
            f.write("\n")

    print("head files written to:", out_dir)


if __name__ == "__main__":
    main()