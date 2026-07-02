#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
00b_add_fasta_to_manifest.py

Purpose:
Append FASTA file path information to an existing manifest.csv.

Input:
- Existing manifest.csv
- One or more directories, for example:
  example_data/raw/web/260304
  example_data/raw/web/260305

Output:
- New manifest_with_fasta.csv
- Added fields:
  fasta_path
  has_fasta
  fasta_matched_by

Matching priority:
1. Exact filename match to cfdb_id / entry_id
2. Exact filename match to pdb_id
3. Path contains cfdb_id / entry_id
4. Path contains pdb_id
"""

import argparse
import csv
import os
from pathlib import Path
from typing import Dict, List, Tuple


FASTA_EXTS = (
    ".fasta",
    ".fa",
    ".faa",
    ".fas",
    ".seq",
    ".fasta.gz",
    ".fa.gz",
    ".faa.gz",
    ".fas.gz",
)


def normalize_id(x: str) -> str:
    if x is None:
        return ""
    return str(x).strip().lower()


def strip_fasta_suffix(filename: str) -> str:
    name = filename
    for ext in [
        ".fasta.gz", ".fa.gz", ".faa.gz", ".fas.gz",
        ".fasta", ".fa", ".faa", ".fas", ".seq"
    ]:
        if name.lower().endswith(ext):
            return name[: -len(ext)]
    return Path(name).stem


def is_fasta_file(path: Path) -> bool:
    lower = path.name.lower()
    return any(lower.endswith(ext) for ext in FASTA_EXTS)


def index_fasta_files(fasta_dirs: List[str]) -> Tuple[Dict[str, List[str]], List[str]]:
    """
    Scan directories and build a FASTA index:
    stem.lower() -> [full_path]
    """

    index: Dict[str, List[str]] = {}
    all_files: List[str] = []

    for d in fasta_dirs:
        root = Path(d)

        if not root.exists():
            print(f"[WARN] FASTA directory not found: {d}")
            continue

        print(f"[INFO] Scanning FASTA directory: {d}")

        for path in root.rglob("*"):
            if not path.is_file():
                continue

            if not is_fasta_file(path):
                continue

            full_path = str(path)
            all_files.append(full_path)

            stem = normalize_id(strip_fasta_suffix(path.name))
            if stem:
                index.setdefault(stem, []).append(full_path)

    print(f"[INFO] Total FASTA files indexed: {len(all_files)}")
    print(f"[INFO] Unique FASTA stems indexed: {len(index)}")

    return index, all_files


def find_fasta_path(
    entry_id: str,
    cfdb_id: str,
    pdb_id: str,
    fasta_index: Dict[str, List[str]],
    all_fasta_files: List[str],
) -> Tuple[str, str]:
    """
    Find FASTA files according to priority.
    Returns:
    fasta_path, matched_by
    """

    entry_key = normalize_id(entry_id)
    cfdb_key = normalize_id(cfdb_id)
    pdb_key = normalize_id(pdb_id)

    # 1. Exact filename match to CFDB id
    if cfdb_key and cfdb_key in fasta_index:
        return fasta_index[cfdb_key][0], "exact_cfdb"

    # 2. Exact filename match to entry_id
    if entry_key and entry_key in fasta_index:
        return fasta_index[entry_key][0], "exact_entry"

    # 3. Exact filename match to PDB id
    if pdb_key and pdb_key in fasta_index:
        return fasta_index[pdb_key][0], "exact_pdb"

    # 4. Path contains CFDB id
    if cfdb_key:
        for p in all_fasta_files:
            if cfdb_key in p.lower():
                return p, "contains_cfdb"

    # 5. Path contains entry_id
    if entry_key:
        for p in all_fasta_files:
            if entry_key in p.lower():
                return p, "contains_entry"

    # 6. Path contains PDB id
    if pdb_key:
        for p in all_fasta_files:
            if pdb_key in p.lower():
                return p, "contains_pdb"

    return "", "not_found"


def add_fasta_to_manifest(manifest_in: str, fasta_dirs: List[str], manifest_out: str) -> None:
    if not os.path.exists(manifest_in):
        raise FileNotFoundError(f"Cannot find manifest: {manifest_in}")

    os.makedirs(os.path.dirname(manifest_out) or ".", exist_ok=True)

    fasta_index, all_fasta_files = index_fasta_files(fasta_dirs)

    rows_out = []

    with open(manifest_in, "r", encoding="utf-8", errors="ignore", newline="") as fin:
        reader = csv.DictReader(fin)

        if reader.fieldnames is None:
            raise RuntimeError(f"No header found in manifest: {manifest_in}")

        fieldnames = reader.fieldnames

        total = 0
        found = 0
        missing = 0

        for row in reader:
            total += 1

            entry_id = row.get("entry_id", "")
            cfdb_id = row.get("cfdb_id", "")
            pdb_id = row.get("pdb_id", "")

            fasta_path, matched_by = find_fasta_path(
                entry_id=entry_id,
                cfdb_id=cfdb_id,
                pdb_id=pdb_id,
                fasta_index=fasta_index,
                all_fasta_files=all_fasta_files,
            )

            if fasta_path:
                row["fasta_path"] = fasta_path
                row["has_fasta"] = "yes"
                row["fasta_matched_by"] = matched_by
                found += 1
            else:
                row["fasta_path"] = ""
                row["has_fasta"] = "no"
                row["fasta_matched_by"] = "not_found"
                missing += 1

                old_note = row.get("note", "")
                if old_note:
                    row["note"] = old_note + "; No matching FASTA file found"
                else:
                    row["note"] = "No matching FASTA file found"

            rows_out.append(row)

    # Output fields: keep original fields and append FASTA fields
    out_fields = list(fieldnames)

    for col in ["fasta_path", "has_fasta", "fasta_matched_by"]:
        if col not in out_fields:
            out_fields.append(col)

    with open(manifest_out, "w", encoding="utf-8", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(rows_out)

    print("========== FASTA Manifest Summary ==========")
    print(f"Input manifest rows    : {total}")
    print(f"Rows with FASTA found  : {found}")
    print(f"Rows missing FASTA     : {missing}")
    print(f"Output manifest        : {manifest_out}")
    print("===========================================")


def main():
    parser = argparse.ArgumentParser(
        description="Add FASTA paths to an existing manifest.csv."
    )

    parser.add_argument(
        "--manifest",
        required=True,
        help="Path to the existing manifest.csv"
    )

    parser.add_argument(
        "--fasta-dirs",
        nargs="+",
        required=True,
        help="One or more directories containing FASTA files"
    )

    parser.add_argument(
        "--out",
        required=True,
        help="Output path for the new manifest_with_fasta.csv"
    )

    args = parser.parse_args()

    add_fasta_to_manifest(
        manifest_in=args.manifest,
        fasta_dirs=args.fasta_dirs,
        manifest_out=args.out,
    )


if __name__ == "__main__":
    main()