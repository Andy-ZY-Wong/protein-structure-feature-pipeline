#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import math
from pathlib import Path

manifest = Path("example_outputs/manifest/manifest_skeleton.csv")
out_dir = Path("example_outputs/manifests_3000")
chunk_size = 3000

out_dir.mkdir(parents=True, exist_ok=True)

rows = []
seen = set()

with open(manifest, newline="", encoding="utf-8", errors="ignore") as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames

    for r in reader:
        entry_id = r.get("entry_id", "").strip()
        if not entry_id or entry_id in seen:
            continue

        has_structure = r.get("has_structure", "").lower()
        skeleton_status = r.get("skeleton_status", "").lower()

        if has_structure and has_structure not in {"yes", "true", "1"}:
            continue
        if skeleton_status and skeleton_status != "done":
            continue

        rows.append(r)
        seen.add(entry_id)

n = len(rows)
num_chunks = math.ceil(n / chunk_size)

chunk_paths = []

for i in range(num_chunks):
    chunk = rows[i * chunk_size : (i + 1) * chunk_size]
    out = out_dir / f"manifest_chunk_{i:03d}.csv"
    chunk_paths.append(out)

    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(chunk)

list_path = out_dir / "chunks_list.txt"
with open(list_path, "w", encoding="utf-8") as f:
    for p in chunk_paths:
        f.write(str(p) + "\n")

print("total entries:", n)
print("chunk size:", chunk_size)
print("num chunks:", num_chunks)
print("chunk list:", list_path)

for p in chunk_paths:
    print(p)