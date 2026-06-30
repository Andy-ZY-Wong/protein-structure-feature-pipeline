#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import os
import re
import subprocess
from pathlib import Path
from collections import Counter, defaultdict


def safe_decode(b):
    for enc in ["utf-8", "latin-1", "gbk"]:
        try:
            return b.decode(enc, errors="replace")
        except Exception:
            pass
    return b.decode("utf-8", errors="replace")


def get_suffix(path: Path):
    if path.suffixes:
        return "".join(path.suffixes).lstrip(".")
    return "[no_suffix]"


def normalize_line(s):
    s = s.strip()
    s = re.sub(r"\s+", " ", s)
    s = re.sub(r"[-+]?\d+\.\d+", "<FLOAT>", s)
    s = re.sub(r"[-+]?\d+", "<NUM>", s)

    # 很长的纯序列行压缩成 <SEQ>
    if re.fullmatch(r"[A-Za-z\-*\.]{40,}", s):
        return "<SEQ>"

    if len(s) > 160:
        s = s[:160] + "..."

    return s


def first_nonempty_lines(text, n=8):
    out = []
    for line in text.splitlines():
        s = line.strip()
        if s:
            out.append(s)
        if len(out) >= n:
            break
    return out


def detect_kind(path: Path, text: str, raw: bytes):
    name = path.name.lower()
    suffix = get_suffix(path).lower()
    lines = first_nonempty_lines(text, n=20)
    first = lines[0] if lines else ""

    if raw.startswith(b"\x1f\x8b"):
        return "gzip_compressed"

    if not lines:
        return "empty_or_binary"

    first_upper = first.upper()

    # mmCIF
    if first.startswith("data_") or "_atom_site." in text[:5000] or "_entry.id" in text[:5000]:
        return "mmCIF_structure"

    # PDB
    pdb_prefixes = ("HEADER", "TITLE", "ATOM", "HETATM", "SEQRES", "MODEL", "CRYST1", "REMARK")
    if first_upper.startswith(pdb_prefixes) or "\nATOM  " in text[:5000] or "\nHETATM" in text[:5000]:
        return "PDB_structure"

    # FASTA
    if first.startswith(">"):
        return "FASTA_sequence"

    # JSON
    if first.startswith("{") or first.startswith("["):
        return "JSON"

    # CSV / TSV
    if "\t" in first:
        return "TSV_or_tabular_text"
    if "," in first and len(first.split(",")) >= 3:
        return "CSV_or_comma_table"

    # 常见日志
    log_words = ["error", "warning", "traceback", "finished", "started", "processing", "job", "slurm"]
    low = text[:2000].lower()
    if any(w in low for w in log_words) and suffix in {"log", "txt", "out", "err"}:
        return "log_or_runtime_text"

    # 对齐/序列处理类文件特征
    if any(x in low for x in ["align", "alignment", "delete", "deleted", "full", "chain", "residue"]):
        return "custom_sequence_or_alignment_text"

    return "plain_text_unknown"


def file_command(path: Path):
    try:
        r = subprocess.run(
            ["file", "-b", str(path)],
            text=True,
            capture_output=True,
            timeout=5
        )
        if r.returncode == 0:
            return r.stdout.strip()
    except Exception:
        pass
    return ""


def make_signature(text):
    lines = first_nonempty_lines(text, n=6)
    if not lines:
        return "[empty]"

    norm = [normalize_line(x) for x in lines]
    sig = " || ".join(norm)
    if len(sig) > 500:
        sig = sig[:500] + "..."
    return sig


def short_head(text, max_lines=40):
    lines = text.splitlines()[:max_lines]
    return "\n".join(lines)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", default="example_data/raw/batch")
    p.add_argument("--out-dir", default="example_outputs/qc/batch_content_inspect")
    p.add_argument("--read-bytes", type=int, default=20000)
    p.add_argument("--sample-lines", type=int, default=60)
    p.add_argument("--max-samples-per-suffix", type=int, default=12)
    args = p.parse_args()

    root = Path(args.root)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    file_audit_csv = out_dir / "batch_file_content_audit.csv"
    summary_csv = out_dir / "batch_suffix_content_summary.csv"
    samples_md = out_dir / "batch_suffix_samples.md"

    suffix_counter = Counter()
    suffix_kind_counter = defaultdict(Counter)
    suffix_sig_counter = defaultdict(Counter)
    suffix_filetype_counter = defaultdict(Counter)

    examples = defaultdict(dict)

    rows = []

    files = [p for p in root.rglob("*") if p.is_file()]
    total = len(files)

    print(f"[INFO] root: {root}")
    print(f"[INFO] files found: {total}")
    print(f"[INFO] out_dir: {out_dir}")

    for i, path in enumerate(files, start=1):
        if i % 5000 == 0:
            print(f"[INFO] processed {i}/{total}")

        rel = str(path.relative_to(root))
        suffix = get_suffix(path)

        try:
            size = path.stat().st_size
        except Exception:
            size = 0

        try:
            with open(path, "rb") as f:
                raw = f.read(args.read_bytes)
        except Exception as e:
            row = {
                "path": str(path),
                "relative_path": rel,
                "suffix": suffix,
                "size_bytes": size,
                "detected_kind": "read_error",
                "first_line": "",
                "signature": "",
                "file_command": "",
                "error": repr(e),
            }
            rows.append(row)
            continue

        text = safe_decode(raw)
        lines = first_nonempty_lines(text, n=8)
        first_line = lines[0] if lines else ""

        detected = detect_kind(path, text, raw)
        sig = make_signature(text)
        fcmd = file_command(path)

        suffix_counter[suffix] += 1
        suffix_kind_counter[suffix][detected] += 1
        suffix_sig_counter[suffix][sig] += 1
        suffix_filetype_counter[suffix][fcmd] += 1

        # 每个 suffix 下保留若干种代表性签名的样本
        if len(examples[suffix]) < args.max_samples_per_suffix:
            if sig not in examples[suffix]:
                examples[suffix][sig] = {
                    "path": str(path),
                    "relative_path": rel,
                    "detected_kind": detected,
                    "file_command": fcmd,
                    "head": short_head(text, args.sample_lines),
                }

        rows.append({
            "path": str(path),
            "relative_path": rel,
            "suffix": suffix,
            "size_bytes": size,
            "detected_kind": detected,
            "first_line": first_line[:300],
            "signature": sig,
            "file_command": fcmd,
            "error": "",
        })

    # 每个文件一行的详细表
    with open(file_audit_csv, "w", newline="", encoding="utf-8") as f:
        fields = [
            "path",
            "relative_path",
            "suffix",
            "size_bytes",
            "detected_kind",
            "first_line",
            "signature",
            "file_command",
            "error",
        ]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    # 每个 suffix 的汇总表
    summary_rows = []
    for suffix, count in suffix_counter.most_common():
        kind_summary = "; ".join(
            f"{k}:{v}" for k, v in suffix_kind_counter[suffix].most_common()
        )

        filetype_summary = "; ".join(
            f"{k}:{v}" for k, v in suffix_filetype_counter[suffix].most_common(5)
        )

        top_sigs = []
        for sig, n in suffix_sig_counter[suffix].most_common(5):
            top_sigs.append(f"[{n}] {sig}")

        summary_rows.append({
            "suffix": suffix,
            "count": count,
            "detected_kind_summary": kind_summary,
            "file_command_top": filetype_summary,
            "top_signatures": " || ".join(top_sigs),
        })

    with open(summary_csv, "w", newline="", encoding="utf-8") as f:
        fields = [
            "suffix",
            "count",
            "detected_kind_summary",
            "file_command_top",
            "top_signatures",
        ]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(summary_rows)

    # 每种 suffix 的代表内容样本
    with open(samples_md, "w", encoding="utf-8") as f:
        f.write("# Batch file content samples\n\n")
        f.write(f"Root: `{root}`\n\n")
        f.write(f"Total files scanned: `{total}`\n\n")
        f.write(f"Read bytes per file: `{args.read_bytes}`\n\n")

        for suffix, count in suffix_counter.most_common():
            f.write(f"\n\n## EXT: {suffix}  count={count}\n\n")
            f.write("Detected kind summary:\n\n")
            for k, v in suffix_kind_counter[suffix].most_common():
                f.write(f"- {k}: {v}\n")

            f.write("\nRepresentative samples:\n\n")

            for j, (sig, ex) in enumerate(examples[suffix].items(), start=1):
                f.write(f"\n### Sample {j}\n\n")
                f.write(f"- path: `{ex['path']}`\n")
                f.write(f"- relative_path: `{ex['relative_path']}`\n")
                f.write(f"- detected_kind: `{ex['detected_kind']}`\n")
                f.write(f"- file_command: `{ex['file_command']}`\n")
                f.write(f"- signature: `{sig}`\n\n")
                f.write("```text\n")
                f.write(ex["head"])
                f.write("\n```\n")

    print("[DONE]")
    print("Detailed file audit:", file_audit_csv)
    print("Suffix summary:", summary_csv)
    print("Content samples:", samples_md)


if __name__ == "__main__":
    main()