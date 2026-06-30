#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
04_reparse_qscore_raw.py

功能：
从已经保存好的 ChimeraX QScore raw 文件重新解析 Q-score，
并回填到已有 chain JSON 中。

不会重新运行 ChimeraX。
不会重新打开 map。
不会重新计算 Q-score。

默认 raw 文件命名：
  example_outputs/qc/qscore_raw_all/cf1a0a.qscore.raw.txt

默认 JSON 目录：
  example_outputs/json/cf1a0a/chains/A.json
"""

import argparse
import json
import os
import re
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple
from collections import Counter


def now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def read_json(path: Path) -> Dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def write_json(path: Path, data: Dict) -> None:
    tmp = str(path) + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    os.replace(tmp, path)


def append_log(data: Dict, step: str, status: str, message: str) -> None:
    data.setdefault("logs", [])
    data["logs"].append({
        "time": now_iso(),
        "step": step,
        "status": status,
        "message": message
    })


def normalize_chain_id(x: str) -> str:
    x = str(x or "").strip()
    if "/" in x:
        x = x.split("/")[-1]
    if "#" in x:
        x = x.split("#")[-1]
    return x.strip()


def parse_residue_number(text: str) -> Optional[Tuple[int, str]]:
    """
    支持：
      50
      50A
      A:50
      /A:50
    """
    text = str(text or "").strip()
    if not text:
        return None

    m = re.search(r"(-?\d+)([A-Za-z]?)", text)
    if not m:
        return None

    return int(m.group(1)), m.group(2).strip()


def parse_qscore_raw(raw_path: Path) -> Dict[str, Dict[int, float]]:
    """
    解析 ChimeraX QScore raw 文件。

    识别格式：

      ChainId A - Protein
          Name    Number    Q_backbone    Q_side_chain    Q_residue    Q_expected@3.00A
          ALA     50        0.6998        -0.0971         0.5006       0.5804

    返回：
      {
        "A": {
          50: 0.5006,
          51: 0.4972,
          ...
        },
        "B": {...}
      }

    注意：
    - 使用 Q_residue 作为 residue-level Q-score；
    - 没有 residue name / 没有 Q_residue 的空行会跳过；
    - Q_side_chain 为 N/A 不影响，因为我们不用它。
    """

    result: Dict[str, Dict[int, float]] = {}

    current_chain = None
    in_protein_table = False
    header = None

    with open(raw_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue

            # 进入某条链的 residue-level 表格
            m = re.match(r"^ChainId\s+(\S+)\s+-\s+Protein", s)
            if m:
                current_chain = normalize_chain_id(m.group(1))
                result.setdefault(current_chain, {})
                in_protein_table = True
                header = None
                continue

            # 遇到 ChainId A 这种 summary 行，不是 residue 表格
            if s.startswith("ChainId ") and "- Protein" not in s:
                current_chain = None
                in_protein_table = False
                header = None
                continue

            if not in_protein_table or not current_chain:
                continue

            parts = s.split()
            low = [p.lower() for p in parts]

            # 表头：Name Number Q_backbone Q_side_chain Q_residue Q_expected@...
            if "name" in low and "number" in low and "q_residue" in low:
                header = low
                continue

            if header is None:
                continue

            try:
                name_i = header.index("name")
                number_i = header.index("number")
                q_i = header.index("q_residue")
            except ValueError:
                continue

            if len(parts) <= max(name_i, number_i, q_i):
                continue

            resname = parts[name_i].strip()
            resnum_text = parts[number_i].strip()
            q_text = parts[q_i].strip()

            # 跳过空残基行，比如：
            # 118    0.5804
            if not resname or len(resname) < 3:
                continue

            if q_text.upper() == "N/A":
                continue

            parsed = parse_residue_number(resnum_text)
            if not parsed:
                continue

            try:
                resseq, _icode = parsed
                q = float(q_text)
            except Exception:
                continue

            result[current_chain][resseq] = round(q, 4)

    return result


def build_residue_mapping(chain_json: Dict) -> Dict[int, int]:
    """
    建立：
      pdb_resseq -> chain_index

    如果有重复编号，保留第一次出现。
    """

    mapping: Dict[int, int] = {}

    for item in chain_json.get("residue_mapping", []):
        try:
            chain_index = int(item.get("chain_index"))

            if "pdb_resseq" in item:
                resseq = int(item.get("pdb_resseq"))
            elif "resseq" in item:
                resseq = int(item.get("resseq"))
            elif "pdb_residue_number" in item:
                resseq = int(item.get("pdb_residue_number"))
            else:
                continue

            mapping.setdefault(resseq, chain_index)

        except Exception:
            continue

    return mapping


def update_chain_qscore(
    chain_json_path: Path,
    raw_scores_for_chain: Dict[int, float],
    min_coverage: float = 0.2
) -> Tuple[str, int, int, float]:
    """
    把某条链的 raw qscore 写回 chain JSON。

    返回：
      status, mapped_count, length, coverage
    """

    d = read_json(chain_json_path)

    chain_id = str(d.get("chain", {}).get("chain_id", chain_json_path.stem))
    length = int(d.get("chain", {}).get("length", 0))

    if length <= 0:
        status = "failed"
        reason = "Invalid chain length."

        d.setdefault("computed_features", {})
        d["computed_features"]["qscore"] = {
            "status": status,
            "type": "values",
            "source": {"method": "ChimeraX QScore"},
            "coordinate_system": "chain_index_1_based",
            "reason": reason,
            "values": [],
            "mapped_count": 0,
            "coverage": 0.0,
            "mean": None
        }

        write_json(chain_json_path, d)
        return status, 0, length, 0.0

    resseq_to_idx = build_residue_mapping(d)

    values: List[Optional[float]] = [None] * length

    for resseq, q in raw_scores_for_chain.items():
        idx = resseq_to_idx.get(resseq)
        if idx is None:
            continue

        if 1 <= idx <= length:
            values[idx - 1] = q

    valid = [v for v in values if v is not None]
    mapped_count = len(valid)
    coverage = mapped_count / length if length > 0 else 0.0

    if mapped_count == 0:
        status = "failed"
        reason = "Raw Q-score file parsed, but no residue could be mapped to chain_index."
    elif coverage < min_coverage:
        status = "partial"
        reason = f"Low Q-score coverage: {coverage:.4f}"
    else:
        status = "done"
        reason = ""

    qscore_obj = {
        "status": status,
        "type": "values",
        "source": {
            "method": "ChimeraX QScore",
            "value_column": "Q_residue",
            "note": "Re-parsed from saved ChimeraX QScore raw output."
        },
        "coordinate_system": "chain_index_1_based",
        "values": values,
        "mapped_count": mapped_count,
        "coverage": round(coverage, 4),
        "mean": round(sum(valid) / len(valid), 4) if valid else None
    }

    if reason:
        qscore_obj["reason"] = reason

    d.setdefault("computed_features", {})
    d["computed_features"]["qscore"] = qscore_obj

    d.setdefault("status", {})
    d["status"]["qscore"] = status

    append_log(
        d,
        step="qscore_reparse",
        status=status,
        message=(
            f"Re-parsed Q-score from raw file. "
            f"chain={chain_id}; mapped_count={mapped_count}; "
            f"length={length}; coverage={coverage:.4f}"
        )
    )

    write_json(chain_json_path, d)

    return status, mapped_count, length, coverage


def update_entry_status(json_root: Path, entry_id: str) -> None:
    """
    根据所有 chain qscore 状态更新 entry.json。
    """

    entry_json_path = json_root / entry_id / "entry.json"
    chain_dir = json_root / entry_id / "chains"

    if not entry_json_path.exists() or not chain_dir.is_dir():
        return

    statuses = []

    for jf in sorted(chain_dir.glob("*.json")):
        try:
            d = read_json(jf)
            st = d.get("computed_features", {}).get("qscore", {}).get("status", "missing")
            statuses.append(st)
        except Exception:
            statuses.append("read_error")

    if not statuses:
        final = "failed"
    elif all(s == "done" for s in statuses):
        final = "done"
    elif any(s in {"done", "partial"} for s in statuses):
        final = "partial"
    else:
        final = "failed"

    try:
        entry = read_json(entry_json_path)
        entry.setdefault("status", {})
        entry["status"]["qscore"] = final

        append_log(
            entry,
            step="qscore_reparse",
            status=final,
            message=f"Entry-level Q-score status updated from chain statuses: {dict(Counter(statuses))}"
        )

        write_json(entry_json_path, entry)

    except Exception:
        return


def entry_id_from_raw(raw_path: Path) -> str:
    """
    cf1a0a.qscore.raw.txt -> cf1a0a
    """
    name = raw_path.name
    if name.endswith(".qscore.raw.txt"):
        return name.replace(".qscore.raw.txt", "")
    if name.endswith(".raw.txt"):
        return name.replace(".raw.txt", "")
    return raw_path.stem


def process_one_raw(
    raw_path: Path,
    json_root: Path,
    min_coverage: float = 0.2,
    verbose: bool = False
) -> Dict:
    entry_id = entry_id_from_raw(raw_path)
    entry_dir = json_root / entry_id
    chain_dir = entry_dir / "chains"

    summary = {
        "entry_id": entry_id,
        "raw_path": str(raw_path),
        "status": "failed",
        "chain_count": 0,
        "done": 0,
        "partial": 0,
        "failed": 0,
        "missing_chain_json": 0,
        "note": ""
    }

    if not chain_dir.is_dir():
        summary["note"] = "chain JSON directory not found"
        return summary

    try:
        parsed = parse_qscore_raw(raw_path)
    except Exception as e:
        summary["note"] = f"failed to parse raw file: {repr(e)}"
        return summary

    if not parsed:
        summary["note"] = "no Q-score residue rows parsed from raw file"
        return summary

    chain_jsons = {
        p.stem: p
        for p in sorted(chain_dir.glob("*.json"))
    }

    for chain_id, chain_path in chain_jsons.items():
        summary["chain_count"] += 1

        raw_scores = parsed.get(chain_id)

        if raw_scores is None:
            summary["missing_chain_json"] += 1

            # 这个链没有 raw qscore，写成 failed
            try:
                d = read_json(chain_path)
                d.setdefault("computed_features", {})
                d["computed_features"]["qscore"] = {
                    "status": "failed",
                    "type": "values",
                    "source": {
                        "method": "ChimeraX QScore",
                        "value_column": "Q_residue",
                        "note": "Re-parsed from saved ChimeraX QScore raw output."
                    },
                    "coordinate_system": "chain_index_1_based",
                    "reason": "This chain was not found in Q-score raw output.",
                    "values": [],
                    "mapped_count": 0,
                    "coverage": 0.0,
                    "mean": None
                }
                d.setdefault("status", {})
                d["status"]["qscore"] = "failed"
                write_json(chain_path, d)
            except Exception:
                pass

            summary["failed"] += 1
            continue

        st, mapped_count, length, coverage = update_chain_qscore(
            chain_json_path=chain_path,
            raw_scores_for_chain=raw_scores,
            min_coverage=min_coverage
        )

        summary[st] = summary.get(st, 0) + 1

        if verbose:
            print(
                f"[{entry_id}] chain={chain_id} "
                f"status={st} mapped={mapped_count}/{length} coverage={coverage:.4f}"
            )

    update_entry_status(json_root, entry_id)

    if summary["failed"] == 0 and summary["partial"] == 0 and summary["done"] > 0:
        summary["status"] = "done"
    elif summary["done"] > 0 or summary["partial"] > 0:
        summary["status"] = "partial"
    else:
        summary["status"] = "failed"

    return summary


def write_summary_csv(path: Path, rows: List[Dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "entry_id",
        "raw_path",
        "status",
        "chain_count",
        "done",
        "partial",
        "failed",
        "missing_chain_json",
        "note"
    ]

    import csv

    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main():
    p = argparse.ArgumentParser(
        description="Re-parse saved ChimeraX QScore raw outputs and update chain JSON files."
    )

    p.add_argument(
        "--raw-dir",
        default="example_outputs/qc/qscore_raw_all",
        help="保存 *.qscore.raw.txt 的目录"
    )
    p.add_argument(
        "--json-root",
        default="example_outputs/json",
        help="JSON 根目录"
    )
    p.add_argument(
        "--summary-out",
        default="example_outputs/qc/qscore_reparse_summary.csv",
        help="输出 reparse 汇总表"
    )
    p.add_argument(
        "--limit",
        type=int,
        default=0,
        help="只处理前 N 个 raw 文件，用于测试"
    )
    p.add_argument(
        "--entry",
        default="",
        help="只处理某一个 entry，例如 cf1a0a"
    )
    p.add_argument(
        "--min-coverage",
        type=float,
        default=0.2,
        help="coverage 低于该阈值时标记为 partial"
    )
    p.add_argument(
        "--verbose",
        action="store_true"
    )

    args = p.parse_args()

    raw_dir = Path(args.raw_dir)
    json_root = Path(args.json_root)
    summary_out = Path(args.summary_out)

    if args.entry:
        raw_files = [raw_dir / f"{args.entry}.qscore.raw.txt"]
    else:
        raw_files = sorted(raw_dir.glob("*.qscore.raw.txt"))

    if args.limit and args.limit > 0:
        raw_files = raw_files[:args.limit]

    print("========== QScore Raw Reparse ==========")
    print(f"Raw dir      : {raw_dir}")
    print(f"JSON root    : {json_root}")
    print(f"Summary out  : {summary_out}")
    print(f"Raw files    : {len(raw_files)}")
    print(f"Min coverage : {args.min_coverage}")
    print("========================================")

    rows = []
    total = Counter()

    for i, raw_path in enumerate(raw_files, start=1):
        print(f"[INFO] ({i}/{len(raw_files)}) Reparsing: {raw_path.name}")

        if not raw_path.exists():
            row = {
                "entry_id": entry_id_from_raw(raw_path),
                "raw_path": str(raw_path),
                "status": "failed",
                "chain_count": 0,
                "done": 0,
                "partial": 0,
                "failed": 0,
                "missing_chain_json": 0,
                "note": "raw file not found"
            }
        else:
            row = process_one_raw(
                raw_path=raw_path,
                json_root=json_root,
                min_coverage=args.min_coverage,
                verbose=args.verbose
            )

        rows.append(row)
        total[row["status"]] += 1

    write_summary_csv(summary_out, rows)

    print("========== Summary ==========")
    print(f"Processed raw files : {len(rows)}")
    print(f"Entry status        : {dict(total)}")
    print(f"Wrote summary       : {summary_out}")
    print("=============================")


if __name__ == "__main__":
    main()