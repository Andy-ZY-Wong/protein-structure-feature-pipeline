#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
01_build_skeleton.py

功能：
根据 manifest.csv，为每个蛋白结构条目生成每条链的 JSON 骨架。

输入：
- manifest.csv / manifest_with_fasta.csv
  必须至少包含：
    entry_id
    cfdb_id
    emdb_id
    pdb_id
    structure_path
    has_structure
  可选包含：
    fasta_path
    has_fasta

输出：
- output/json/{entry_id}/entry.json
- output/json/{entry_id}/chains/{chain_id}.json
- 更新后的 manifest_skeleton.csv

每条链 JSON 包含：
- entry_id / cfdb_id / emdb_id / pdb_id
- chain_id
- chain sequence
- sequence length
- sequence hash
- sequence group
- residue_mapping: chain_index -> pdb residue number
- FASTA 参考信息（如果有）
- tracks 空容器
- 各阶段状态
"""

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional

from Bio import SeqIO
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1


# 一些常见修饰氨基酸映射
CUSTOM_AA_MAP = {
    "MSE": "M",
    "SEC": "U",
    "PYL": "O",
    "UNK": "X",
}


def now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def safe_filename(x: str) -> str:
    """
    链 ID 有时可能包含特殊字符，用于文件名时做简单清理。
    """
    x = str(x)
    for ch in ["/", "\\", " ", ":", ";", ",", "|"]:
        x = x.replace(ch, "_")
    return x


def sha256_short(seq: str, n: int = 12) -> str:
    return hashlib.sha256(seq.encode("utf-8")).hexdigest()[:n]


def open_maybe_gzip_to_temp(path: str) -> Tuple[str, Optional[str]]:
    """
    Bio.PDB 对 gzip 文件兼容性不稳定，因此如果是 .gz，先解压到临时文件。
    返回：
      parse_path, temp_path
    如果不是 gz：
      temp_path = None
    """
    if not path.lower().endswith(".gz"):
        return path, None

    suffix = ""
    lower = path.lower()
    if lower.endswith(".pdb.gz"):
        suffix = ".pdb"
    elif lower.endswith(".cif.gz"):
        suffix = ".cif"
    elif lower.endswith(".ent.gz"):
        suffix = ".ent"
    else:
        suffix = ".tmp"

    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    tmp.close()

    with gzip.open(path, "rb") as fin, open(tmp.name, "wb") as fout:
        shutil.copyfileobj(fin, fout)

    return tmp.name, tmp.name


def parse_structure(structure_path: str, entry_id: str):
    """
    读取 PDB/mmCIF 结构。
    """
    parse_path, tmp_path = open_maybe_gzip_to_temp(structure_path)

    try:
        lower = parse_path.lower()
        if lower.endswith(".cif"):
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)

        structure = parser.get_structure(entry_id, parse_path)
        return structure, tmp_path

    except Exception:
        if tmp_path and os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise


def residue_to_one_letter(res) -> str:
    """
    把 Bio.PDB residue 转换成单字母氨基酸。
    """
    resname = res.get_resname().strip().upper()
    try:
        aa = seq1(resname, custom_map=CUSTOM_AA_MAP, undef_code="X")
    except Exception:
        aa = "X"

    if len(aa) != 1:
        aa = "X"
    return aa


def extract_protein_chains(structure) -> List[Dict]:
    """
    从结构中提取所有蛋白链。

    只读取第一个 model。
    每条链返回：
    {
      chain_id,
      sequence,
      length,
      residue_mapping
    }
    """

    chains_out = []

    models = list(structure.get_models())
    if not models:
        return chains_out

    model = models[0]

    for chain in model:
        chain_id = str(chain.id).strip()
        if not chain_id:
            chain_id = "_blank"

        seq_chars = []
        residue_mapping = []
        chain_index = 0

        for res in chain:
            # res.id = (hetfield, resseq, icode)
            # 只保留氨基酸残基。standard=False 可以保留 MSE 等修饰氨基酸。
            if not is_aa(res, standard=False):
                continue

            aa = residue_to_one_letter(res)
            chain_index += 1
            seq_chars.append(aa)

            hetfield, resseq, icode = res.get_id()

            residue_mapping.append({
                "chain_index": chain_index,
                "pdb_resseq": int(resseq),
                "pdb_icode": str(icode).strip(),
                "pdb_hetfield": str(hetfield).strip(),
                "resname": res.get_resname().strip(),
                "one_letter": aa
            })

        seq = "".join(seq_chars)

        # 跳过没有蛋白残基的链，比如水、配体链等
        if len(seq) == 0:
            continue

        chains_out.append({
            "chain_id": chain_id,
            "sequence": seq,
            "length": len(seq),
            "residue_mapping": residue_mapping,
        })

    return chains_out


def read_fasta_records(fasta_path: str) -> List[Dict]:
    """
    读取 FASTA，支持 .gz。
    返回多个 record。
    """

    if not fasta_path:
        return []

    if not os.path.exists(fasta_path):
        return []

    records = []

    if fasta_path.lower().endswith(".gz"):
        handle = gzip.open(fasta_path, "rt", encoding="utf-8", errors="ignore")
    else:
        handle = open(fasta_path, "r", encoding="utf-8", errors="ignore")

    with handle:
        for rec in SeqIO.parse(handle, "fasta"):
            seq = str(rec.seq).replace("*", "").upper()
            records.append({
                "id": rec.id,
                "description": rec.description,
                "sequence": seq,
                "length": len(seq),
                "sequence_hash": sha256_short(seq),
            })

    return records


def simple_fasta_match(chain_seq: str, fasta_records: List[Dict]) -> Dict:
    """
    为一条 PDB chain sequence 找最可能对应的 FASTA record。

    这里不做复杂 alignment，只做 skeleton 阶段的轻量匹配：
    1. 完全一致
    2. chain_seq 是 FASTA 的子串
    3. FASTA 是 chain_seq 的子串
    4. 如果只有一个 FASTA record，就作为候选参考
    5. 否则返回 unmatched

    真正精确的 UniProt / reference mapping 放到后续步骤做。
    """

    if not fasta_records:
        return {
            "status": "no_fasta",
            "method": "",
            "fasta_id": "",
            "fasta_description": "",
            "fasta_length": None,
            "reference_sequence": "",
            "note": "No FASTA file or no FASTA records."
        }

    chain_seq = chain_seq.upper()

    # 1. 完全一致
    for rec in fasta_records:
        if chain_seq == rec["sequence"]:
            return {
                "status": "matched",
                "method": "exact_sequence_match",
                "fasta_id": rec["id"],
                "fasta_description": rec["description"],
                "fasta_length": rec["length"],
                "reference_sequence": rec["sequence"],
                "note": ""
            }

    # 2. chain 是 fasta 子串
    for rec in fasta_records:
        if chain_seq in rec["sequence"]:
            return {
                "status": "matched",
                "method": "chain_is_subsequence_of_fasta",
                "fasta_id": rec["id"],
                "fasta_description": rec["description"],
                "fasta_length": rec["length"],
                "reference_sequence": rec["sequence"],
                "note": "PDB chain sequence is contained in FASTA sequence."
            }

    # 3. fasta 是 chain 子串
    for rec in fasta_records:
        if rec["sequence"] in chain_seq:
            return {
                "status": "matched",
                "method": "fasta_is_subsequence_of_chain",
                "fasta_id": rec["id"],
                "fasta_description": rec["description"],
                "fasta_length": rec["length"],
                "reference_sequence": rec["sequence"],
                "note": "FASTA sequence is contained in PDB chain sequence."
            }

    # 4. 如果只有一个 FASTA record，则保存为候选
    if len(fasta_records) == 1:
        rec = fasta_records[0]
        return {
            "status": "candidate",
            "method": "single_fasta_record_no_exact_match",
            "fasta_id": rec["id"],
            "fasta_description": rec["description"],
            "fasta_length": rec["length"],
            "reference_sequence": rec["sequence"],
            "note": "Only one FASTA record found, but no exact/subsequence match. Detailed alignment should be done later."
        }

    # 5. 多个 FASTA 且不能简单匹配
    return {
        "status": "unmatched",
        "method": "no_simple_match",
        "fasta_id": "",
        "fasta_description": "",
        "fasta_length": None,
        "reference_sequence": "",
        "note": "Multiple FASTA records found, but no simple exact/subsequence match."
    }


def assign_sequence_groups(chains: List[Dict]) -> Dict[str, str]:
    """
    按 chain sequence 分组。
    返回：
      sequence_hash -> group_id
    """

    unique_hashes = []
    seen = set()

    for ch in chains:
        h = sha256_short(ch["sequence"])
        if h not in seen:
            seen.add(h)
            unique_hashes.append(h)

    hash_to_group = {}
    for i, h in enumerate(unique_hashes, start=1):
        hash_to_group[h] = f"group_{i:03d}"

    return hash_to_group


def build_chain_json(
    row: Dict,
    chain: Dict,
    sequence_group: str,
    fasta_match: Dict
) -> Dict:
    """
    构建单条链 JSON 骨架。
    """

    entry_id = row.get("entry_id", "")
    cfdb_id = row.get("cfdb_id", "")
    emdb_id = row.get("emdb_id", "")
    pdb_id = row.get("pdb_id", "")

    chain_seq = chain["sequence"]
    chain_hash = sha256_short(chain_seq)

    chain_json = {
        "schema_version": "sequence_tracks_v1",
        "generated_at": now_iso(),

        "entry": {
            "entry_id": entry_id,
            "cfdb_id": cfdb_id,
            "emdb_id": emdb_id,
            "pdb_id": pdb_id,
        },

        "source_files": {
            "structure_path": row.get("structure_path", ""),
            "fasta_path": row.get("fasta_path", ""),
        },

        "chain": {
            "chain_id": chain["chain_id"],
            "sequence": chain_seq,
            "length": chain["length"],
            "sequence_hash": chain_hash,
            "sequence_group": sequence_group,
            "coordinate_system": "chain_index_1_based"
        },

        "reference_fasta": fasta_match,

        "residue_mapping": chain["residue_mapping"],

        "tracks": {},

        "computed_features": {},

        "status": {
            "skeleton": "done",
            "local_tracks": "pending",
            "uniprot_mapping": "pending",
            "external_annotations": "pending",
            "frontend_tracks": "pending"
        },

        "logs": [
            {
                "time": now_iso(),
                "step": "skeleton",
                "status": "done",
                "message": "Skeleton JSON generated from structure file."
            }
        ]
    }

    return chain_json


def write_json(path: str, data: Dict) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)


def process_one_entry(row: Dict, out_dir: str, overwrite: bool = False) -> Tuple[str, Dict]:
    """
    处理 manifest 中的一行。
    返回：
      status, updated_row
    """

    updated = dict(row)

    entry_id = row.get("entry_id", "").strip()
    structure_path = row.get("structure_path", "").strip()
    has_structure = row.get("has_structure", "").strip().lower()

    if not entry_id:
        updated["skeleton_status"] = "failed"
        updated["overall_status"] = "failed"
        updated["note"] = append_note(updated.get("note", ""), "Missing entry_id")
        return "failed", updated

    if has_structure != "yes" or not structure_path or not os.path.exists(structure_path):
        updated["skeleton_status"] = "failed"
        updated["overall_status"] = "missing_structure"
        updated["note"] = append_note(updated.get("note", ""), "Structure file missing")
        return "failed", updated

    entry_out_dir = os.path.join(out_dir, entry_id)
    chains_out_dir = os.path.join(entry_out_dir, "chains")
    entry_json_path = os.path.join(entry_out_dir, "entry.json")

    if os.path.exists(entry_json_path) and not overwrite:
        updated["skeleton_status"] = "done"
        updated["overall_status"] = updated.get("overall_status", "pending") or "pending"
        updated["note"] = append_note(updated.get("note", ""), "Skeleton already exists; skipped")
        return "skipped", updated

    tmp_path = None

    try:
        structure, tmp_path = parse_structure(structure_path, entry_id)
        chains = extract_protein_chains(structure)

        if tmp_path and os.path.exists(tmp_path):
            os.remove(tmp_path)

        if not chains:
            updated["skeleton_status"] = "failed"
            updated["overall_status"] = "failed"
            updated["note"] = append_note(updated.get("note", ""), "No protein chains found")
            return "failed", updated

        # 读取 FASTA，如果有
        fasta_path = row.get("fasta_path", "").strip()
        fasta_records = read_fasta_records(fasta_path) if fasta_path else []

        # 分组
        hash_to_group = assign_sequence_groups(chains)

        chain_summaries = []
        groups_summary: Dict[str, Dict] = {}

        for chain in chains:
            chain_id = chain["chain_id"]
            chain_safe = safe_filename(chain_id)

            h = sha256_short(chain["sequence"])
            group_id = hash_to_group[h]

            fasta_match = simple_fasta_match(chain["sequence"], fasta_records)

            chain_json = build_chain_json(
                row=row,
                chain=chain,
                sequence_group=group_id,
                fasta_match=fasta_match
            )

            chain_json_path = os.path.join(chains_out_dir, f"{chain_safe}.json")
            write_json(chain_json_path, chain_json)

            chain_summaries.append({
                "chain_id": chain_id,
                "length": chain["length"],
                "sequence_hash": h,
                "sequence_group": group_id,
                "json_path": chain_json_path,
                "fasta_match_status": fasta_match.get("status", "")
            })

            if group_id not in groups_summary:
                groups_summary[group_id] = {
                    "sequence_hash": h,
                    "length": chain["length"],
                    "chains": []
                }
            groups_summary[group_id]["chains"].append(chain_id)

        entry_json = {
            "schema_version": "sequence_tracks_entry_v1",
            "generated_at": now_iso(),
            "entry": {
                "entry_id": row.get("entry_id", ""),
                "cfdb_id": row.get("cfdb_id", ""),
                "emdb_id": row.get("emdb_id", ""),
                "pdb_id": row.get("pdb_id", ""),
            },
            "source_files": {
                "structure_path": structure_path,
                "fasta_path": fasta_path,
            },
            "chain_count": len(chains),
            "chains": chain_summaries,
            "sequence_groups": [
                {
                    "group_id": gid,
                    "sequence_hash": info["sequence_hash"],
                    "length": info["length"],
                    "chains": info["chains"]
                }
                for gid, info in groups_summary.items()
            ],
            "status": {
                "skeleton": "done",
                "local_tracks": "pending",
                "uniprot_mapping": "pending",
                "external_annotations": "pending",
                "frontend_tracks": "pending"
            },
            "logs": [
                {
                    "time": now_iso(),
                    "step": "skeleton",
                    "status": "done",
                    "message": f"Generated skeleton for {len(chains)} chains."
                }
            ]
        }

        write_json(entry_json_path, entry_json)

        updated["skeleton_status"] = "done"
        if not updated.get("overall_status") or updated.get("overall_status") in ["pending", "missing_structure"]:
            updated["overall_status"] = "pending"
        updated["note"] = append_note(updated.get("note", ""), f"Skeleton generated; chains={len(chains)}")
        return "done", updated

    except Exception as e:
        if tmp_path and os.path.exists(tmp_path):
            os.remove(tmp_path)

        updated["skeleton_status"] = "failed"
        updated["overall_status"] = "failed"
        updated["note"] = append_note(updated.get("note", ""), f"Skeleton error: {repr(e)}")
        return "failed", updated


def append_note(old: str, new: str) -> str:
    old = old or ""
    if not old:
        return new
    if new in old:
        return old
    return old + "; " + new


def build_skeleton_from_manifest(
    manifest_path: str,
    out_dir: str,
    manifest_out: str,
    limit: int = 0,
    overwrite: bool = False
) -> None:
    """
    主流程。
    """

    if not os.path.exists(manifest_path):
        raise FileNotFoundError(f"Cannot find manifest: {manifest_path}")

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.dirname(manifest_out) or ".", exist_ok=True)

    rows = []
    with open(manifest_path, "r", encoding="utf-8", errors="ignore", newline="") as fin:
        reader = csv.DictReader(fin)
        if reader.fieldnames is None:
            raise RuntimeError(f"No header found in manifest: {manifest_path}")
        fieldnames = list(reader.fieldnames)
        for row in reader:
            rows.append(row)

    if limit and limit > 0:
        rows_to_process = rows[:limit]
        rows_remaining = rows[limit:]
    else:
        rows_to_process = rows
        rows_remaining = []

    stats = {
        "done": 0,
        "skipped": 0,
        "failed": 0,
    }

    updated_rows = []

    for i, row in enumerate(rows_to_process, start=1):
        entry_id = row.get("entry_id", "")
        print(f"[INFO] ({i}/{len(rows_to_process)}) Processing entry: {entry_id}")

        status, updated = process_one_entry(
            row=row,
            out_dir=out_dir,
            overwrite=overwrite
        )

        stats[status] = stats.get(status, 0) + 1
        updated_rows.append(updated)

    # limit 模式下，剩余的原样保留
    updated_rows.extend(rows_remaining)

    # 保证新增字段存在
    out_fields = list(fieldnames)
    for col in ["skeleton_status", "overall_status", "note"]:
        if col not in out_fields:
            out_fields.append(col)

    with open(manifest_out, "w", encoding="utf-8", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(updated_rows)

    print("========== Skeleton Summary ==========")
    print(f"Manifest input       : {manifest_path}")
    print(f"Output JSON dir      : {out_dir}")
    print(f"Manifest output      : {manifest_out}")
    print(f"Entries processed    : {len(rows_to_process)}")
    print(f"Done                 : {stats.get('done', 0)}")
    print(f"Skipped              : {stats.get('skipped', 0)}")
    print(f"Failed               : {stats.get('failed', 0)}")
    print("======================================")


def main():
    parser = argparse.ArgumentParser(
        description="Build per-entry/per-chain JSON skeletons from manifest.csv."
    )

    parser.add_argument(
        "--manifest",
        required=True,
        help="manifest.csv 或 manifest_with_fasta.csv 路径"
    )

    parser.add_argument(
        "--out-dir",
        required=True,
        help="输出 JSON 根目录，例如 example_outputs/json"
    )

    parser.add_argument(
        "--manifest-out",
        required=True,
        help="输出更新后的 manifest，例如 manifest_skeleton.csv"
    )

    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="只处理前 N 行，用于测试；默认 0 表示处理全部"
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="如果 JSON 已存在，是否覆盖重建"
    )

    args = parser.parse_args()

    build_skeleton_from_manifest(
        manifest_path=args.manifest,
        out_dir=args.out_dir,
        manifest_out=args.manifest_out,
        limit=args.limit,
        overwrite=args.overwrite
    )


if __name__ == "__main__":
    main()