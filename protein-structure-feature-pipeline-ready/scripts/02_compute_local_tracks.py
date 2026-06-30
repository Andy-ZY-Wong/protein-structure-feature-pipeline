#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
02_compute_local_tracks.py

功能：
读取 01_build_skeleton.py 生成的 chain-level JSON，
批量计算本地序列/结构指标并填入 JSON 的 computed_features 字段。

本脚本计算：
1. Secondary Structure（二级结构）：DSSP / mkdssp
2. Buried Residues / RSA：freesasa-python
3. Hydropathy：Kyte-Doolittle
4. Disorder：可选 IUPred2A；如果没有工具则标记 skipped
5. Atomic clashes：内部 VDW overlap 近似算法
6. Bond outliers：内部主链键长异常近似算法
7. Cis-peptide：omega 二面角计算

说明：
- Atomic clashes 和 Bond outliers 是本地近似版本，不等同于 wwPDB/RCSB 官方 validation report。
- Q-score 不在本脚本中计算，后续单独基于 map + PDB 计算后再补入 JSON。
- 批处理按 manifest 行分组，默认 3000 个 entry 为一批。
- batch-index 从 0 开始：
    batch-index 0 -> 第 1 到 3000 条
    batch-index 1 -> 第 3001 到 6000 条
"""

import argparse
import csv
import gzip
import json
import math
import os
import shutil
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np

from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch
from Bio.PDB.DSSP import DSSP
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.vectors import calc_dihedral


# ============================================================
# 基础工具函数
# ============================================================

def now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def read_json(path: str) -> Dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def write_json(path: str, data: Dict) -> None:
    tmp = path + ".tmp"
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


def append_note(old: str, new: str) -> str:
    old = old or ""
    if not old:
        return new
    if new in old:
        return old
    return old + "; " + new


def safe_float(x, ndigits: int = 4):
    if x is None:
        return None
    try:
        x = float(x)
        if math.isnan(x) or math.isinf(x):
            return None
        return round(x, ndigits)
    except Exception:
        return None


def open_maybe_gzip_to_temp(path: str) -> Tuple[str, Optional[str]]:
    """
    如果结构文件是 .gz，临时解压。
    返回 parse_path, temp_path。
    """
    if not path.lower().endswith(".gz"):
        return path, None

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


def cleanup_temp(temp_path: Optional[str]) -> None:
    if temp_path and os.path.exists(temp_path):
        os.remove(temp_path)


def parse_structure(structure_path: str, entry_id: str):
    """
    读取 PDB/mmCIF 结构。
    返回 structure, parse_path, temp_path。
    """
    parse_path, temp_path = open_maybe_gzip_to_temp(structure_path)
    lower = parse_path.lower()

    if lower.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure(entry_id, parse_path)
    return structure, parse_path, temp_path


def get_first_model(structure):
    models = list(structure.get_models())
    if not models:
        return None
    return models[0]


def residue_key_from_mapping_item(x: Dict) -> Tuple[str, int, str]:
    het = str(x.get("pdb_hetfield", "")).strip()
    resseq = int(x.get("pdb_resseq"))
    icode = str(x.get("pdb_icode", "")).strip()
    return het, resseq, icode


def build_residue_key_to_chain_index(chain_json: Dict) -> Dict[Tuple[str, int, str], int]:
    out = {}
    for item in chain_json.get("residue_mapping", []):
        key = residue_key_from_mapping_item(item)
        out[key] = int(item["chain_index"])
    return out


# ============================================================
# 1. Hydropathy：Kyte-Doolittle
# ============================================================

KD = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8,
    "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6,
    "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
    "X": 0.0, "U": 0.0, "O": 0.0
}


def compute_hydropathy(seq: str, window: int = 9) -> List[float]:
    vals = np.array([KD.get(aa.upper(), 0.0) for aa in seq], dtype=float)

    w = max(3, int(window))
    if w % 2 == 0:
        w += 1

    pad = w // 2
    padded = np.pad(vals, (pad, pad), mode="edge")
    out = np.convolve(padded, np.ones(w) / w, mode="valid")

    return [safe_float(v, 4) for v in out]


# ============================================================
# 2. Secondary Structure：DSSP
# ============================================================

def normalize_dssp_ss(ss: str) -> str:
    """
    DSSP 原始二级结构符号：
      H = alpha helix
      B = beta bridge
      E = extended strand
      G = 3-10 helix
      I = pi-helix
      T = turn
      S = bend
      blank = coil / other

    统一成：
      H = helix
      E = strand
      C = coil / other
    """
    ss = (ss or "").strip()

    if ss in {"H", "G", "I"}:
        return "H"

    if ss in {"E", "B"}:
        return "E"

    return "C"


def segments_from_position_labels(pos_to_label: Dict[int, str], length: int) -> List[Dict]:
    """
    把每个残基位置的 H/E/C 标签合并成连续区间。

    输入：
      pos_to_label:
        {
          1: "C",
          2: "H",
          3: "H",
          ...
        }

    输出：
      [
        {"start": 2, "end": 10, "label": "H"},
        {"start": 25, "end": 30, "label": "E"}
      ]

    注意：
      C 不输出，只输出 H/E，方便前端画轨道。
    """

    items = []

    current_label = None
    start = None
    prev = None

    for pos in range(1, length + 1):
        label = pos_to_label.get(pos, "C")

        if label == "C":
            if current_label in {"H", "E"}:
                items.append({
                    "start": start,
                    "end": prev,
                    "label": current_label
                })

            current_label = None
            start = None
            prev = None
            continue

        if current_label is None:
            current_label = label
            start = pos
            prev = pos

        elif label == current_label and pos == prev + 1:
            prev = pos

        else:
            items.append({
                "start": start,
                "end": prev,
                "label": current_label
            })

            current_label = label
            start = pos
            prev = pos

    if current_label in {"H", "E"}:
        items.append({
            "start": start,
            "end": prev,
            "label": current_label
        })

    return items


def make_no_seqres_pdb_for_dssp(original_pdb: str) -> str:
    """
    生成一个去掉 SEQRES 记录的临时 PDB，专门用于 DSSP。

    背景：
      有些 PDB 文件中 SEQRES 和 ATOM 记录不一致。
      例如 cf100j.pdb 中 mkdssp 报错：
        A residue found in the ATOM records was not found in the SEQRES records

    处理：
      - 删除 SEQRES 行；
      - 保留 HEADER / TITLE / ATOM / HETATM / TER / END 等其它内容；
      - 不修改原始 PDB；
      - 返回临时 PDB 路径。
    """

    tmp = tempfile.NamedTemporaryFile(
        delete=False,
        suffix=".noseqres.pdb",
        mode="w"
    )

    tmp_path = tmp.name

    with open(original_pdb, "r", encoding="utf-8", errors="ignore") as fin, tmp:
        for line in fin:
            if line.startswith("SEQRES"):
                continue

            tmp.write(line)

    return tmp_path


def run_mkdssp_to_file(input_pdb: str, dssp_exe: str = "auto") -> str:
    """
    直接调用 mkdssp 生成 DSSP 文件。

    你的环境里 mkdssp 不支持：
      mkdssp -i input.pdb -o output.dssp

    当前已验证可用的方式是：
      mkdssp input.pdb output.dssp

    本函数会优先尝试：
      1. mkdssp input.pdb output.dssp
      2. mkdssp input.pdb > output.dssp

    返回：
      dssp 输出文件路径。
    """

    if dssp_exe and dssp_exe != "auto":
        exe = dssp_exe
    else:
        exe = shutil.which("mkdssp") or shutil.which("dssp")

    if exe is None:
        raise RuntimeError("Cannot find mkdssp or dssp executable.")

    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".dssp")
    dssp_out = tmp.name
    tmp.close()

    commands = [
        [exe, input_pdb, dssp_out],
        [exe, input_pdb],
    ]

    last_error = None

    for cmd in commands:
        try:
            proc = subprocess.run(
                cmd,
                text=True,
                capture_output=True,
                timeout=300
            )

            # 情况 1：mkdssp input.pdb output.dssp
            if proc.returncode == 0:
                if os.path.exists(dssp_out) and os.path.getsize(dssp_out) > 0:
                    return dssp_out

                # 情况 2：mkdssp input.pdb，把结果输出到 stdout
                if proc.stdout and proc.stdout.strip():
                    with open(dssp_out, "w", encoding="utf-8") as f:
                        f.write(proc.stdout)

                    if os.path.getsize(dssp_out) > 0:
                        return dssp_out

            last_error = (
                f"cmd={' '.join(cmd)}; "
                f"returncode={proc.returncode}; "
                f"stdout={proc.stdout[:500] if proc.stdout else ''}; "
                f"stderr={proc.stderr[:1000] if proc.stderr else ''}"
            )

        except Exception as e:
            last_error = repr(e)

    if os.path.exists(dssp_out):
        os.remove(dssp_out)

    raise RuntimeError(f"mkdssp failed. Last error: {last_error}")


def parse_dssp_file_to_chain_results(
    dssp_path: str,
    chain_jsons: Dict[str, Dict]
) -> Dict[str, Dict]:
    """
    解析 mkdssp 生成的 .dssp 文件，并映射回 chain_index 坐标。

    DSSP 经典文本格式中，数据区一般在这一行之后：
      "  #  RESIDUE AA STRUCTURE ..."

    常用固定列：
      line[5:10]  = PDB residue number
      line[10]    = insertion code
      line[11]    = chain ID
      line[13]    = amino acid
      line[16]    = secondary structure symbol

    输出：
      {
        "A": {
          "per_position": ["C", "H", "H", ...],
          "segments": [
            {"start": 2, "end": 10, "label": "H"}
          ]
        }
      }
    """

    result = {}

    for chain_id, cj in chain_jsons.items():
        L = int(cj["chain"]["length"])
        result[chain_id] = {
            "per_position": ["C"] * L,
            "segments": []
        }

    key_maps = {
        chain_id: build_residue_key_to_chain_index(cj)
        for chain_id, cj in chain_jsons.items()
    }

    in_data = False

    with open(dssp_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("  #  RESIDUE"):
                in_data = True
                continue

            if not in_data:
                continue

            if len(line) < 17:
                continue

            # DSSP 中 ! 通常表示链断裂或缺失残基
            if len(line) > 13 and line[13] == "!":
                continue

            try:
                resseq_text = line[5:10].strip()

                if not resseq_text:
                    continue

                resseq = int(resseq_text)
                icode = line[10].strip()
                chain_id = line[11].strip()
                ss_raw = line[16].strip()

            except Exception:
                continue

            if chain_id not in chain_jsons:
                continue

            ss = normalize_dssp_ss(ss_raw)

            # 优先匹配普通 ATOM 残基：hetfield 为空
            key = ("", resseq, icode)
            chain_index = key_maps[chain_id].get(key)

            # 如果没匹配到，放宽条件：只按 resseq + icode 匹配
            if chain_index is None:
                for k, idx in key_maps[chain_id].items():
                    if k[1] == resseq and k[2] == icode:
                        chain_index = idx
                        break

            if chain_index is None:
                continue

            result[chain_id]["per_position"][chain_index - 1] = ss

    for chain_id, cj in chain_jsons.items():
        L = int(cj["chain"]["length"])

        pos_to_label = {
            i + 1: ss
            for i, ss in enumerate(result[chain_id]["per_position"])
        }

        result[chain_id]["segments"] = segments_from_position_labels(pos_to_label, L)

    return result


def compute_dssp_for_entry(
    structure,
    parse_path: str,
    chain_jsons: Dict[str, Dict],
    dssp_exe: str = "auto"
) -> Dict[str, Dict]:
    """
    对一个 entry 跑一次 DSSP，然后拆分到每条链。

    当前策略：
      1. 先用原始 PDB 跑 mkdssp；
      2. 如果失败，生成去掉 SEQRES 的临时 PDB；
      3. 用 no-SEQRES PDB 再跑 mkdssp；
      4. 解析 .dssp 文本结果；
      5. 映射回每条 chain JSON 的 chain_index 坐标。

    这样可以解决 cf100j 这类：
      SEQRES 与 ATOM 不一致导致 mkdssp 失败
    的问题。
    """

    dssp_path = None
    noseqres_pdb = None

    try:
        # 1. 优先尝试原始结构文件
        try:
            dssp_path = run_mkdssp_to_file(
                input_pdb=parse_path,
                dssp_exe=dssp_exe
            )

        except Exception as e1:
            print(f"[WARN] DSSP failed on original file: {repr(e1)}")
            print("[INFO] Trying no-SEQRES PDB fallback for DSSP...")

            # 2. 如果原始 PDB 失败，删除 SEQRES 后重试
            noseqres_pdb = make_no_seqres_pdb_for_dssp(parse_path)

            dssp_path = run_mkdssp_to_file(
                input_pdb=noseqres_pdb,
                dssp_exe=dssp_exe
            )

            print(f"[INFO] DSSP fallback succeeded with no-SEQRES PDB: {noseqres_pdb}")

        # 3. 解析 DSSP 输出并映射回 chain_index
        result = parse_dssp_file_to_chain_results(
            dssp_path=dssp_path,
            chain_jsons=chain_jsons
        )

        return result

    finally:
        if dssp_path and os.path.exists(dssp_path):
            os.remove(dssp_path)

        if noseqres_pdb and os.path.exists(noseqres_pdb):
            os.remove(noseqres_pdb)

# ============================================================
# 3. RSA / Buried Residues：FreeSASA
# ============================================================

def compute_freesasa_for_entry(
    structure_path: str,
    chain_jsons: Dict[str, Dict],
    buried_threshold: float = 0.20
) -> Dict[str, Dict]:

    try:
        import freesasa
    except Exception as e:
        raise RuntimeError(f"Cannot import freesasa: {repr(e)}")

    fs_structure = freesasa.Structure(structure_path)
    fs_result = freesasa.calc(fs_structure)
    areas = fs_result.residueAreas()

    out = {}

    for chain_id, cj in chain_jsons.items():
        L = int(cj["chain"]["length"])
        out[chain_id] = {
            "values": [None] * L,
            "buried_positions": [],
            "mapped_count": 0
        }

    # FreeSASA 通常按 chain + residue number 返回
    resnum_maps = {}

    for chain_id, cj in chain_jsons.items():
        m = {}
        for item in cj.get("residue_mapping", []):
            resseq = int(item["pdb_resseq"])
            chain_index = int(item["chain_index"])
            # 对插入码的处理简化为同 resseq 取首次出现
            m.setdefault(resseq, chain_index)
        resnum_maps[chain_id] = m

    for chain_id, res_map in areas.items():
        chain_id = str(chain_id).strip()

        if chain_id not in out:
            continue

        for resnum, area in res_map.items():
            try:
                resseq = int(str(resnum).strip())
            except Exception:
                continue

            chain_index = resnum_maps.get(chain_id, {}).get(resseq)

            if chain_index is None:
                continue

            rel = getattr(area, "relativeTotal", None)

            if rel is None:
                continue

            rel = float(rel)
            out[chain_id]["values"][chain_index - 1] = safe_float(rel, 4)
            out[chain_id]["mapped_count"] += 1

            if rel < buried_threshold:
                out[chain_id]["buried_positions"].append(chain_index)

    return out


# ============================================================
# 4. Disorder：IUPred2A 可选
# ============================================================

def run_iupred2a(seq: str, iupred_cmd: str, mode: str = "long") -> Optional[List[float]]:
    """
    运行 IUPred2A。
    命令示例：
      --iupred-cmd "python /path/to/iupred2a.py"

    实际执行：
      python /path/to/iupred2a.py temp.fasta long
    """

    if not iupred_cmd:
        return None

    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta", mode="w")
    tmp.write(">query\n")
    tmp.write(seq + "\n")
    tmp.close()

    try:
        cmd = iupred_cmd.split() + [tmp.name, mode]

        proc = subprocess.run(
            cmd,
            text=True,
            capture_output=True,
            timeout=180
        )

        if proc.returncode != 0:
            return None

        scores = []

        for line in proc.stdout.splitlines():
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split()

            if len(parts) < 3:
                continue

            try:
                pos = int(parts[0])
                score = float(parts[2])
            except Exception:
                continue

            scores.append((pos, score))

        if not scores:
            return None

        scores = sorted(scores, key=lambda x: x[0])
        values = [safe_float(s, 4) for _, s in scores]

        if len(values) != len(seq):
            return None

        return values

    finally:
        if os.path.exists(tmp.name):
            os.remove(tmp.name)


# ============================================================
# 5. Atomic Clashes：近似 VDW overlap
# ============================================================

VDW_RADII = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "P": 1.80,
    "F": 1.47,
    "CL": 1.75,
    "BR": 1.85,
    "I": 1.98,
}


def get_atom_element(atom) -> str:
    e = getattr(atom, "element", "") or ""
    e = str(e).strip().upper()

    if e:
        return e

    name = atom.get_name().strip().upper()

    if len(name) >= 2 and name[:2] in VDW_RADII:
        return name[:2]

    if len(name) >= 1:
        return name[0]

    return "C"


def atom_radius(atom) -> float:
    return VDW_RADII.get(get_atom_element(atom), 1.70)


def is_hydrogen(atom) -> bool:
    return get_atom_element(atom) == "H"


def atom_chain_id(atom) -> str:
    return str(atom.get_parent().get_parent().id).strip()


def atom_residue_key(atom) -> Tuple[str, int, str]:
    res = atom.get_parent()
    het, resseq, icode = res.get_id()
    return str(het).strip(), int(resseq), str(icode).strip()


def compute_atomic_clashes_approx(
    structure,
    chain_jsons: Dict[str, Dict],
    overlap_cutoff: float = 0.40,
    search_radius: float = 4.0
) -> Dict[str, Dict]:
    """
    近似 Atomic clashes：
    如果 overlap = r1 + r2 - distance >= overlap_cutoff，则认为有 clash。

    注意：
    - 排除氢原子
    - 排除同一残基
    - 排除同一链相邻残基，减少共价键误判
    - 这不是官方 MolProbity clashscore
    """

    model = get_first_model(structure)

    if model is None:
        raise RuntimeError("No model found.")

    atoms = [a for a in model.get_atoms() if not is_hydrogen(a)]

    if not atoms:
        raise RuntimeError("No atoms found.")

    ns = NeighborSearch(atoms)

    out = {}

    for chain_id, cj in chain_jsons.items():
        out[chain_id] = {
            "positions": set(),
            "items": [],
            "clash_pair_count": 0
        }

    res_maps = {}

    for chain_id, cj in chain_jsons.items():
        m = {}
        for item in cj.get("residue_mapping", []):
            key = residue_key_from_mapping_item(item)
            m[key] = int(item["chain_index"])
        res_maps[chain_id] = m

    seen_pairs = set()

    for atom1 in atoms:
        close_atoms = ns.search(atom1.get_coord(), search_radius, level="A")

        for atom2 in close_atoms:
            if atom1 is atom2:
                continue

            pair_key = tuple(sorted((id(atom1), id(atom2))))

            if pair_key in seen_pairs:
                continue

            seen_pairs.add(pair_key)

            ch1 = atom_chain_id(atom1)
            ch2 = atom_chain_id(atom2)
            key1 = atom_residue_key(atom1)
            key2 = atom_residue_key(atom2)

            # 排除同一残基
            if ch1 == ch2 and key1 == key2:
                continue

            # 排除同一链相邻残基，减少共价连接误判
            if ch1 == ch2:
                if abs(key1[1] - key2[1]) <= 1:
                    continue

            d = float(np.linalg.norm(atom1.get_coord() - atom2.get_coord()))
            overlap = atom_radius(atom1) + atom_radius(atom2) - d

            if overlap < overlap_cutoff:
                continue

            for ch, key, atom, other_atom in [
                (ch1, key1, atom1, atom2),
                (ch2, key2, atom2, atom1),
            ]:
                if ch not in out:
                    continue

                chain_index = res_maps.get(ch, {}).get(key)

                if chain_index is None:
                    key_alt = ("", key[1], key[2])
                    chain_index = res_maps.get(ch, {}).get(key_alt)

                if chain_index is None:
                    continue

                other_key = atom_residue_key(other_atom)

                out[ch]["positions"].add(chain_index)
                out[ch]["items"].append({
                    "pos": chain_index,
                    "pdb_resseq": key[1],
                    "atom": atom.get_name().strip(),
                    "other_chain": atom_chain_id(other_atom),
                    "other_pdb_resseq": other_key[1],
                    "other_atom": other_atom.get_name().strip(),
                    "distance": safe_float(d, 3),
                    "overlap": safe_float(overlap, 3)
                })
                out[ch]["clash_pair_count"] += 1

    for ch in out:
        out[ch]["positions"] = sorted(out[ch]["positions"])

    return out


# ============================================================
# 6. Bond Outliers：近似主链键长异常
# ============================================================

BACKBONE_BONDS = [
    ("N", "CA", 1.459),
    ("CA", "C", 1.525),
    ("C", "O", 1.229),
]


def atom_by_name(residue, atom_name: str):
    atom_name = atom_name.strip()
    if atom_name in residue:
        return residue[atom_name]
    return None


def compute_bond_outliers_approx(
    structure,
    chain_jsons: Dict[str, Dict],
    tolerance: float = 0.08,
    peptide_bond_max_distance: float = 1.9
) -> Dict[str, Dict]:
    """
    近似 Bond outliers：
    - 检查 N-CA、CA-C、C-O
    - 检查相邻残基 C(prev)-N(current)
    - 如果偏离 ideal bond length 超过 tolerance，则标记 outlier
    - 如果 C-N 距离大于 peptide_bond_max_distance，认为不是连续肽键，跳过，避免缺失残基误报
    """

    model = get_first_model(structure)

    if model is None:
        raise RuntimeError("No model found.")

    out = {}

    for chain_id, cj in chain_jsons.items():
        out[chain_id] = {
            "positions": set(),
            "items": [],
            "outlier_count": 0
        }

    for chain_id, cj in chain_jsons.items():
        if chain_id not in model:
            continue

        chain = model[chain_id]

        resseq_to_index = {}

        for item in cj.get("residue_mapping", []):
            resseq_to_index[int(item["pdb_resseq"])] = int(item["chain_index"])

        residues = [res for res in chain if is_aa(res, standard=False)]
        residues = sorted(
            residues,
            key=lambda r: (int(r.get_id()[1]), str(r.get_id()[2]).strip())
        )

        prev_res = None

        for res in residues:
            het, resseq, icode = res.get_id()
            resseq = int(resseq)
            chain_index = resseq_to_index.get(resseq)

            if chain_index is None:
                prev_res = res
                continue

            # 残基内主链键
            for a1_name, a2_name, ideal in BACKBONE_BONDS:
                a1 = atom_by_name(res, a1_name)
                a2 = atom_by_name(res, a2_name)

                if a1 is None or a2 is None:
                    continue

                d = float(np.linalg.norm(a1.get_coord() - a2.get_coord()))
                delta = abs(d - ideal)

                if delta > tolerance:
                    out[chain_id]["positions"].add(chain_index)
                    out[chain_id]["items"].append({
                        "pos": chain_index,
                        "pdb_resseq": resseq,
                        "bond": f"{a1_name}-{a2_name}",
                        "observed": safe_float(d, 3),
                        "ideal": ideal,
                        "delta": safe_float(delta, 3)
                    })
                    out[chain_id]["outlier_count"] += 1

            # 相邻残基 peptide bond: C(prev)-N(current)
            if prev_res is not None:
                prev_c = atom_by_name(prev_res, "C")
                cur_n = atom_by_name(res, "N")

                if prev_c is not None and cur_n is not None:
                    cn_dist = float(np.linalg.norm(prev_c.get_coord() - cur_n.get_coord()))

                    # 跳过明显不是连续肽键的情况
                    if cn_dist <= peptide_bond_max_distance:
                        ideal = 1.329
                        delta = abs(cn_dist - ideal)

                        if delta > tolerance:
                            out[chain_id]["positions"].add(chain_index)
                            out[chain_id]["items"].append({
                                "pos": chain_index,
                                "pdb_resseq": resseq,
                                "bond": "C(prev)-N",
                                "observed": safe_float(cn_dist, 3),
                                "ideal": ideal,
                                "delta": safe_float(delta, 3)
                            })
                            out[chain_id]["outlier_count"] += 1

            prev_res = res

    for ch in out:
        out[ch]["positions"] = sorted(out[ch]["positions"])

    return out


# ============================================================
# 7. Cis-peptide：omega 二面角
# ============================================================

def compute_cis_peptides(
    structure,
    chain_jsons: Dict[str, Dict],
    cis_cutoff_degree: float = 30.0,
    peptide_bond_max_distance: float = 1.9
) -> Dict[str, Dict]:
    """
    计算 cis-peptide。

    omega 定义：
      CA(i-1) - C(i-1) - N(i) - CA(i)

    判定：
      |omega| <= cis_cutoff_degree
      且 C(i-1)-N(i) 距离 <= peptide_bond_max_distance

    输出位置使用当前残基 i 的 chain_index。
    """

    model = get_first_model(structure)

    if model is None:
        raise RuntimeError("No model found.")

    out = {}

    for chain_id, cj in chain_jsons.items():
        out[chain_id] = {
            "positions": [],
            "items": [],
            "cis_count": 0
        }

        if chain_id not in model:
            continue

        chain = model[chain_id]

        resseq_to_index = {}

        for item in cj.get("residue_mapping", []):
            resseq_to_index[int(item["pdb_resseq"])] = int(item["chain_index"])

        residues = [res for res in chain if is_aa(res, standard=False)]
        residues = sorted(
            residues,
            key=lambda r: (int(r.get_id()[1]), str(r.get_id()[2]).strip())
        )

        for prev_res, cur_res in zip(residues[:-1], residues[1:]):
            prev_resseq = int(prev_res.get_id()[1])
            cur_resseq = int(cur_res.get_id()[1])

            prev_idx = resseq_to_index.get(prev_resseq)
            cur_idx = resseq_to_index.get(cur_resseq)

            if prev_idx is None or cur_idx is None:
                continue

            if "CA" not in prev_res or "C" not in prev_res:
                continue

            if "N" not in cur_res or "CA" not in cur_res:
                continue

            prev_ca = prev_res["CA"]
            prev_c = prev_res["C"]
            cur_n = cur_res["N"]
            cur_ca = cur_res["CA"]

            cn_dist = float(np.linalg.norm(prev_c.get_coord() - cur_n.get_coord()))

            if cn_dist > peptide_bond_max_distance:
                continue

            omega_rad = calc_dihedral(
                prev_ca.get_vector(),
                prev_c.get_vector(),
                cur_n.get_vector(),
                cur_ca.get_vector()
            )

            omega_deg = math.degrees(float(omega_rad))

            if abs(omega_deg) <= cis_cutoff_degree:
                out[chain_id]["positions"].append(cur_idx)
                out[chain_id]["items"].append({
                    "pos": cur_idx,
                    "prev_pos": prev_idx,
                    "pdb_resseq": cur_resseq,
                    "prev_pdb_resseq": prev_resseq,
                    "omega_degree": safe_float(omega_deg, 3),
                    "peptide_cn_distance": safe_float(cn_dist, 3),
                    "current_resname": cur_res.get_resname().strip(),
                    "previous_resname": prev_res.get_resname().strip()
                })
                out[chain_id]["cis_count"] += 1

    return out


# ============================================================
# JSON 更新
# ============================================================

def update_chain_json_local_features(
    chain_json_path: str,
    dssp_result: Optional[Dict],
    rsa_result: Optional[Dict],
    hydropathy_values: List[float],
    disorder_values: Optional[List[float]],
    clash_result: Optional[Dict],
    bond_result: Optional[Dict],
    cis_result: Optional[Dict],
    args
) -> str:
    data = read_json(chain_json_path)

    computed = data.setdefault("computed_features", {})

    computed["hydropathy"] = {
        "status": "done",
        "type": "values",
        "source": {
            "method": "Kyte-Doolittle",
            "window": args.hydro_window
        },
        "coordinate_system": "chain_index_1_based",
        "values": hydropathy_values
    }

    if dssp_result is not None:
        computed["secondary_structure"] = {
            "status": "done",
            "type": "segments",
            "source": {
                "method": "DSSP",
                "executable": args.dssp_exe
            },
            "coordinate_system": "chain_index_1_based",
            "per_position": dssp_result.get("per_position", []),
            "items": dssp_result.get("segments", [])
        }
    else:
        computed["secondary_structure"] = {
            "status": "failed",
            "type": "segments",
            "source": {"method": "DSSP"},
            "coordinate_system": "chain_index_1_based",
            "per_position": [],
            "items": []
        }

    if rsa_result is not None:
        computed["rsa"] = {
            "status": "done",
            "type": "values",
            "source": {
                "method": "FreeSASA"
            },
            "coordinate_system": "chain_index_1_based",
            "values": rsa_result.get("values", []),
            "mapped_count": rsa_result.get("mapped_count", 0)
        }

        computed["buried_residues"] = {
            "status": "done",
            "type": "points",
            "source": {
                "method": "FreeSASA",
                "threshold": args.buried_threshold
            },
            "coordinate_system": "chain_index_1_based",
            "positions": rsa_result.get("buried_positions", [])
        }
    else:
        computed["rsa"] = {
            "status": "failed",
            "type": "values",
            "source": {"method": "FreeSASA"},
            "coordinate_system": "chain_index_1_based",
            "values": []
        }

        computed["buried_residues"] = {
            "status": "failed",
            "type": "points",
            "source": {"method": "FreeSASA"},
            "coordinate_system": "chain_index_1_based",
            "positions": []
        }

    if disorder_values is not None:
        computed["disorder"] = {
            "status": "done",
            "type": "values",
            "source": {
                "method": "IUPred2A",
                "mode": args.iupred_mode
            },
            "coordinate_system": "chain_index_1_based",
            "values": disorder_values
        }
    else:
        computed["disorder"] = {
            "status": "skipped",
            "type": "values",
            "source": {
                "method": "IUPred2A"
            },
            "coordinate_system": "chain_index_1_based",
            "reason": "IUPred2A command not provided or prediction failed.",
            "values": []
        }

    if clash_result is not None:
        computed["atomic_clashes"] = {
            "status": "done",
            "type": "points",
            "source": {
                "method": "internal_vdw_overlap_approximation",
                "overlap_cutoff": args.clash_overlap_cutoff
            },
            "note": "Approximate clash detection; not official MolProbity clashscore.",
            "coordinate_system": "chain_index_1_based",
            "positions": clash_result.get("positions", []),
            "items": clash_result.get("items", []),
            "clash_pair_count": clash_result.get("clash_pair_count", 0)
        }
    else:
        computed["atomic_clashes"] = {
            "status": "failed",
            "type": "points",
            "coordinate_system": "chain_index_1_based",
            "positions": [],
            "items": []
        }

    if bond_result is not None:
        computed["bond_outliers"] = {
            "status": "done",
            "type": "points",
            "source": {
                "method": "internal_backbone_bond_length_approximation",
                "tolerance": args.bond_tolerance
            },
            "note": "Approximate bond-length outlier detection; not official wwPDB/Phenix validation.",
            "coordinate_system": "chain_index_1_based",
            "positions": bond_result.get("positions", []),
            "items": bond_result.get("items", []),
            "outlier_count": bond_result.get("outlier_count", 0)
        }
    else:
        computed["bond_outliers"] = {
            "status": "failed",
            "type": "points",
            "coordinate_system": "chain_index_1_based",
            "positions": [],
            "items": []
        }

    if cis_result is not None:
        computed["cis_peptide"] = {
            "status": "done",
            "type": "points",
            "source": {
                "method": "internal_omega_dihedral",
                "cis_cutoff_degree": args.cis_cutoff_degree,
                "peptide_bond_max_distance": args.peptide_bond_max_distance
            },
            "coordinate_system": "chain_index_1_based",
            "positions": cis_result.get("positions", []),
            "items": cis_result.get("items", []),
            "cis_count": cis_result.get("cis_count", 0)
        }
    else:
        computed["cis_peptide"] = {
            "status": "failed",
            "type": "points",
            "coordinate_system": "chain_index_1_based",
            "positions": [],
            "items": []
        }

    statuses = [
        computed["hydropathy"]["status"],
        computed["secondary_structure"]["status"],
        computed["rsa"]["status"],
        computed["buried_residues"]["status"],
        computed["disorder"]["status"],
        computed["atomic_clashes"]["status"],
        computed["bond_outliers"]["status"],
        computed["cis_peptide"]["status"],
    ]

    if all(s in {"done", "skipped"} for s in statuses):
        final_status = "done"
    elif any(s == "done" for s in statuses):
        final_status = "partial"
    else:
        final_status = "failed"

    data.setdefault("status", {})
    data["status"]["local_tracks"] = final_status

    append_log(
        data,
        step="local_tracks",
        status=final_status,
        message="Local sequence/structure tracks computed and written."
    )

    write_json(chain_json_path, data)

    return final_status


# ============================================================
# Entry 处理
# ============================================================

def find_chain_jsons(json_root: str, entry_id: str) -> Dict[str, str]:
    chain_dir = os.path.join(json_root, entry_id, "chains")

    if not os.path.isdir(chain_dir):
        return {}

    out = {}

    for p in sorted(Path(chain_dir).glob("*.json")):
        try:
            data = read_json(str(p))
            chain_id = data["chain"]["chain_id"]
            out[chain_id] = str(p)
        except Exception:
            continue

    return out


def load_chain_jsons(chain_paths: Dict[str, str]) -> Dict[str, Dict]:
    return {
        chain_id: read_json(path)
        for chain_id, path in chain_paths.items()
    }


def update_entry_json_status(json_root: str, entry_id: str, status: str) -> None:
    entry_path = os.path.join(json_root, entry_id, "entry.json")

    if not os.path.exists(entry_path):
        return

    try:
        data = read_json(entry_path)
        data.setdefault("status", {})
        data["status"]["local_tracks"] = status

        append_log(
            data,
            step="local_tracks",
            status=status,
            message="Entry-level local track status updated."
        )

        write_json(entry_path, data)

    except Exception:
        return


def process_one_entry(row: Dict, json_root: str, args) -> Tuple[str, Dict]:
    updated = dict(row)

    entry_id = row.get("entry_id", "").strip()
    structure_path = row.get("structure_path", "").strip()

    if not entry_id:
        updated["local_tracks_status"] = "failed"
        updated["note"] = append_note(updated.get("note", ""), "Missing entry_id")
        return "failed", updated

    if not structure_path or not os.path.exists(structure_path):
        updated["local_tracks_status"] = "failed"
        updated["note"] = append_note(updated.get("note", ""), "Missing structure_path")
        return "failed", updated

    chain_paths = find_chain_jsons(json_root, entry_id)

    if not chain_paths:
        updated["local_tracks_status"] = "failed"
        updated["note"] = append_note(updated.get("note", ""), "No chain JSON found")
        return "failed", updated

    chain_jsons = load_chain_jsons(chain_paths)

    temp_path = None

    try:
        structure, parse_path, temp_path = parse_structure(structure_path, entry_id)

        dssp_all = None
        try:
            dssp_all = compute_dssp_for_entry(
                structure=structure,
                parse_path=parse_path,
                chain_jsons=chain_jsons,
                dssp_exe=args.dssp_exe
            )
        except Exception as e:
            print(f"[WARN] DSSP failed for {entry_id}: {repr(e)}")
            dssp_all = None

        rsa_all = None
        try:
            rsa_all = compute_freesasa_for_entry(
                structure_path=parse_path,
                chain_jsons=chain_jsons,
                buried_threshold=args.buried_threshold
            )
        except Exception as e:
            print(f"[WARN] FreeSASA failed for {entry_id}: {repr(e)}")
            rsa_all = None

        clashes_all = None
        try:
            clashes_all = compute_atomic_clashes_approx(
                structure=structure,
                chain_jsons=chain_jsons,
                overlap_cutoff=args.clash_overlap_cutoff,
                search_radius=args.clash_search_radius
            )
        except Exception as e:
            print(f"[WARN] Atomic clashes failed for {entry_id}: {repr(e)}")
            clashes_all = None

        bonds_all = None
        try:
            bonds_all = compute_bond_outliers_approx(
                structure=structure,
                chain_jsons=chain_jsons,
                tolerance=args.bond_tolerance,
                peptide_bond_max_distance=args.peptide_bond_max_distance
            )
        except Exception as e:
            print(f"[WARN] Bond outliers failed for {entry_id}: {repr(e)}")
            bonds_all = None

        cis_all = None
        try:
            cis_all = compute_cis_peptides(
                structure=structure,
                chain_jsons=chain_jsons,
                cis_cutoff_degree=args.cis_cutoff_degree,
                peptide_bond_max_distance=args.peptide_bond_max_distance
            )
        except Exception as e:
            print(f"[WARN] Cis-peptide failed for {entry_id}: {repr(e)}")
            cis_all = None

        hydro_cache = {}
        disorder_cache = {}

        chain_statuses = []

        for chain_id, chain_path in chain_paths.items():
            cj = chain_jsons[chain_id]
            seq = cj["chain"]["sequence"]
            seq_hash = cj["chain"]["sequence_hash"]

            if seq_hash not in hydro_cache:
                hydro_cache[seq_hash] = compute_hydropathy(
                    seq=seq,
                    window=args.hydro_window
                )

            if seq_hash not in disorder_cache:
                if args.iupred_cmd:
                    disorder_cache[seq_hash] = run_iupred2a(
                        seq=seq,
                        iupred_cmd=args.iupred_cmd,
                        mode=args.iupred_mode
                    )
                else:
                    disorder_cache[seq_hash] = None

            status = update_chain_json_local_features(
                chain_json_path=chain_path,
                dssp_result=dssp_all.get(chain_id) if dssp_all else None,
                rsa_result=rsa_all.get(chain_id) if rsa_all else None,
                hydropathy_values=hydro_cache[seq_hash],
                disorder_values=disorder_cache[seq_hash],
                clash_result=clashes_all.get(chain_id) if clashes_all else None,
                bond_result=bonds_all.get(chain_id) if bonds_all else None,
                cis_result=cis_all.get(chain_id) if cis_all else None,
                args=args
            )

            chain_statuses.append(status)

        cleanup_temp(temp_path)

        if all(s == "done" for s in chain_statuses):
            entry_status = "done"
        elif any(s in {"done", "partial"} for s in chain_statuses):
            entry_status = "partial"
        else:
            entry_status = "failed"

        updated["local_tracks_status"] = entry_status
        updated["note"] = append_note(
            updated.get("note", ""),
            f"Local tracks computed; status={entry_status}"
        )

        update_entry_json_status(json_root, entry_id, entry_status)

        return entry_status, updated

    except Exception as e:
        cleanup_temp(temp_path)

        updated["local_tracks_status"] = "failed"
        updated["note"] = append_note(
            updated.get("note", ""),
            f"Local tracks error: {repr(e)}"
        )

        update_entry_json_status(json_root, entry_id, "failed")

        return "failed", updated


# ============================================================
# Manifest / Batch
# ============================================================

def read_manifest(path: str) -> Tuple[List[Dict], List[str]]:
    with open(path, "r", encoding="utf-8", errors="ignore", newline="") as f:
        reader = csv.DictReader(f)

        if reader.fieldnames is None:
            raise RuntimeError(f"No header found in manifest: {path}")

        return list(reader), list(reader.fieldnames)


def write_manifest(path: str, rows: List[Dict], fieldnames: List[str]) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    out_fields = list(fieldnames)

    for col in ["local_tracks_status", "note"]:
        if col not in out_fields:
            out_fields.append(col)

    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(rows)


def select_batch(rows: List[Dict], batch_size: int, batch_index: int) -> Tuple[List[Dict], int, int]:
    start = batch_index * batch_size
    end = min(start + batch_size, len(rows))
    return rows[start:end], start, end


def run_batch(args) -> None:
    rows, fieldnames = read_manifest(args.manifest)

    eligible_indices = []

    for i, row in enumerate(rows):
        if row.get("has_structure", "").lower() == "yes" and row.get("skeleton_status", "").lower() == "done":
            eligible_indices.append(i)

    eligible_rows = [rows[i] for i in eligible_indices]

    batch_rows, bstart, bend = select_batch(
        eligible_rows,
        batch_size=args.batch_size,
        batch_index=args.batch_index
    )

    print("========== Batch Info ==========")
    print(f"Manifest rows total      : {len(rows)}")
    print(f"Eligible rows            : {len(eligible_rows)}")
    print(f"Batch size               : {args.batch_size}")
    print(f"Batch index              : {args.batch_index}")
    print(f"Eligible slice           : [{bstart}, {bend})")
    print(f"Rows in this batch       : {len(batch_rows)}")
    print("================================")

    stats = {
        "done": 0,
        "partial": 0,
        "failed": 0,
    }

    entry_to_global_index = {}

    for i, row in enumerate(rows):
        entry_to_global_index.setdefault(row.get("entry_id", ""), i)

    for n, row in enumerate(batch_rows, start=1):
        entry_id = row.get("entry_id", "")

        print(f"[INFO] ({n}/{len(batch_rows)}) Processing local tracks: {entry_id}")

        status, updated = process_one_entry(
            row=row,
            json_root=args.json_root,
            args=args
        )

        stats[status] = stats.get(status, 0) + 1

        global_i = entry_to_global_index.get(entry_id)

        if global_i is not None:
            rows[global_i] = updated

    write_manifest(args.manifest_out, rows, fieldnames)

    print("========== Local Tracks Summary ==========")
    print(f"Manifest input      : {args.manifest}")
    print(f"JSON root           : {args.json_root}")
    print(f"Manifest output     : {args.manifest_out}")
    print(f"Processed in batch  : {len(batch_rows)}")
    print(f"Done                : {stats.get('done', 0)}")
    print(f"Partial             : {stats.get('partial', 0)}")
    print(f"Failed              : {stats.get('failed', 0)}")
    print("==========================================")


def main():
    parser = argparse.ArgumentParser(
        description="Compute local sequence/structure tracks and write them into chain-level JSON files."
    )

    parser.add_argument(
        "--manifest",
        required=True,
        help="输入 manifest_skeleton.csv"
    )

    parser.add_argument(
        "--json-root",
        required=True,
        help="JSON 根目录，例如 example_outputs/json"
    )

    parser.add_argument(
        "--manifest-out",
        required=True,
        help="输出更新后的 manifest"
    )

    parser.add_argument(
        "--batch-size",
        type=int,
        default=3000,
        help="每批处理多少个 entry，默认 3000"
    )

    parser.add_argument(
        "--batch-index",
        type=int,
        default=0,
        help="第几批，从 0 开始"
    )

    parser.add_argument(
        "--dssp-exe",
        default="auto",
        help="DSSP 可执行程序，默认 auto，可选 mkdssp/dssp/完整路径"
    )

    parser.add_argument(
        "--buried-threshold",
        type=float,
        default=0.20,
        help="RSA 小于该值视为 buried，默认 0.20"
    )

    parser.add_argument(
        "--hydro-window",
        type=int,
        default=9,
        help="Hydropathy 滑窗大小，默认 9"
    )

    parser.add_argument(
        "--iupred-cmd",
        default="",
        help="IUPred2A 命令，例如 'python /path/to/iupred2a.py'。不提供则 disorder 标记 skipped。"
    )

    parser.add_argument(
        "--iupred-mode",
        default="long",
        choices=["long", "short", "glob"],
        help="IUPred2A 模式，默认 long"
    )

    parser.add_argument(
        "--clash-overlap-cutoff",
        type=float,
        default=0.40,
        help="Atomic clash 近似判断的 VDW overlap 阈值，默认 0.40 Å"
    )

    parser.add_argument(
        "--clash-search-radius",
        type=float,
        default=4.0,
        help="Atomic clash 邻近搜索半径，默认 4.0 Å"
    )

    parser.add_argument(
        "--bond-tolerance",
        type=float,
        default=0.08,
        help="Bond outlier 近似判断的键长偏差阈值，默认 0.08 Å"
    )

    parser.add_argument(
        "--cis-cutoff-degree",
        type=float,
        default=30.0,
        help="Cis-peptide 判定阈值，|omega| 小于该角度认为是 cis，默认 30°"
    )

    parser.add_argument(
        "--peptide-bond-max-distance",
        type=float,
        default=1.9,
        help="判断相邻残基是否存在 peptide bond 的 C-N 最大距离，默认 1.9 Å"
    )

    args = parser.parse_args()

    run_batch(args)


if __name__ == "__main__":
    main()