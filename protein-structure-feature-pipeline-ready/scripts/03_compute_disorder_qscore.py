#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
03_compute_disorder_qscore.py

功能：
在已有 chain-level JSON 的基础上，补充两个 feature：

1. disorder
   - 基于序列
   - 默认使用 IUPred2A
   - 输出每个 residue 的 disorder score
   - 生成 disorder segments

2. qscore
   - 基于 PDB + cryo-EM map
   - 通过 ChimeraX + QScore 命令计算
   - 输出每个 residue 的 Q-score
   - 映射回每条链 JSON 的 chain_index 坐标

重要修改：
- QScore 命令使用：
    qscore #1 toVolume #2 ...
  而不是：
    qscore residues #1 toVolume #2 ...
- 增强了 QScore 输出解析，支持 Qavg 表格格式。
"""

import argparse
import csv
import json
import os
import re
import shlex
import shutil
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ============================================================
# 基础工具
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


def read_manifest(path: str) -> Tuple[List[Dict], List[str]]:
    with open(path, "r", encoding="utf-8", errors="ignore", newline="") as f:
        r = csv.DictReader(f)
        return list(r), list(r.fieldnames or [])


def write_manifest(path: str, rows: List[Dict], fieldnames: List[str]) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    out_fields = list(fieldnames)

    for col in [
        "disorder_qscore_status",
        "disorder_status",
        "qscore_status",
        "note"
    ]:
        if col not in out_fields:
            out_fields.append(col)

    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=out_fields,
            extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(rows)


def append_note(old: str, new: str) -> str:
    old = old or ""
    if not new:
        return old
    if not old:
        return new
    if new in old:
        return old
    return old + "; " + new


def sniff_delimiter(path: str) -> str:
    if path.lower().endswith(".tsv"):
        return "\t"
    return ","


def find_chain_jsons(json_root: str, entry_id: str) -> Dict[str, str]:
    chain_dir = Path(json_root) / entry_id / "chains"
    if not chain_dir.is_dir():
        return {}

    out = {}

    for p in sorted(chain_dir.glob("*.json")):
        try:
            d = read_json(str(p))
            chain_id = str(d["chain"]["chain_id"])
            out[chain_id] = str(p)
        except Exception:
            continue

    return out


def load_chain_jsons(chain_paths: Dict[str, str]) -> Dict[str, Dict]:
    return {
        chain_id: read_json(path)
        for chain_id, path in chain_paths.items()
    }


def update_entry_status(
    json_root: str,
    entry_id: str,
    disorder_status: str,
    qscore_status: str
) -> None:
    entry_path = Path(json_root) / entry_id / "entry.json"
    if not entry_path.exists():
        return

    try:
        d = read_json(str(entry_path))
        d.setdefault("status", {})

        d["status"]["disorder"] = disorder_status
        d["status"]["qscore"] = qscore_status

        if disorder_status in {"done", "skipped"} and qscore_status in {"done", "skipped"}:
            final = "done"
        elif disorder_status == "failed" and qscore_status == "failed":
            final = "failed"
        else:
            final = "partial"

        d["status"]["disorder_qscore"] = final

        append_log(
            d,
            step="disorder_qscore",
            status=final,
            message=f"Disorder status={disorder_status}; Q-score status={qscore_status}."
        )

        write_json(str(entry_path), d)

    except Exception:
        return


# ============================================================
# Disorder: IUPred2A
# ============================================================

def build_segments_from_values(
    values: List[Optional[float]],
    threshold: float = 0.5
) -> List[Dict]:
    """
    把连续 score >= threshold 的区域合并成 segment。
    """

    segs = []
    start = None

    for i, v in enumerate(values, start=1):
        is_disordered = (v is not None and v >= threshold)

        if is_disordered and start is None:
            start = i

        if not is_disordered and start is not None:
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


def run_iupred2a(
    seq: str,
    iupred_cmd: str,
    mode: str = "long",
    timeout: int = 180
) -> Optional[List[float]]:
    """
    运行 IUPred2A。

    典型命令：
      --iupred-cmd "python /path/to/iupred2a.py"

    实际执行：
      python iupred2a.py temp.fasta long

    解析格式：
      pos aa score
    """

    if not iupred_cmd:
        return None

    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta", mode="w")
    tmp.write(">query\n")
    tmp.write(seq + "\n")
    tmp.close()

    try:
        cmd = shlex.split(iupred_cmd) + [tmp.name, mode]

        proc = subprocess.run(
            cmd,
            text=True,
            capture_output=True,
            timeout=timeout
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
        values = [round(float(s), 4) for _, s in scores]

        if len(values) != len(seq):
            return None

        return values

    finally:
        if os.path.exists(tmp.name):
            os.remove(tmp.name)


def write_disorder_to_chain_json(
    chain_path: str,
    values: Optional[List[float]],
    mode: str,
    threshold: float,
    reason: str = ""
) -> str:
    d = read_json(chain_path)
    L = int(d["chain"]["length"])

    computed = d.setdefault("computed_features", {})

    if values is not None and len(values) == L:
        computed["disorder"] = {
            "status": "done",
            "type": "values",
            "source": {
                "method": "IUPred2A",
                "mode": mode,
                "threshold": threshold
            },
            "coordinate_system": "chain_index_1_based",
            "values": values,
            "segments": build_segments_from_values(values, threshold=threshold)
        }
        status = "done"

    else:
        computed["disorder"] = {
            "status": "failed" if reason else "skipped",
            "type": "values",
            "source": {
                "method": "IUPred2A",
                "mode": mode,
                "threshold": threshold
            },
            "coordinate_system": "chain_index_1_based",
            "reason": reason or "IUPred2A command not provided or disorder was skipped.",
            "values": [],
            "segments": []
        }
        status = computed["disorder"]["status"]

    d.setdefault("status", {})
    d["status"]["disorder"] = status

    append_log(
        d,
        step="disorder",
        status=status,
        message=f"Disorder computed with IUPred2A mode={mode}."
    )

    write_json(chain_path, d)
    return status


def get_existing_feature_status(chain_path: str, feature_name: str, default: str = "skipped") -> str:
    try:
        d = read_json(chain_path)
        return d.get("computed_features", {}).get(feature_name, {}).get("status", default)
    except Exception:
        return default


# ============================================================
# Q-score: ChimeraX + map
# ============================================================

def load_map_table(path: str) -> Dict[str, str]:
    """
    读取 map table。

    支持 CSV/TSV，推荐字段：
      entry_id,map_path

    也兼容：
      cfdb_id,map_path
      map / map_file
    """

    if not path:
        return {}

    out = {}
    delimiter = sniff_delimiter(path)

    with open(path, newline="", encoding="utf-8", errors="ignore") as f:
        r = csv.DictReader(f, delimiter=delimiter)

        for row in r:
            key = row.get("entry_id") or row.get("cfdb_id")
            mp = row.get("map_path") or row.get("map") or row.get("map_file")

            if key and mp and os.path.exists(mp):
                out[str(key).strip()] = str(mp).strip()

    return out


def resolve_map_path(row: Dict, args, map_table: Dict[str, str]) -> Optional[str]:
    """
    多种方式寻找 map 路径：

    1. map table
    2. manifest 里的 map_path / map_file / em_map_path
    3. --map-root + --map-pattern
    """

    entry_id = row.get("entry_id", "").strip()
    cfdb_id = row.get("cfdb_id", "").strip()
    emdb_id = row.get("emdb_id", "").strip()
    pdb_id = row.get("pdb_id", "").strip()

    for key in [entry_id, cfdb_id]:
        if key in map_table:
            p = map_table[key]
            if os.path.exists(p):
                return p

    for col in ["map_path", "map_file", "em_map_path", "density_map_path"]:
        p = row.get(col, "")
        if p and os.path.exists(p):
            return p

    if args.map_root and args.map_pattern:
        try:
            fname = args.map_pattern.format(
                entry_id=entry_id,
                cfdb_id=cfdb_id,
                emdb_id=emdb_id,
                pdb_id=pdb_id
            )
            p = os.path.join(args.map_root, fname)
            if os.path.exists(p):
                return p
        except Exception:
            pass

    return None


def find_chimerax(chimerax_cmd: str) -> Optional[str]:
    if chimerax_cmd:
        if os.path.exists(chimerax_cmd):
            return chimerax_cmd

        found = shutil.which(chimerax_cmd)
        if found:
            return found

        return chimerax_cmd

    for name in ["chimerax", "ChimeraX", "UCSF-ChimeraX"]:
        p = shutil.which(name)
        if p:
            return p

    return None


def make_chimerax_qscore_script(
    pdb_path: str,
    map_path: str,
    output_path: str,
    qscore_command: str
) -> str:
    """
    生成 ChimeraX command script。

    注意：
    正确形式是：
      qscore #1 toVolume #2 ...

    不是：
      qscore residues #1 toVolume #2 ...
    """

    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".cxc", mode="w")
    tmp_path = tmp.name

    tmp.write(f'open "{pdb_path}"\n')
    tmp.write(f'open "{map_path}"\n')
    tmp.write(qscore_command.format(output_path=output_path) + "\n")
    tmp.write("exit\n")
    tmp.close()

    return tmp_path


def stdout_has_qscore_table(text: str) -> bool:
    for line in text.splitlines():
        low = line.lower()
        if "chain" in low and "number" in low and "qavg" in low:
            return True
    return False


def run_chimerax_qscore(
    pdb_path: str,
    map_path: str,
    chimerax_cmd: str,
    timeout: int = 1800
) -> Tuple[Optional[str], str]:
    """
    运行 ChimeraX QScore，返回 QScore 输出文件路径。

    会依次尝试多个 QScore 命令写法。
    """

    chimerax = find_chimerax(chimerax_cmd)

    if not chimerax:
        return None, "Cannot find ChimeraX executable."

    command_templates = [
        # 推荐主命令
        'qscore #1 toVolume #2 useGui false logDetails true outputFile "{output_path}"',

        # 某些版本可能不支持 useGui
        'qscore #1 toVolume #2 logDetails true outputFile "{output_path}"',

        # 更简化版本
        'qscore #1 toVolume #2 outputFile "{output_path}"',

        # 最简版本：如果没有 outputFile，则从 stdout 里尝试解析
        'qscore #1 toVolume #2',
    ]

    last_error = ""

    for qcmd in command_templates:
        qscore_out = tempfile.NamedTemporaryFile(delete=False, suffix=".qscore.txt")
        qscore_out_path = qscore_out.name
        qscore_out.close()

        cxc_path = make_chimerax_qscore_script(
            pdb_path=pdb_path,
            map_path=map_path,
            output_path=qscore_out_path,
            qscore_command=qcmd
        )

        try:
            cmd = [chimerax, "--nogui", cxc_path]

            proc = subprocess.run(
                cmd,
                text=True,
                capture_output=True,
                timeout=timeout
            )

            stdout = proc.stdout or ""
            stderr = proc.stderr or ""

            if proc.returncode == 0:
                if os.path.exists(qscore_out_path) and os.path.getsize(qscore_out_path) > 0:
                    return qscore_out_path, ""

                if stdout_has_qscore_table(stdout):
                    with open(qscore_out_path, "w", encoding="utf-8") as f:
                        f.write(stdout)

                    if os.path.getsize(qscore_out_path) > 0:
                        return qscore_out_path, ""

            last_error = (
                f"ChimeraX failed or produced no QScore output. "
                f"qcmd={qcmd}; "
                f"returncode={proc.returncode}; "
                f"stdout_tail={stdout[-2000:]}; "
                f"stderr_tail={stderr[-3000:]}"
            )

        except Exception as e:
            last_error = f"ChimeraX exception with qcmd={qcmd}: {repr(e)}"

        finally:
            if os.path.exists(cxc_path):
                os.remove(cxc_path)

            if os.path.exists(qscore_out_path) and os.path.getsize(qscore_out_path) == 0:
                os.remove(qscore_out_path)

    return None, last_error


def normalize_chain_id(x: str) -> str:
    """
    兼容 ChimeraX 输出中的 chain 表达形式。
    """

    x = str(x or "").strip()

    if not x:
        return ""

    if "/" in x:
        x = x.split("/")[-1]

    if "#" in x:
        x = x.split("#")[-1]

    return x.strip()


def parse_residue_number(text: str) -> Optional[Tuple[int, str]]:
    """
    解析 residue number 和 insertion code。

    支持：
      117
      117A
      A:117
      /A:117
    """

    text = str(text or "").strip()

    if not text:
        return None

    m = re.search(r"(-?\d+)([A-Za-z]?)", text)

    if not m:
        return None

    return int(m.group(1)), m.group(2).strip()


def parse_qscore_output(qscore_path: str) -> List[Dict]:
    """
    解析 ChimeraX QScore 输出。

    优先兼容标准表格：
      Chain Number Name Qavg Qworst Qbb Qsc
      A     1      MET  0.198 ...

    返回：
      [
        {"chain_id": "A", "resseq": 1, "icode": "", "qscore": 0.198}
      ]
    """

    rows = []

    with open(qscore_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = [x.rstrip("\n") for x in f if x.strip()]

    if not lines:
        return rows

    # ------------------------------------------------------------
    # 1. 空白分隔表格：Chain Number Name Qavg ...
    # ------------------------------------------------------------
    header_idx = None
    header_parts = None

    for i, line in enumerate(lines):
        parts = line.strip().split()
        low = [p.lower() for p in parts]

        if "chain" in low and "number" in low and "qavg" in low:
            header_idx = i
            header_parts = low
            break

    if header_idx is not None and header_parts:
        try:
            chain_i = header_parts.index("chain")
            number_i = header_parts.index("number")
            qavg_i = header_parts.index("qavg")
        except Exception:
            chain_i = 0
            number_i = 1
            qavg_i = 3

        for line in lines[header_idx + 1:]:
            s = line.strip()

            if not s:
                continue

            if set(s.replace(" ", "")) <= {"-"}:
                continue

            parts = s.split()

            if len(parts) <= max(chain_i, number_i, qavg_i):
                continue

            chain_id = normalize_chain_id(parts[chain_i])
            parsed = parse_residue_number(parts[number_i])

            if not parsed:
                continue

            resseq, icode = parsed

            try:
                qavg_text = parts[qavg_i]
                if qavg_text.upper() == "N/A":
                    continue
                qavg = float(qavg_text)
            except Exception:
                continue

            rows.append({
                "chain_id": chain_id,
                "resseq": resseq,
                "icode": icode,
                "qscore": round(qavg, 4)
            })

        if rows:
            return rows

    # ------------------------------------------------------------
    # 2. CSV / TSV 格式
    # ------------------------------------------------------------
    header_line = None
    delimiter = None

    for line in lines[:50]:
        if "\t" in line:
            parts = line.split("\t")
            delim = "\t"
        elif "," in line:
            parts = line.split(",")
            delim = ","
        else:
            continue

        low = [p.strip().lower() for p in parts]

        has_chain = any("chain" == x or "chain" in x for x in low)
        has_score = any(
            "qavg" in x
            or "qscore" in x
            or "q-score" in x
            or x == "score"
            for x in low
        )

        if has_chain and has_score:
            header_line = line
            delimiter = delim
            break

    if header_line and delimiter:
        start = lines.index(header_line)
        reader = csv.DictReader(lines[start:], delimiter=delimiter)

        for r in reader:
            keys = {k.lower().strip(): k for k in r.keys() if k is not None}

            chain_key = None
            res_key = None
            score_key = None

            for lk, ok in keys.items():
                if chain_key is None and "chain" in lk:
                    chain_key = ok

                if res_key is None and (
                    lk in {"res", "residue", "residue number", "number", "resnum", "resseq"}
                    or "residue" in lk
                ):
                    res_key = ok

                if score_key is None and (
                    "qavg" in lk
                    or "qscore" in lk
                    or "q-score" in lk
                    or lk == "score"
                ):
                    score_key = ok

            if not chain_key or not res_key or not score_key:
                continue

            try:
                chain_id = normalize_chain_id(r[chain_key])
                parsed = parse_residue_number(r[res_key])
                if not parsed:
                    continue

                resseq, icode = parsed
                q = float(str(r[score_key]).strip())

                rows.append({
                    "chain_id": chain_id,
                    "resseq": resseq,
                    "icode": icode,
                    "qscore": round(q, 4)
                })

            except Exception:
                continue

        if rows:
            return rows

    # ------------------------------------------------------------
    # 3. fallback：逐行正则
    # ------------------------------------------------------------
    line_regex = re.compile(
        r"(?P<chain>[A-Za-z0-9])\s+"
        r"(?P<resnum>-?\d+)(?P<icode>[A-Za-z]?)\s+"
        r"(?P<resname>[A-Za-z]{3})?\s*"
        r"(?P<score>[+-]?\d+\.\d+)"
    )

    for line in lines:
        m = line_regex.search(line)

        if not m:
            continue

        try:
            rows.append({
                "chain_id": normalize_chain_id(m.group("chain")),
                "resseq": int(m.group("resnum")),
                "icode": (m.group("icode") or "").strip(),
                "qscore": round(float(m.group("score")), 4)
            })
        except Exception:
            continue

    return rows


def build_residue_maps(chain_jsons: Dict[str, Dict]) -> Dict[str, Dict[Tuple[int, str], int]]:
    """
    chain_id -> (pdb_resseq, pdb_icode) -> chain_index
    """

    out = {}

    for chain_id, d in chain_jsons.items():
        m = {}

        for item in d.get("residue_mapping", []):
            try:
                resseq = int(item.get("pdb_resseq"))
                icode = str(item.get("pdb_icode", "")).strip()
                idx = int(item.get("chain_index"))
                m[(resseq, icode)] = idx
            except Exception:
                continue

        out[str(chain_id)] = m

    return out


def map_qscores_to_chains(
    qscore_rows: List[Dict],
    chain_jsons: Dict[str, Dict]
) -> Dict[str, List[Optional[float]]]:
    """
    把 qscore 输出映射到每条链的 chain_index。
    """

    maps = build_residue_maps(chain_jsons)

    out = {}

    for chain_id, d in chain_jsons.items():
        L = int(d["chain"]["length"])
        out[str(chain_id)] = [None] * L

    for row in qscore_rows:
        chain_id = normalize_chain_id(row.get("chain_id", ""))

        if chain_id not in out:
            continue

        try:
            resseq = int(row["resseq"])
            icode = str(row.get("icode", "")).strip()
            q = float(row["qscore"])
        except Exception:
            continue

        idx = maps.get(chain_id, {}).get((resseq, icode))

        if idx is None:
            # 放宽：只按 resseq 匹配
            for (rseq, ic), ci in maps.get(chain_id, {}).items():
                if rseq == resseq:
                    idx = ci
                    break

        if idx is None:
            continue

        if 1 <= idx <= len(out[chain_id]):
            out[chain_id][idx - 1] = round(q, 4)

    return out


def write_qscore_to_chain_json(
    chain_path: str,
    values: Optional[List[Optional[float]]],
    reason: str = ""
) -> str:
    d = read_json(chain_path)
    L = int(d["chain"]["length"])

    computed = d.setdefault("computed_features", {})

    if values is not None and len(values) == L and any(v is not None for v in values):
        valid = [v for v in values if v is not None]

        computed["qscore"] = {
            "status": "done",
            "type": "values",
            "source": {
                "method": "ChimeraX QScore"
            },
            "coordinate_system": "chain_index_1_based",
            "values": values,
            "mapped_count": len(valid),
            "mean": round(sum(valid) / len(valid), 4) if valid else None
        }

        status = "done"

    else:
        computed["qscore"] = {
            "status": "failed" if reason else "skipped",
            "type": "values",
            "source": {
                "method": "ChimeraX QScore"
            },
            "coordinate_system": "chain_index_1_based",
            "reason": reason or "Q-score not requested or missing map/ChimeraX.",
            "values": [],
            "mapped_count": 0,
            "mean": None
        }

        status = computed["qscore"]["status"]

    d.setdefault("status", {})
    d["status"]["qscore"] = status

    append_log(
        d,
        step="qscore",
        status=status,
        message="Q-score written to chain JSON."
    )

    write_json(chain_path, d)
    return status


# ============================================================
# 主处理逻辑
# ============================================================

def process_one_entry(row: Dict, json_root: str, map_table: Dict[str, str], args) -> Dict:
    updated = dict(row)

    entry_id = row.get("entry_id", "").strip()
    structure_path = row.get("structure_path", "").strip()

    chain_paths = find_chain_jsons(json_root, entry_id)

    if not chain_paths:
        updated["disorder_qscore_status"] = "failed"
        updated["disorder_status"] = "failed" if not args.skip_disorder else "skipped"
        updated["qscore_status"] = "failed" if not args.skip_qscore else "skipped"
        updated["note"] = append_note(updated.get("note", ""), "No chain JSON found")
        return updated

    chain_jsons = load_chain_jsons(chain_paths)

    # -------------------------
    # Disorder
    # -------------------------
    disorder_statuses = []

    if args.skip_disorder:
        for chain_id, chain_path in chain_paths.items():
            st = get_existing_feature_status(chain_path, "disorder", default="skipped")
            disorder_statuses.append(st)

    else:
        seq_cache = {}

        for chain_id, chain_path in chain_paths.items():
            cj = chain_jsons[chain_id]
            seq = cj["chain"]["sequence"]
            seq_hash = cj["chain"].get("sequence_hash", seq)

            if seq_hash not in seq_cache:
                values = run_iupred2a(
                    seq=seq,
                    iupred_cmd=args.iupred_cmd,
                    mode=args.iupred_mode,
                    timeout=args.iupred_timeout
                )
                seq_cache[seq_hash] = values

            values = seq_cache[seq_hash]

            if values is None:
                reason = "IUPred2A failed or command not provided." if args.iupred_cmd else ""
            else:
                reason = ""

            st = write_disorder_to_chain_json(
                chain_path=chain_path,
                values=values,
                mode=args.iupred_mode,
                threshold=args.disorder_threshold,
                reason=reason
            )

            disorder_statuses.append(st)

    if all(s == "done" for s in disorder_statuses):
        disorder_status = "done"
    elif all(s == "skipped" for s in disorder_statuses):
        disorder_status = "skipped"
    elif any(s == "done" for s in disorder_statuses):
        disorder_status = "partial"
    else:
        disorder_status = "failed"

    # -------------------------
    # Q-score
    # -------------------------
    qscore_statuses = []

    if args.skip_qscore:
        for chain_id, chain_path in chain_paths.items():
            st = get_existing_feature_status(chain_path, "qscore", default="skipped")
            qscore_statuses.append(st)

    else:
        map_path = resolve_map_path(row, args, map_table)
        qscore_values_by_chain = None
        reason = ""

        if not structure_path or not os.path.exists(structure_path):
            reason = "Missing structure_path."

        elif not map_path or not os.path.exists(map_path):
            reason = "Missing map path."

        else:
            qscore_out, reason = run_chimerax_qscore(
                pdb_path=structure_path,
                map_path=map_path,
                chimerax_cmd=args.chimerax_cmd,
                timeout=args.qscore_timeout
            )

            if qscore_out:
                try:
                    qrows = parse_qscore_output(qscore_out)

                    if not qrows:
                        reason = "QScore output parsed but no residue rows were found."
                        qscore_values_by_chain = None

                    else:
                        qscore_values_by_chain = map_qscores_to_chains(qrows, chain_jsons)
                        reason = ""

                        if args.keep_qscore_raw:
                            raw_dir = Path(args.keep_qscore_raw)
                            raw_dir.mkdir(parents=True, exist_ok=True)
                            raw_copy = raw_dir / f"{entry_id}.qscore.raw.txt"
                            shutil.copyfile(qscore_out, raw_copy)

                except Exception as e:
                    reason = f"Failed to parse/map QScore output: {repr(e)}"
                    qscore_values_by_chain = None

                finally:
                    if os.path.exists(qscore_out):
                        os.remove(qscore_out)

        for chain_id, chain_path in chain_paths.items():
            values = None

            if qscore_values_by_chain and chain_id in qscore_values_by_chain:
                values = qscore_values_by_chain[chain_id]

            st = write_qscore_to_chain_json(
                chain_path=chain_path,
                values=values,
                reason=reason if values is None else ""
            )

            qscore_statuses.append(st)

    if all(s == "done" for s in qscore_statuses):
        qscore_status = "done"
    elif all(s == "skipped" for s in qscore_statuses):
        qscore_status = "skipped"
    elif any(s == "done" for s in qscore_statuses):
        qscore_status = "partial"
    else:
        qscore_status = "failed"

    # -------------------------
    # entry status
    # -------------------------
    update_entry_status(json_root, entry_id, disorder_status, qscore_status)

    updated["disorder_status"] = disorder_status
    updated["qscore_status"] = qscore_status

    if disorder_status in {"done", "skipped"} and qscore_status in {"done", "skipped"}:
        updated["disorder_qscore_status"] = "done"
    elif disorder_status == "failed" and qscore_status == "failed":
        updated["disorder_qscore_status"] = "failed"
    else:
        updated["disorder_qscore_status"] = "partial"

    return updated


def select_batch(
    rows: List[Dict],
    batch_size: int,
    batch_index: int,
    limit: int = 0
) -> Tuple[List[Dict], int, int]:
    start = batch_index * batch_size
    end = min(start + batch_size, len(rows))

    selected = rows[start:end]

    if limit and limit > 0:
        selected = selected[:limit]

    return selected, start, end


def run(args) -> None:
    rows, fieldnames = read_manifest(args.manifest)

    eligible = []

    for row in rows:
        if (
            row.get("has_structure", "").lower() == "yes"
            and row.get("skeleton_status", "").lower() == "done"
        ):
            eligible.append(row)

    selected, start, end = select_batch(
        eligible,
        batch_size=args.batch_size,
        batch_index=args.batch_index,
        limit=args.limit
    )

    map_table = load_map_table(args.map_table)

    print("========== Disorder + Q-score Batch Info ==========")
    print(f"Manifest rows total : {len(rows)}")
    print(f"Eligible rows       : {len(eligible)}")
    print(f"Batch size          : {args.batch_size}")
    print(f"Batch index         : {args.batch_index}")
    print(f"Eligible slice      : [{start}, {end})")
    print(f"Limit               : {args.limit}")
    print(f"Rows selected       : {len(selected)}")
    print(f"skip_disorder       : {args.skip_disorder}")
    print(f"skip_qscore         : {args.skip_qscore}")
    print(f"map table entries   : {len(map_table)}")
    print("===================================================")

    entry_to_global = {}

    for i, row in enumerate(rows):
        entry_to_global.setdefault(row.get("entry_id", ""), i)

    stats = {
        "done": 0,
        "partial": 0,
        "failed": 0
    }

    for n, row in enumerate(selected, start=1):
        entry_id = row.get("entry_id", "")
        print(f"[INFO] ({n}/{len(selected)}) Processing disorder/qscore: {entry_id}")

        updated = process_one_entry(
            row=row,
            json_root=args.json_root,
            map_table=map_table,
            args=args
        )

        st = updated.get("disorder_qscore_status", "failed")
        stats[st] = stats.get(st, 0) + 1

        gi = entry_to_global.get(entry_id)

        if gi is not None:
            rows[gi] = updated

    write_manifest(args.manifest_out, rows, fieldnames)

    print("========== Disorder + Q-score Summary ==========")
    print(f"Manifest input   : {args.manifest}")
    print(f"Manifest output  : {args.manifest_out}")
    print(f"Processed        : {len(selected)}")
    print(f"Done             : {stats.get('done', 0)}")
    print(f"Partial          : {stats.get('partial', 0)}")
    print(f"Failed           : {stats.get('failed', 0)}")
    print("================================================")


def main():
    p = argparse.ArgumentParser(
        description="Compute disorder and Q-score and write them into chain-level JSON files."
    )

    p.add_argument("--manifest", required=True)
    p.add_argument("--json-root", required=True)
    p.add_argument("--manifest-out", required=True)

    p.add_argument("--batch-size", type=int, default=3000)
    p.add_argument("--batch-index", type=int, default=0)
    p.add_argument(
        "--limit",
        type=int,
        default=0,
        help="只处理当前 batch 中前 N 个 entry，用于测试。"
    )

    # disorder
    p.add_argument("--skip-disorder", action="store_true")
    p.add_argument(
        "--iupred-cmd",
        default="",
        help='例如 "python /path/to/iupred2a.py"'
    )
    p.add_argument("--iupred-mode", default="long", choices=["long", "short", "glob"])
    p.add_argument("--iupred-timeout", type=int, default=180)
    p.add_argument("--disorder-threshold", type=float, default=0.5)

    # qscore
    p.add_argument("--skip-qscore", action="store_true")
    p.add_argument(
        "--chimerax-cmd",
        default="",
        help="ChimeraX 可执行程序路径，例如 chimerax。空则自动搜索。"
    )
    p.add_argument("--qscore-timeout", type=int, default=1800)

    p.add_argument(
        "--map-table",
        default="",
        help="可选 CSV/TSV，包含 entry_id,map_path。"
    )
    p.add_argument("--map-root", default="", help="map 根目录。")
    p.add_argument(
        "--map-pattern",
        default="",
        help='例如 "{entry_id}.map" 或 "emd_{emdb_id}.map"。'
    )
    p.add_argument(
        "--keep-qscore-raw",
        default="",
        help="可选，保存 QScore 原始输出目录。"
    )

    args = p.parse_args()
    run(args)


if __name__ == "__main__":
    main()