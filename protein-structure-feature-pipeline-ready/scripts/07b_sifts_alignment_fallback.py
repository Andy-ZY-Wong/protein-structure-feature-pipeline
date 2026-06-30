#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
07b_sifts_alignment_fallback.py

用途：
在已有 chain JSON 上重建 UniProt mapping 与外部注释。

策略：
1. 先尝试 direct SIFTS residue mapping。
   如果 coverage 和 identity 同时达标，直接使用。
2. 如果 direct SIFTS 不可靠，不再依赖本地 pdb_resseq。
   改用当前 JSON 的 chain.sequence 与候选 UniProt sequence 做 alignment fallback。
3. 用最终得到的 chain_index -> uniprot_pos 映射，把 UniProt features、SIFTS Pfam、SIFTS InterPro 转回 chain_index。

注意：
- 不重新计算本地结构指标。
- 会覆盖 JSON 中已有的 uniprot_mapping / external_annotations / tracks.pfam / tracks.interpro / tracks.uniprot_features。
- 会保留 direct SIFTS 的 QC 结果到 sifts_residue_mapping_qc。
"""

import argparse
import csv
import importlib.util
import json
import os
import re
import time
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


DEFAULT_FEATURE_TYPES = {
    "Initiator methionine", "Signal peptide", "Transit peptide", "Propeptide", "Chain", "Peptide",
    "Domain", "Repeat", "Region", "Coiled coil", "Compositional bias", "Motif", "Zinc finger",
    "DNA binding", "Nucleotide binding", "Calcium binding", "Topological domain", "Transmembrane",
    "Intramembrane", "Active site", "Binding site", "Metal binding", "Site", "Modified residue",
    "Lipidation", "Glycosylation", "Disulfide bond", "Cross-link",
}


def now_iso():
    return datetime.now().isoformat(timespec="seconds")


def norm_seq(s):
    return re.sub(r"[^A-Za-z]", "", str(s or "")).upper().replace("*", "")


def read_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def to_builtin_for_json(obj):
    """
    Convert numpy scalars / arrays / sets / Path objects into JSON-serializable
    Python built-in types.
    """
    if obj is None or isinstance(obj, (str, int, float, bool)):
        return obj

    if isinstance(obj, dict):
        return {str(to_builtin_for_json(k)): to_builtin_for_json(v) for k, v in obj.items()}

    if isinstance(obj, (list, tuple)):
        return [to_builtin_for_json(x) for x in obj]

    if isinstance(obj, set):
        return [to_builtin_for_json(x) for x in sorted(obj, key=lambda x: str(x))]

    if isinstance(obj, Path):
        return str(obj)

    # numpy scalar, e.g. numpy.int64 / numpy.float64
    if hasattr(obj, "item"):
        try:
            return to_builtin_for_json(obj.item())
        except Exception:
            pass

    # numpy array
    if hasattr(obj, "tolist"):
        try:
            return to_builtin_for_json(obj.tolist())
        except Exception:
            pass

    return str(obj)


def write_json_atomic(path, data):
    tmp = str(path) + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(to_builtin_for_json(data), f, ensure_ascii=False, indent=2)
    os.replace(tmp, path)


def append_log(d, step, status, message):
    d.setdefault("logs", [])
    d["logs"].append({"time": now_iso(), "step": step, "status": status, "message": message})


def load_base_module(path):
    spec = importlib.util.spec_from_file_location("base_sifts", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def download_url(url, out_path, retries=3, sleep_sec=1.0):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists() and out_path.stat().st_size > 0:
        return
    last = ""
    for _ in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "sequence-remap-fallback/1.0"})
            with urllib.request.urlopen(req, timeout=60) as resp:
                data = resp.read()
            with open(out_path, "wb") as f:
                f.write(data)
            return
        except Exception as e:
            last = repr(e)
            time.sleep(sleep_sec)
    raise RuntimeError(f"Failed to download {url}: {last}")


def fetch_uniprot_json(acc, cache_dir, base=None):
    if base is not None and hasattr(base, "fetch_uniprot_json"):
        p = base.fetch_uniprot_json(acc, Path(cache_dir))
        return read_json(Path(p))

    acc = str(acc or "").strip()
    out = Path(cache_dir) / "uniprot_json" / f"{acc}.json"
    if not out.exists() or out.stat().st_size == 0:
        download_url(f"https://rest.uniprot.org/uniprotkb/{acc}.json", out)
    return read_json(out)


def get_uniprot_sequence(acc, cache_dir, base=None):
    d = fetch_uniprot_json(acc, cache_dir, base=base)
    seq = d.get("sequence", {})
    if isinstance(seq, dict):
        return norm_seq(seq.get("value", ""))
    return ""


def get_entries_from_manifest(path, limit_entries=0):
    rows = []
    with open(path, newline="", encoding="utf-8", errors="ignore") as f:
        rows = list(csv.DictReader(f))
    out, seen = [], set()
    for r in rows:
        eid = r.get("entry_id", "").strip()
        if not eid or eid in seen:
            continue
        has_structure = r.get("has_structure", "").lower()
        skeleton_status = r.get("skeleton_status", "").lower()
        if has_structure and has_structure not in {"yes", "true", "1"}:
            continue
        if skeleton_status and skeleton_status != "done":
            continue
        out.append(eid)
        seen.add(eid)
        if limit_entries and len(out) >= limit_entries:
            break
    return out


def get_entries_from_json_root(json_root, limit_entries=0):
    out = []
    for p in sorted(Path(json_root).iterdir()):
        if p.is_dir() and (p / "chains").is_dir():
            out.append(p.name)
        if limit_entries and len(out) >= limit_entries:
            break
    return out


def find_chain_jsons(json_root, entry_id):
    p = Path(json_root) / entry_id / "chains"
    return sorted(p.glob("*.json")) if p.is_dir() else []



def parse_sifts_xref_rows_standalone(xml_gz_path):
    """
    Standalone parser for SIFTS Pfam / InterPro crossRefDb records.
    This is used when the base script does not return xref_rows.
    """
    import gzip
    import xml.etree.ElementTree as ET

    def strip_ns_local(tag):
        return tag.split("}", 1)[-1] if "}" in tag else tag

    def attr(elem, names, default=""):
        for n in names:
            if n in elem.attrib:
                return elem.attrib.get(n, default)
        lower = {k.lower(): v for k, v in elem.attrib.items()}
        for n in names:
            if n.lower() in lower:
                return lower[n.lower()]
        return default

    def parse_num(x):
        s = str(x or "").strip()
        if not s or s.lower() in {".", "?", "none", "null", "nan"}:
            return None
        m = re.match(r"^(-?\d+)", s)
        if not m:
            return None
        return int(m.group(1))

    rows = []

    with gzip.open(xml_gz_path, "rb") as f:
        root = ET.parse(f).getroot()

    for elem in root.iter():
        if strip_ns_local(elem.tag) != "residue":
            continue

        uniprot_refs = []
        xrefs = []

        for child in list(elem):
            if strip_ns_local(child.tag) != "crossRefDb":
                continue

            src = attr(child, ["dbSource"]).strip()
            coord_sys = attr(child, ["dbCoordSys", "dbCoordSystem"])
            accession = attr(child, ["dbAccessionId", "dbAccessionID", "dbAccession"]).strip()
            resnum = attr(child, ["dbResNum", "dbResnum", "db_res_num"]).strip()
            evidence = attr(child, ["dbEvidence", "evidence"]).strip()

            if src == "UniProt":
                upos = parse_num(resnum)
                if accession and upos is not None:
                    uniprot_refs.append({
                        "uniprot_acc": accession,
                        "uniprot_pos": upos,
                    })

            elif src in {"Pfam", "InterPro"}:
                xrefs.append({
                    "database": src,
                    "accession": accession,
                    "coord_sys": coord_sys,
                    "resnum": resnum,
                    "evidence": evidence,
                    "raw": dict(child.attrib),
                })

        if not xrefs:
            continue

        context_acc = uniprot_refs[0]["uniprot_acc"] if uniprot_refs else ""
        context_upos = uniprot_refs[0]["uniprot_pos"] if uniprot_refs else None

        for x in xrefs:
            coord = str(x.get("coord_sys", "")).lower()
            upos = None

            if coord.startswith("uniprot"):
                upos = parse_num(x.get("resnum"))
            elif context_upos is not None:
                upos = context_upos

            if upos is None:
                continue

            rows.append({
                "database": x["database"],
                "accession": x["accession"],
                "coord_sys": x["coord_sys"],
                "evidence": x["evidence"],
                "uniprot_acc": context_acc,
                "uniprot_pos": upos,
                "pdb_chain": "",
                "pdb_resseq": None,
                "pdb_icode": "",
                "raw": x["raw"],
            })

    return rows

def load_sifts_index(pdb_id, cache_dir, base, cache):
    pdb_id = str(pdb_id or "").strip().lower()
    if pdb_id in cache:
        return cache[pdb_id]

    xml_gz = base.fetch_sifts_xml_gz(pdb_id, Path(cache_dir))
    parsed = base.parse_sifts_xml(Path(xml_gz))

    if isinstance(parsed, dict):
        mapping_rows = parsed.get("mapping_rows", [])
        xref_rows = parsed.get("xref_rows", [])
    else:
        mapping_rows = parsed
        xref_rows = []

    if not xref_rows:
        try:
            xref_rows = parse_sifts_xref_rows_standalone(Path(xml_gz))
        except Exception as e:
            print(f"[WARN] standalone SIFTS xref parse failed for {pdb_id}: {repr(e)}")
            xref_rows = []

    try:
        idx = base.build_sifts_index(mapping_rows, xref_rows)
    except TypeError:
        idx = base.build_sifts_index(mapping_rows)
        idx["xref_rows"] = xref_rows

    idx.setdefault("xref_rows", xref_rows)
    idx["uniprot_accessions"] = sorted(set(r.get("uniprot_acc") for r in mapping_rows if r.get("uniprot_acc")))
    cache[pdb_id] = idx
    return idx


def quality_status(cov, ident, high_cov, high_ident, ok_cov, ok_ident, partial_cov):
    if cov >= high_cov and ident >= high_ident:
        return "HIGH", "done"
    if cov >= ok_cov and ident >= ok_ident:
        return "ACCEPTABLE", "done"
    if cov >= partial_cov and ident >= ok_ident:
        return "PARTIAL_USABLE", "partial"
    return "PROBLEM", "problem"


def score_one_sifts_chain(d, idx, sifts_chain):
    chain_seq = norm_seq(d.get("chain", {}).get("sequence", ""))
    chain_len = int(d.get("chain", {}).get("length", len(chain_seq)))
    mapping = []
    matched = 0

    for item in d.get("residue_mapping", []):
        try:
            cidx = int(item.get("chain_index"))
            resseq = int(item.get("pdb_resseq"))
        except Exception:
            continue
        icode = str(item.get("pdb_icode", "") or "").strip()
        aa = str(item.get("one_letter", "") or "").upper()

        srow = idx["by_chain_res"].get((sifts_chain, resseq, icode))
        if not srow:
            srow = idx["by_chain_res_noicode"].get((sifts_chain, resseq))
        if not srow:
            continue

        uni_aa = srow.get("uniprot_aa", "")
        ok = bool(aa and uni_aa and aa == uni_aa)
        if ok:
            matched += 1
        mapping.append({
            "chain_index": cidx,
            "pdb_chain": sifts_chain,
            "pdb_resseq": resseq,
            "pdb_icode": icode,
            "chain_aa": aa,
            "uniprot_acc": srow.get("uniprot_acc"),
            "uniprot_pos": srow.get("uniprot_pos"),
            "uniprot_aa": uni_aa,
            "aa_match": ok,
        })

    mapping = sorted(mapping, key=lambda x: x["chain_index"])
    mapped = len(mapping)
    cov = mapped / chain_len if chain_len else 0.0
    ident = matched / mapped if mapped else 0.0
    return mapped, matched, cov, ident, mapping


def best_direct_sifts_mapping(d, idx, thresholds):
    high_cov, high_ident, ok_cov, ok_ident, partial_cov = thresholds
    chain_seq = norm_seq(d.get("chain", {}).get("sequence", ""))
    chain_len = int(d.get("chain", {}).get("length", len(chain_seq)))
    local_chain = str(d.get("chain", {}).get("chain_id", ""))

    candidates = []
    for ch in sorted(idx.get("chains", [])):
        mapped, matched, cov, ident, mapping = score_one_sifts_chain(d, idx, ch)
        candidates.append({
            "sifts_chain_id": ch,
            "mapped_count": mapped,
            "matched_count": matched,
            "coverage": cov,
            "identity": ident,
            "mapping": mapping,
        })
    candidates.sort(key=lambda x: (x["identity"], x["coverage"], x["mapped_count"]), reverse=True)

    best = candidates[0] if candidates else {"sifts_chain_id": "", "mapped_count": 0, "coverage": 0, "identity": 0, "mapping": []}
    q, st = quality_status(best["coverage"], best["identity"], high_cov, high_ident, ok_cov, ok_ident, partial_cov)
    if best["mapped_count"] == 0:
        q, st = "FAILED", "failed"

    accs = sorted(set(x.get("uniprot_acc") for x in best.get("mapping", []) if x.get("uniprot_acc")))
    mapped_idx = set(x["chain_index"] for x in best.get("mapping", []))
    unmapped = sorted(set(range(1, chain_len + 1)) - mapped_idx)

    reason = []
    if best["coverage"] < ok_cov:
        reason.append(f"Low SIFTS coverage: {best['coverage']:.4f}")
    if best["identity"] < ok_ident:
        reason.append(f"Low SIFTS identity: {best['identity']:.4f}")
    if best.get("sifts_chain_id") and best.get("sifts_chain_id") != local_chain:
        reason.append(f"Best SIFTS chain differs: {local_chain}->{best.get('sifts_chain_id')}")

    return {
        "status": st,
        "quality": q,
        "source": "SIFTS",
        "method": "direct_sifts_residue_mapping_best_chain",
        "coordinate_system": "chain_index_1_based",
        "chain_id": local_chain,
        "sifts_chain_id": best.get("sifts_chain_id", ""),
        "mapped_count": best.get("mapped_count", 0),
        "chain_length": chain_len,
        "coverage": round(best.get("coverage", 0), 4),
        "identity": round(best.get("identity", 0), 4),
        "uniprot_accessions": accs,
        "unmapped_chain_index": unmapped,
        "mapping": best.get("mapping", []),
        "candidate_scores": [
            {"sifts_chain_id": c["sifts_chain_id"], "mapped_count": c["mapped_count"], "coverage": round(c["coverage"], 4), "identity": round(c["identity"], 4)}
            for c in candidates[:20]
        ],
        "reason": "; ".join(reason),
    }


def exact_subsequence_mapping(chain_seq, uni_seq):
    start = uni_seq.find(chain_seq)
    if start < 0:
        return None
    return [(i, start + i, True) for i in range(1, len(chain_seq) + 1)]


def pairwise_alignment_mapping(chain_seq, uni_seq):
    try:
        from Bio.Align import PairwiseAligner
    except Exception as e:
        raise RuntimeError("Biopython is required. Run: conda install -c conda-forge biopython") from e

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    alns = aligner.align(chain_seq, uni_seq)
    if len(alns) == 0:
        return [], 0.0

    aln = alns[0]
    mapping = []
    for (c0, c1), (u0, u1) in zip(aln.aligned[0], aln.aligned[1]):
        n = min(c1 - c0, u1 - u0)
        for k in range(n):
            ci = c0 + k
            ui = u0 + k
            if ci < len(chain_seq) and ui < len(uni_seq):
                mapping.append((ci + 1, ui + 1, chain_seq[ci] == uni_seq[ui]))
    return sorted(set(mapping), key=lambda x: x[0]), float(aln.score)


def alignment_mapping_object(d, acc, uni_seq, thresholds, source_note):
    high_cov, high_ident, ok_cov, ok_ident, partial_cov = thresholds
    chain = d.get("chain", {})
    chain_seq = norm_seq(chain.get("sequence", ""))
    chain_len = int(chain.get("length", len(chain_seq)))
    local_chain = str(chain.get("chain_id", ""))

    raw = exact_subsequence_mapping(chain_seq, uni_seq)
    method = "sequence_alignment_exact_subsequence"
    score = len(chain_seq) * 2.0
    if raw is None:
        raw, score = pairwise_alignment_mapping(chain_seq, uni_seq)
        method = "sequence_alignment_pairwise_local"

    residue_lookup = {}
    for item in d.get("residue_mapping", []):
        try:
            residue_lookup[int(item.get("chain_index"))] = item
        except Exception:
            pass

    mapping = []
    matched = 0
    for cidx, upos, aa_match in raw:
        c_aa = chain_seq[cidx - 1] if 1 <= cidx <= len(chain_seq) else ""
        u_aa = uni_seq[upos - 1] if 1 <= upos <= len(uni_seq) else ""
        if aa_match:
            matched += 1
        rm = residue_lookup.get(cidx, {})
        mapping.append({
            "chain_index": cidx,
            "pdb_resseq": rm.get("pdb_resseq"),
            "pdb_icode": rm.get("pdb_icode", ""),
            "pdb_hetfield": rm.get("pdb_hetfield", ""),
            "chain_aa": c_aa,
            "uniprot_acc": acc,
            "uniprot_pos": upos,
            "uniprot_aa": u_aa,
            "aa_match": bool(aa_match),
        })

    mapped = len(mapping)
    cov = mapped / chain_len if chain_len else 0.0
    ident = matched / mapped if mapped else 0.0
    q, st = quality_status(cov, ident, high_cov, high_ident, ok_cov, ok_ident, partial_cov)
    if mapped == 0:
        q, st = "FAILED", "failed"

    unmapped = sorted(set(range(1, chain_len + 1)) - set(x["chain_index"] for x in mapping))
    reason = []
    if cov < ok_cov:
        reason.append(f"Low alignment coverage: {cov:.4f}")
    if ident < ok_ident:
        reason.append(f"Low alignment identity: {ident:.4f}")

    return {
        "status": st,
        "quality": q,
        "source": "UniProt sequence alignment",
        "method": method,
        "fallback_from": source_note,
        "coordinate_system": "chain_index_1_based",
        "chain_id": local_chain,
        "sifts_chain_id": "",
        "mapped_count": mapped,
        "chain_length": chain_len,
        "coverage": round(cov, 4),
        "identity": round(ident, 4),
        "alignment_score": round(score, 4),
        "uniprot_accessions": [acc] if acc else [],
        "unmapped_chain_index": unmapped,
        "mapping": mapping,
        "reason": "; ".join(reason),
    }


def best_alignment_fallback(d, accs, cache_dir, base, thresholds, max_candidates=80):
    candidates = []
    seen = set()
    for acc in accs:
        acc = str(acc or "").strip()
        if not acc or acc in seen:
            continue
        seen.add(acc)
        if len(seen) > max_candidates:
            break
        try:
            uni_seq = get_uniprot_sequence(acc, cache_dir, base=base)
            if not uni_seq:
                continue
            obj = alignment_mapping_object(d, acc, uni_seq, thresholds, "direct_sifts_residue_mapping_failed")
            candidates.append(obj)
        except Exception:
            continue

    candidates.sort(key=lambda x: (x.get("identity", 0), x.get("coverage", 0), x.get("mapped_count", 0)), reverse=True)
    if not candidates:
        chain_seq = norm_seq(d.get("chain", {}).get("sequence", ""))
        chain_len = int(d.get("chain", {}).get("length", len(chain_seq)))
        return {
            "status": "failed",
            "quality": "FAILED",
            "source": "UniProt sequence alignment",
            "method": "sequence_alignment_fallback",
            "coordinate_system": "chain_index_1_based",
            "chain_id": str(d.get("chain", {}).get("chain_id", "")),
            "sifts_chain_id": "",
            "mapped_count": 0,
            "chain_length": chain_len,
            "coverage": 0.0,
            "identity": 0.0,
            "uniprot_accessions": [],
            "unmapped_chain_index": list(range(1, chain_len + 1)),
            "mapping": [],
            "reason": "No candidate UniProt accession could be aligned.",
        }

    best = candidates[0]
    best["candidate_alignment_scores"] = [
        {
            "uniprot_acc": c.get("uniprot_accessions", [""])[0] if c.get("uniprot_accessions") else "",
            "method": c.get("method", ""),
            "mapped_count": c.get("mapped_count", 0),
            "coverage": c.get("coverage", 0),
            "identity": c.get("identity", 0),
            "quality": c.get("quality", ""),
        }
        for c in candidates[:20]
    ]
    return best


def final_mapping(d, idx, cache_dir, base, thresholds):
    direct = best_direct_sifts_mapping(d, idx, thresholds)
    if direct.get("quality") in {"HIGH", "ACCEPTABLE", "PARTIAL_USABLE"}:
        direct["final_mapping_method"] = "direct_sifts_residue_mapping"
        return direct, direct

    accs = list(idx.get("uniprot_accessions", []))
    fb = best_alignment_fallback(d, accs, cache_dir, base, thresholds)
    fb["final_mapping_method"] = "sequence_alignment_fallback"
    return fb, direct


def feature_position(feature_part):
    if not isinstance(feature_part, dict):
        return None
    try:
        return int(feature_part.get("value"))
    except Exception:
        return None


def extract_uniprot_features(uj, feature_types, include_all=False):
    out = []
    for f in uj.get("features", []):
        if not isinstance(f, dict):
            continue
        ftype = str(f.get("type", "")).strip()
        if not include_all and feature_types and ftype not in feature_types:
            continue
        loc = f.get("location", {})
        if not isinstance(loc, dict):
            continue
        start = feature_position(loc.get("start", {}))
        end = feature_position(loc.get("end", {}))
        if start is None and end is not None:
            start = end
        if end is None and start is not None:
            end = start
        if start is None or end is None:
            continue
        if start > end:
            start, end = end, start
        out.append({
            "type": ftype,
            "description": str(f.get("description", "") or "").strip(),
            "uniprot_start": start,
            "uniprot_end": end,
            "raw_feature": f,
        })
    return out


def build_uniprot_pos_to_chain(mapping_obj, acc=None):
    out = defaultdict(list)
    for m in mapping_obj.get("mapping", []):
        if acc is not None and m.get("uniprot_acc") != acc:
            continue
        try:
            out[int(m["uniprot_pos"])].append(int(m["chain_index"]))
        except Exception:
            pass
    return out


def fragments(positions):
    positions = sorted(set(int(x) for x in positions))
    if not positions:
        return []
    out = []
    s = p = positions[0]
    for x in positions[1:]:
        if x == p + 1:
            p = x
        else:
            out.append({"chain_start": s, "chain_end": p})
            s = p = x
    out.append({"chain_start": s, "chain_end": p})
    return out


def convert_uniprot_features(mapping_obj, feature_by_acc):
    items = []
    for acc in mapping_obj.get("uniprot_accessions", []):
        pos_to_chain = build_uniprot_pos_to_chain(mapping_obj, acc=acc)
        for f in feature_by_acc.get(acc, []):
            positions = []
            for upos in range(f["uniprot_start"], f["uniprot_end"] + 1):
                positions.extend(pos_to_chain.get(upos, []))
            positions = sorted(set(positions))
            if not positions:
                continue
            items.append({
                "source": "UniProt",
                "database": "UniProt",
                "uniprot_acc": acc,
                "feature_type": f["type"],
                "description": f["description"],
                "uniprot_start": f["uniprot_start"],
                "uniprot_end": f["uniprot_end"],
                "chain_start": min(positions),
                "chain_end": max(positions),
                "mapped_chain_positions": positions,
                "fragments": fragments(positions),
                "mapped_count": len(positions),
            })
    return {"status": "done" if items else "empty", "source": "UniProt REST JSON + final mapping", "coordinate_system": "chain_index_1_based", "items": items, "item_count": len(items)}


def convert_sifts_xrefs(mapping_obj, xref_rows, database):
    pos_to_chain = build_uniprot_pos_to_chain(mapping_obj, acc=None)
    final_accs = set(mapping_obj.get("uniprot_accessions", []))
    groups = {}
    for r in xref_rows:
        if r.get("database") != database:
            continue
        try:
            upos = int(r.get("uniprot_pos"))
        except Exception:
            continue
        xacc = str(r.get("uniprot_acc", "") or "").strip()
        if xacc and final_accs and xacc not in final_accs:
            continue
        positions = pos_to_chain.get(upos, [])
        if not positions:
            continue
        accession = str(r.get("accession", "") or "").strip()
        if not accession:
            continue
        key = (database, accession)
        if key not in groups:
            groups[key] = {"positions": set(), "uniprot_positions": set(), "evidence": set()}
        for p in positions:
            groups[key]["positions"].add(int(p))
        groups[key]["uniprot_positions"].add(upos)
        ev = str(r.get("evidence", "") or "").strip()
        if ev:
            groups[key]["evidence"].add(ev)

    items = []
    for (_, accession), g in sorted(groups.items(), key=lambda kv: kv[0]):
        pos = sorted(g["positions"])
        upos = sorted(g["uniprot_positions"])
        if not pos:
            continue
        item = {
            "source": f"SIFTS {database} crossRefDb",
            "database": database,
            "feature_type": database,
            "accession": accession,
            "uniprot_start": min(upos) if upos else None,
            "uniprot_end": max(upos) if upos else None,
            "chain_start": min(pos),
            "chain_end": max(pos),
            "mapped_chain_positions": pos,
            "fragments": fragments(pos),
            "mapped_count": len(pos),
        }
        ev = sorted(x for x in g["evidence"] if x)
        if ev:
            item["evidence"] = ev
        items.append(item)
    return {"status": "done" if items else "empty", "source": f"SIFTS XML crossRefDb dbSource={database}; converted by final mapping", "coordinate_system": "chain_index_1_based", "items": items, "item_count": len(items)}


def merge_external(uniprot_features, pfam, interpro):
    all_items = []
    for section, obj in [("uniprot_features", uniprot_features), ("pfam", pfam), ("interpro", interpro)]:
        for item in obj.get("items", []):
            x = dict(item)
            x["section"] = section
            all_items.append(x)
    return {
        "status": "done" if all_items else "empty",
        "source": "UniProt REST + SIFTS Pfam/InterPro using final mapping",
        "coordinate_system": "chain_index_1_based",
        "uniprot_features": uniprot_features,
        "pfam": pfam,
        "interpro": interpro,
        "items": all_items,
        "item_count": len(all_items),
        "counts": {"uniprot_features": uniprot_features.get("item_count", 0), "pfam": pfam.get("item_count", 0), "interpro": interpro.get("item_count", 0)},
    }


def update_one_json(json_path, idx, cache_dir, base, feature_types, include_all, thresholds, dry_run=False):
    d = read_json(json_path)
    entry_id = d.get("entry", {}).get("entry_id", json_path.parent.parent.name)
    pdb_id = str(d.get("entry", {}).get("pdb_id", "")).strip().lower()
    chain_id = str(d.get("chain", {}).get("chain_id", json_path.stem))

    fmap, direct_qc = final_mapping(d, idx, cache_dir, base, thresholds)

    trusted_mapping = fmap.get("quality") in {"HIGH", "ACCEPTABLE", "PARTIAL_USABLE"}

    feature_by_acc = {}
    if trusted_mapping:
        for acc in fmap.get("uniprot_accessions", []):
            try:
                uj = fetch_uniprot_json(acc, cache_dir, base=base)
                feature_by_acc[acc] = extract_uniprot_features(uj, feature_types, include_all=include_all)
            except Exception:
                feature_by_acc[acc] = []

        uf = convert_uniprot_features(fmap, feature_by_acc)
        pfam = convert_sifts_xrefs(fmap, idx.get("xref_rows", []), "Pfam")
        interpro = convert_sifts_xrefs(fmap, idx.get("xref_rows", []), "InterPro")
        external = merge_external(uf, pfam, interpro)
    else:
        uf = {
            "status": "unreliable_mapping",
            "source": "not generated because final UniProt mapping is unreliable",
            "coordinate_system": "chain_index_1_based",
            "items": [],
            "item_count": 0,
        }
        pfam = {
            "status": "unreliable_mapping",
            "source": "not generated because final UniProt mapping is unreliable",
            "coordinate_system": "chain_index_1_based",
            "items": [],
            "item_count": 0,
        }
        interpro = {
            "status": "unreliable_mapping",
            "source": "not generated because final UniProt mapping is unreliable",
            "coordinate_system": "chain_index_1_based",
            "items": [],
            "item_count": 0,
        }
        external = {
            "status": "unreliable_mapping",
            "source": "external annotations skipped because final UniProt mapping is unreliable",
            "coordinate_system": "chain_index_1_based",
            "uniprot_features": uf,
            "pfam": pfam,
            "interpro": interpro,
            "items": [],
            "item_count": 0,
            "counts": {
                "uniprot_features": 0,
                "pfam": 0,
                "interpro": 0,
            },
        }

    d["sifts_residue_mapping_qc"] = direct_qc
    d["uniprot_mapping"] = fmap
    d["external_annotations"] = external
    d.setdefault("tracks", {})
    d["tracks"]["uniprot_features"] = {"status": uf["status"], "type": "features", "source": uf["source"], "coordinate_system": "chain_index_1_based", "items": uf["items"]}
    d["tracks"]["pfam"] = {"status": pfam["status"], "type": "features", "source": pfam["source"], "coordinate_system": "chain_index_1_based", "items": pfam["items"]}
    d["tracks"]["interpro"] = {"status": interpro["status"], "type": "features", "source": interpro["source"], "coordinate_system": "chain_index_1_based", "items": interpro["items"]}

    d.setdefault("status", {})
    d["status"]["uniprot_mapping"] = fmap.get("status", "")
    d["status"]["uniprot_mapping_quality"] = fmap.get("quality", "")
    d["status"]["uniprot_mapping_method"] = fmap.get("final_mapping_method", fmap.get("method", ""))
    d["status"]["external_annotations"] = external.get("status", "")
    d["status"]["uniprot_features"] = uf.get("status", "")
    d["status"]["pfam"] = pfam.get("status", "")
    d["status"]["interpro"] = interpro.get("status", "")

    append_log(d, "uniprot_alignment_fallback_remap", external.get("status", ""),
               f"Final mapping method={fmap.get('final_mapping_method')}; quality={fmap.get('quality')}; coverage={fmap.get('coverage')}; identity={fmap.get('identity')}; direct_quality={direct_qc.get('quality')}; direct_cov={direct_qc.get('coverage')}; direct_id={direct_qc.get('identity')}; pfam={pfam.get('item_count')}; interpro={interpro.get('item_count')}")

    if not dry_run:
        write_json_atomic(json_path, d)

    return {
        "entry_id": entry_id,
        "pdb_id": pdb_id,
        "chain_id": chain_id,
        "json_path": str(json_path),
        "status": fmap.get("status", ""),
        "mapping_quality": fmap.get("quality", ""),
        "mapping_method": fmap.get("final_mapping_method", fmap.get("method", "")),
        "uniprot_accessions": ";".join(fmap.get("uniprot_accessions", [])),
        "mapped_count": fmap.get("mapped_count", 0),
        "chain_length": fmap.get("chain_length", ""),
        "coverage": fmap.get("coverage", 0),
        "identity": fmap.get("identity", 0),
        "direct_sifts_status": direct_qc.get("status", ""),
        "direct_sifts_quality": direct_qc.get("quality", ""),
        "direct_sifts_chain": direct_qc.get("sifts_chain_id", ""),
        "direct_sifts_coverage": direct_qc.get("coverage", 0),
        "direct_sifts_identity": direct_qc.get("identity", 0),
        "external_annotation_status": external.get("status", ""),
        "external_annotation_count": external.get("item_count", 0),
        "uniprot_feature_status": uf.get("status", ""),
        "uniprot_feature_count": uf.get("item_count", 0),
        "pfam_status": pfam.get("status", ""),
        "pfam_count": pfam.get("item_count", 0),
        "interpro_status": interpro.get("status", ""),
        "interpro_count": interpro.get("item_count", 0),
        "note": fmap.get("reason", ""),
    }


def write_summary(path, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "entry_id", "pdb_id", "chain_id", "json_path", "status", "mapping_quality", "mapping_method",
        "uniprot_accessions", "mapped_count", "chain_length", "coverage", "identity",
        "direct_sifts_status", "direct_sifts_quality", "direct_sifts_chain", "direct_sifts_coverage", "direct_sifts_identity",
        "external_annotation_status", "external_annotation_count", "uniprot_feature_status", "uniprot_feature_count",
        "pfam_status", "pfam_count", "interpro_status", "interpro_count", "note",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-script", default="scripts/07_sifts_uniprot_annotations.py")
    ap.add_argument("--manifest", default="example_outputs/manifest/manifest_skeleton.csv")
    ap.add_argument("--json-root", default="example_outputs/json")
    ap.add_argument("--cache-dir", default=".cache/sifts_uniprot")
    ap.add_argument("--summary-out", default="example_outputs/qc/sifts_uniprot_alignment_fallback_test100.csv")
    ap.add_argument("--limit-entries", type=int, default=0)
    ap.add_argument("--include-all-uniprot-features", action="store_true")
    ap.add_argument("--feature-types", default="")
    ap.add_argument("--high-coverage", type=float, default=0.95)
    ap.add_argument("--high-identity", type=float, default=0.98)
    ap.add_argument("--ok-coverage", type=float, default=0.80)
    ap.add_argument("--ok-identity", type=float, default=0.90)
    ap.add_argument("--partial-coverage", type=float, default=0.50)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    base = load_base_module(args.base_script)
    thresholds = (args.high_coverage, args.high_identity, args.ok_coverage, args.ok_identity, args.partial_coverage)
    feature_types = {x.strip() for x in args.feature_types.split(",") if x.strip()} if args.feature_types.strip() else DEFAULT_FEATURE_TYPES

    manifest = Path(args.manifest)
    json_root = Path(args.json_root)
    if manifest.exists():
        entries = get_entries_from_manifest(manifest, args.limit_entries)
    else:
        entries = get_entries_from_json_root(json_root, args.limit_entries)

    print("========== UniProt remap with alignment fallback ==========")
    print("entries:", len(entries))
    print("json_root:", json_root)
    print("summary_out:", args.summary_out)
    print("dry_run:", args.dry_run)
    print("==========================================================")

    rows = []
    sifts_cache = {}
    for i, entry in enumerate(entries, start=1):
        chain_jsons = find_chain_jsons(json_root, entry)
        if not chain_jsons:
            continue
        first = read_json(chain_jsons[0])
        pdb_id = str(first.get("entry", {}).get("pdb_id", "")).strip().lower()
        print(f"[ENTRY] {i}/{len(entries)} {entry} pdb={pdb_id} chains={len(chain_jsons)}")

        if not pdb_id or pdb_id.upper() == "NA":
            continue

        try:
            idx = load_sifts_index(pdb_id, Path(args.cache_dir), base, sifts_cache)
        except Exception as e:
            print("  [ERROR] load SIFTS failed:", repr(e))
            continue

        for j, jf in enumerate(chain_jsons, start=1):
            print(f"  [CHAIN] {j}/{len(chain_jsons)} {jf.name}")
            try:
                rows.append(update_one_json(jf, idx, Path(args.cache_dir), base, feature_types, args.include_all_uniprot_features, thresholds, dry_run=args.dry_run))
            except Exception as e:
                d = read_json(jf)
                rows.append({
                    "entry_id": entry,
                    "pdb_id": pdb_id,
                    "chain_id": d.get("chain", {}).get("chain_id", jf.stem),
                    "json_path": str(jf),
                    "status": "failed",
                    "mapping_quality": "FAILED",
                    "mapping_method": "",
                    "uniprot_accessions": "",
                    "mapped_count": 0,
                    "chain_length": d.get("chain", {}).get("length", ""),
                    "coverage": 0,
                    "identity": 0,
                    "direct_sifts_status": "",
                    "direct_sifts_quality": "",
                    "direct_sifts_chain": "",
                    "direct_sifts_coverage": 0,
                    "direct_sifts_identity": 0,
                    "external_annotation_status": "",
                    "external_annotation_count": 0,
                    "uniprot_feature_status": "",
                    "uniprot_feature_count": 0,
                    "pfam_status": "",
                    "pfam_count": 0,
                    "interpro_status": "",
                    "interpro_count": 0,
                    "note": repr(e),
                })
                print("    [ERROR]", repr(e))

    write_summary(args.summary_out, rows)

    c_status = Counter(r["status"] for r in rows)
    c_quality = Counter(r["mapping_quality"] for r in rows)
    c_method = Counter(r["mapping_method"] for r in rows)
    c_ext = Counter(r["external_annotation_status"] for r in rows)
    c_pfam = Counter(r["pfam_status"] for r in rows)
    c_interpro = Counter(r["interpro_status"] for r in rows)

    print("========== Finished ==========")
    print("chains processed:", len(rows))
    print("mapping status:", dict(c_status))
    print("mapping quality:", dict(c_quality))
    print("mapping method:", dict(c_method))
    print("external annotation:", dict(c_ext))
    print("pfam:", dict(c_pfam))
    print("interpro:", dict(c_interpro))
    print("summary:", args.summary_out)
    print("==============================")


if __name__ == "__main__":
    main()
