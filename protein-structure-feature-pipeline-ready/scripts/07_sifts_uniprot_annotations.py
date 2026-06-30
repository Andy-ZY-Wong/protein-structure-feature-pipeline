#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
07_sifts_uniprot_annotations.py

Purpose:
1. Read existing chain JSON files;
2. Use pdb_id / chain_id / residue_mapping from JSON;
3. Download and parse SIFTS XML to build PDB residue -> UniProt residue mappings;
4. Query UniProt JSON and extract UniProt sequence features;
5. Convert UniProt feature coordinates back to chain_index_1_based;
6. Write results back to chain JSON:
   - uniprot_mapping
   - external_annotations

Local structure features are not recomputed.
Files such as .fullLen / .manAlign / .toRCSB are not used.
If the local chain_id does not exactly match the SIFTS chain_id, the script tries to select the best SIFTS chain by residue number.
"""

import argparse
import csv
import gzip
import json
import os
import re
import time
import urllib.request
import urllib.error
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any


DEFAULT_FEATURE_TYPES = {
    # molecule processing
    "Initiator methionine",
    "Signal peptide",
    "Transit peptide",
    "Propeptide",
    "Chain",
    "Peptide",

    # regions
    "Domain",
    "Repeat",
    "Region",
    "Coiled coil",
    "Compositional bias",
    "Motif",
    "Zinc finger",
    "DNA binding",
    "Nucleotide binding",
    "Calcium binding",
    "Topological domain",
    "Transmembrane",
    "Intramembrane",

    # sites / PTM
    "Active site",
    "Binding site",
    "Metal binding",
    "Site",
    "Modified residue",
    "Lipidation",
    "Glycosylation",
    "Disulfide bond",
    "Cross-link",
}


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O", "ASX": "B", "GLX": "Z", "UNK": "X",
}


# ============================================================
# Basic utilities
# ============================================================

def now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def read_json(path: Path) -> Dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def write_json_atomic(path: Path, data: Dict) -> None:
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


def strip_ns(tag: str) -> str:
    return tag.split("}", 1)[-1] if "}" in tag else tag


def attr_get(elem, names: List[str], default: str = "") -> str:
    for n in names:
        if n in elem.attrib:
            return elem.attrib.get(n, default)
    # Handle case differences
    lower_map = {k.lower(): v for k, v in elem.attrib.items()}
    for n in names:
        if n.lower() in lower_map:
            return lower_map[n.lower()]
    return default


def normalize_pdb_id(x: str) -> str:
    return str(x or "").strip().lower()


def normalize_chain_id(x: str) -> str:
    return str(x or "").strip()


def parse_resnum_icode(x: Any) -> Tuple[Optional[int], str]:
    """
    Supports:
      123
      123A
      -1
      123?
    """
    s = str(x or "").strip()
    if not s or s in {".", "?", "None", "null"}:
        return None, ""

    m = re.match(r"^(-?\d+)([A-Za-z]?)$", s)
    if not m:
        m = re.search(r"(-?\d+)([A-Za-z]?)", s)
    if not m:
        return None, ""

    return int(m.group(1)), (m.group(2) or "").strip()


def one_letter_from_resname(x: str) -> str:
    x = str(x or "").strip().upper()
    if len(x) == 1:
        return x
    return AA3_TO_1.get(x, "X")


def get_feature_pos(loc_part: Any) -> Optional[int]:
    """
    UniProt JSON location:
      {"start": {"value": 1, "modifier": "EXACT"}, "end": {...}}
    Sometimes value may be missing or stored as a string.
    """
    if not isinstance(loc_part, dict):
        return None
    val = loc_part.get("value")
    try:
        return int(val)
    except Exception:
        return None


def feature_range_from_uniprot_feature(feature: Dict) -> Tuple[Optional[int], Optional[int]]:
    loc = feature.get("location", {})
    if not isinstance(loc, dict):
        return None, None

    start = get_feature_pos(loc.get("start", {}))
    end = get_feature_pos(loc.get("end", {}))

    if start is None and end is not None:
        start = end
    if end is None and start is not None:
        end = start

    if start is None or end is None:
        return None, None

    if start > end:
        start, end = end, start

    return start, end


# ============================================================
# Read manifest / JSON
# ============================================================

def read_manifest(path: Path) -> List[Dict]:
    with open(path, newline="", encoding="utf-8", errors="ignore") as f:
        return list(csv.DictReader(f))


def get_entries_from_manifest(path: Path, limit_entries: int = 0) -> List[str]:
    rows = read_manifest(path)
    entries = []
    seen = set()

    for row in rows:
        entry_id = row.get("entry_id", "").strip()
        if not entry_id or entry_id in seen:
            continue

        # Prefer structures with existing skeletons
        has_structure = row.get("has_structure", "").lower()
        skeleton_status = row.get("skeleton_status", "").lower()

        if has_structure and has_structure not in {"yes", "true", "1"}:
            continue
        if skeleton_status and skeleton_status != "done":
            continue

        entries.append(entry_id)
        seen.add(entry_id)

        if limit_entries and len(entries) >= limit_entries:
            break

    return entries


def get_entries_from_json_root(json_root: Path, limit_entries: int = 0) -> List[str]:
    entries = []
    for p in sorted(json_root.iterdir()):
        if not p.is_dir():
            continue
        if (p / "chains").is_dir():
            entries.append(p.name)
        if limit_entries and len(entries) >= limit_entries:
            break
    return entries


def find_chain_jsons(json_root: Path, entry_id: str) -> List[Path]:
    chain_dir = json_root / entry_id / "chains"
    if not chain_dir.is_dir():
        return []
    return sorted(chain_dir.glob("*.json"))


# ============================================================
# SIFTS download and parsing
# ============================================================

def download_url(url: str, out_path: Path, retries: int = 3, sleep_sec: float = 1.0) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists() and out_path.stat().st_size > 0:
        return

    last_error = ""
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(
                url,
                headers={"User-Agent": "sequence-track-sifts-uniprot/1.0"}
            )
            with urllib.request.urlopen(req, timeout=60) as resp:
                data = resp.read()
            with open(out_path, "wb") as f:
                f.write(data)
            return
        except Exception as e:
            last_error = repr(e)
            time.sleep(sleep_sec)

    raise RuntimeError(f"Failed to download {url}: {last_error}")


def fetch_sifts_xml_gz(pdb_id: str, cache_dir: Path) -> Path:
    """
    Recommended SIFTS XML path:
      split_xml/{second and third characters}/{pdb}.xml.gz
    fallback：
      xml/{pdb}.xml.gz
    """
    pdb_id = normalize_pdb_id(pdb_id)
    if len(pdb_id) != 4:
        raise RuntimeError(f"Invalid PDB ID for SIFTS: {pdb_id}")

    out = cache_dir / "sifts_xml" / f"{pdb_id}.xml.gz"
    if out.exists() and out.stat().st_size > 0:
        return out

    mid = pdb_id[1:3]
    urls = [
        f"https://ftp.ebi.ac.uk/pubpublic_sifts_mirror/split_xml/{mid}/{pdb_id}.xml.gz",
        f"https://ftp.ebi.ac.uk/pubpublic_sifts_mirror/xml/{pdb_id}.xml.gz",
    ]

    last_error = ""
    for url in urls:
        try:
            download_url(url, out)
            return out
        except Exception as e:
            last_error = repr(e)

    raise RuntimeError(f"No SIFTS XML found for {pdb_id}: {last_error}")


def parse_sifts_xml(xml_gz_path: Path) -> List[Dict]:
    """
    Return residue-level mapping rows:
      [
        {
          "pdb_chain": "A",
          "pdb_resseq": 2,
          "pdb_icode": "",
          "pdb_resname": "SER",
          "uniprot_acc": "...",
          "uniprot_pos": 2,
          "uniprot_resname": "S"
        }
      ]

    SIFTS XML residue nodes usually contain crossRefDb entries with dbSource="PDB" and dbSource="UniProt".
    """
    with gzip.open(xml_gz_path, "rb") as f:
        tree = ET.parse(f)

    root = tree.getroot()
    rows = []

    for elem in root.iter():
        if strip_ns(elem.tag) != "residue":
            continue

        pdb_refs = []
        uniprot_refs = []

        for child in list(elem):
            if strip_ns(child.tag) != "crossRefDb":
                continue

            src = attr_get(child, ["dbSource"]).strip()
            if src == "PDB":
                chain = attr_get(child, ["dbChainId", "dbChainID", "db_chain_id"])
                resnum = attr_get(child, ["dbResNum", "dbResnum", "db_res_num"])
                resname = attr_get(child, ["dbResName", "dbResname", "db_res_name"])
                pdb_refs.append({
                    "chain": normalize_chain_id(chain),
                    "resnum": resnum,
                    "resname": resname
                })

            elif src == "UniProt":
                acc = attr_get(child, ["dbAccessionId", "dbAccessionID", "dbAccession"])
                resnum = attr_get(child, ["dbResNum", "dbResnum", "db_res_num"])
                resname = attr_get(child, ["dbResName", "dbResname", "db_res_name"])
                uniprot_refs.append({
                    "acc": acc.strip(),
                    "resnum": resnum,
                    "resname": resname
                })

        for p in pdb_refs:
            pdb_resseq, pdb_icode = parse_resnum_icode(p["resnum"])
            if pdb_resseq is None or not p["chain"]:
                continue

            for u in uniprot_refs:
                upos, _uicode = parse_resnum_icode(u["resnum"])
                if upos is None or not u["acc"]:
                    continue

                rows.append({
                    "pdb_chain": p["chain"],
                    "pdb_resseq": pdb_resseq,
                    "pdb_icode": pdb_icode,
                    "pdb_resname": p["resname"],
                    "pdb_aa": one_letter_from_resname(p["resname"]),
                    "uniprot_acc": u["acc"],
                    "uniprot_pos": upos,
                    "uniprot_resname": u["resname"],
                    "uniprot_aa": one_letter_from_resname(u["resname"]),
                })

    return rows


def build_sifts_index(rows: List[Dict]) -> Dict:
    """
    Build indexes:
      by_chain_res[(chain, resseq, icode)] -> row
      by_chain_res_noicode[(chain, resseq)] -> row
      chains -> set
    """
    idx = {
        "by_chain_res": {},
        "by_chain_res_noicode": {},
        "chains": set(),
        "rows_by_chain": defaultdict(list),
    }

    for r in rows:
        chain = r["pdb_chain"]
        resseq = r["pdb_resseq"]
        icode = r.get("pdb_icode", "")

        idx["chains"].add(chain)
        idx["rows_by_chain"][chain].append(r)
        idx["by_chain_res"][(chain, resseq, icode)] = r
        idx["by_chain_res_noicode"][(chain, resseq)] = r

    return idx


def infer_best_sifts_chain(chain_json: Dict, sifts_idx: Dict) -> Tuple[str, int]:
    """
    If local chain_id cannot directly match a SIFTS chain, choose the best chain using residue number + amino acid statistics.
    """
    residue_mapping = chain_json.get("residue_mapping", [])
    scores = Counter()

    for item in residue_mapping:
        try:
            resseq = int(item.get("pdb_resseq"))
        except Exception:
            continue

        aa = str(item.get("one_letter", "")).upper()

        for chain in sifts_idx["chains"]:
            r = sifts_idx["by_chain_res_noicode"].get((chain, resseq))
            if not r:
                continue
            score = 1
            if aa and r.get("pdb_aa") and aa == r.get("pdb_aa"):
                score += 2
            scores[chain] += score

    if not scores:
        return "", 0

    best_chain, best_score = scores.most_common(1)[0]
    return best_chain, best_score


# ============================================================
# UniProt JSON download and feature parsing
# ============================================================

def fetch_uniprot_json(acc: str, cache_dir: Path) -> Path:
    acc = str(acc or "").strip()
    if not acc:
        raise RuntimeError("Empty UniProt accession.")

    out = cache_dir / "uniprot_json" / f"{acc}.json"
    if out.exists() and out.stat().st_size > 0:
        return out

    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    download_url(url, out)
    return out


def parse_uniprot_json(path: Path) -> Dict:
    return read_json(path)


def get_uniprot_sequence(uj: Dict) -> str:
    seq = uj.get("sequence", {})
    if isinstance(seq, dict):
        return str(seq.get("value", "")).strip().upper()
    return ""


def get_uniprot_primary_accession(uj: Dict, fallback: str) -> str:
    return str(uj.get("primaryAccession") or fallback).strip()


def extract_uniprot_features(
    uj: Dict,
    feature_types: Optional[set] = None,
    include_all: bool = False
) -> List[Dict]:
    features = uj.get("features", [])
    out = []

    for f in features:
        if not isinstance(f, dict):
            continue

        ftype = str(f.get("type", "")).strip()
        if not include_all:
            if feature_types and ftype not in feature_types:
                continue

        start, end = feature_range_from_uniprot_feature(f)
        if start is None or end is None:
            continue

        desc = str(f.get("description", "") or "").strip()

        out.append({
            "type": ftype,
            "description": desc,
            "uniprot_start": start,
            "uniprot_end": end,
            "raw_feature": f
        })

    return out


# ============================================================
# Map and write back to JSON
# ============================================================

def select_sifts_chain_for_json(chain_json: Dict, sifts_idx: Dict) -> Tuple[str, str]:
    local_chain = normalize_chain_id(chain_json.get("chain", {}).get("chain_id", ""))

    if local_chain in sifts_idx["chains"]:
        return local_chain, "direct_chain_id_match"

    inferred, score = infer_best_sifts_chain(chain_json, sifts_idx)
    if inferred:
        return inferred, f"inferred_by_residue_number_score_{score}"

    return "", "no_sifts_chain_match"


def map_chain_json_to_uniprot(chain_json: Dict, sifts_idx: Dict) -> Tuple[Dict, str]:
    """
    Use SIFTS to build chain_index -> UniProt position mapping for the current chain JSON.
    """
    chain = chain_json.get("chain", {})
    chain_id = normalize_chain_id(chain.get("chain_id", ""))
    chain_seq = str(chain.get("sequence", "")).strip().upper()
    chain_len = int(chain.get("length", len(chain_seq)))

    sifts_chain, chain_method = select_sifts_chain_for_json(chain_json, sifts_idx)

    mapping_items = []

    if not sifts_chain:
        return {
            "status": "failed",
            "source": "SIFTS",
            "method": chain_method,
            "coordinate_system": "chain_index_1_based",
            "chain_id": chain_id,
            "sifts_chain_id": "",
            "mapped_count": 0,
            "coverage": 0.0,
            "identity": 0.0,
            "uniprot_accessions": [],
            "mapping": [],
            "unmapped_chain_index": list(range(1, chain_len + 1)),
            "reason": "Could not match local chain to any SIFTS chain."
        }, ""

    matched_aa = 0

    for item in chain_json.get("residue_mapping", []):
        try:
            chain_index = int(item.get("chain_index"))
            pdb_resseq = int(item.get("pdb_resseq"))
        except Exception:
            continue

        pdb_icode = str(item.get("pdb_icode", "") or "").strip()
        chain_aa = str(item.get("one_letter", "") or "").upper()

        srow = sifts_idx["by_chain_res"].get((sifts_chain, pdb_resseq, pdb_icode))
        if not srow:
            srow = sifts_idx["by_chain_res_noicode"].get((sifts_chain, pdb_resseq))

        if not srow:
            continue

        uniprot_aa = srow.get("uniprot_aa", "")
        aa_match = bool(chain_aa and uniprot_aa and chain_aa == uniprot_aa)

        if aa_match:
            matched_aa += 1

        mapping_items.append({
            "chain_index": chain_index,
            "pdb_chain": sifts_chain,
            "pdb_resseq": pdb_resseq,
            "pdb_icode": pdb_icode,
            "chain_aa": chain_aa,
            "uniprot_acc": srow["uniprot_acc"],
            "uniprot_pos": srow["uniprot_pos"],
            "uniprot_aa": uniprot_aa,
            "aa_match": aa_match
        })

    mapping_items = sorted(mapping_items, key=lambda x: x["chain_index"])

    mapped_count = len(mapping_items)
    coverage = mapped_count / chain_len if chain_len else 0.0
    identity = matched_aa / mapped_count if mapped_count else 0.0
    unmapped = sorted(set(range(1, chain_len + 1)) - set(x["chain_index"] for x in mapping_items))
    accs = sorted(set(x["uniprot_acc"] for x in mapping_items if x.get("uniprot_acc")))

    if mapped_count == 0:
        status = "failed"
        reason = "SIFTS chain was matched, but no residue-level mapping could be written."
    elif coverage >= 0.80:
        status = "done"
        reason = ""
    else:
        status = "partial"
        reason = f"Low SIFTS coverage: {coverage:.4f}"

    obj = {
        "status": status,
        "source": "SIFTS",
        "method": chain_method,
        "coordinate_system": "chain_index_1_based",
        "chain_id": chain_id,
        "sifts_chain_id": sifts_chain,
        "mapped_count": mapped_count,
        "chain_length": chain_len,
        "coverage": round(coverage, 4),
        "identity": round(identity, 4),
        "uniprot_accessions": accs,
        "unmapped_chain_index": unmapped,
        "mapping": mapping_items
    }

    if reason:
        obj["reason"] = reason

    primary_acc = accs[0] if accs else ""
    return obj, primary_acc


def build_uniprot_pos_to_chain_index(uniprot_mapping: Dict, acc: str) -> Dict[int, List[int]]:
    out = defaultdict(list)

    for m in uniprot_mapping.get("mapping", []):
        if m.get("uniprot_acc") != acc:
            continue
        try:
            upos = int(m["uniprot_pos"])
            cidx = int(m["chain_index"])
        except Exception:
            continue
        out[upos].append(cidx)

    return out


def convert_uniprot_features_to_chain_annotations(
    uniprot_mapping: Dict,
    uniprot_feature_by_acc: Dict[str, List[Dict]]
) -> Dict:
    """
    Convert UniProt feature coordinates to chain_index coordinates.
    """
    items = []
    accs = uniprot_mapping.get("uniprot_accessions", [])

    for acc in accs:
        pos_to_chain = build_uniprot_pos_to_chain_index(uniprot_mapping, acc)
        if not pos_to_chain:
            continue

        features = uniprot_feature_by_acc.get(acc, [])

        for f in features:
            start = f["uniprot_start"]
            end = f["uniprot_end"]

            chain_positions = []
            for upos in range(start, end + 1):
                chain_positions.extend(pos_to_chain.get(upos, []))

            chain_positions = sorted(set(chain_positions))
            if not chain_positions:
                continue

            # Check whether the feature is continuous in the current structure chain; keep discontinuous features but mark fragments
            fragments = []
            frag_start = None
            prev = None

            for p in chain_positions:
                if frag_start is None:
                    frag_start = p
                    prev = p
                elif p == prev + 1:
                    prev = p
                else:
                    fragments.append({"chain_start": frag_start, "chain_end": prev})
                    frag_start = p
                    prev = p

            if frag_start is not None:
                fragments.append({"chain_start": frag_start, "chain_end": prev})

            items.append({
                "source": "UniProt",
                "uniprot_acc": acc,
                "feature_type": f["type"],
                "description": f["description"],
                "uniprot_start": start,
                "uniprot_end": end,
                "chain_start": min(chain_positions),
                "chain_end": max(chain_positions),
                "mapped_chain_positions": chain_positions,
                "fragments": fragments
            })

    status = "done" if items else "empty"

    return {
        "status": status,
        "source": "UniProt REST JSON + SIFTS coordinate mapping",
        "coordinate_system": "chain_index_1_based",
        "items": items,
        "item_count": len(items)
    }


def update_chain_json(
    json_path: Path,
    sifts_idx: Dict,
    cache_dir: Path,
    feature_types: Optional[set],
    include_all_features: bool,
) -> Dict:
    d = read_json(json_path)

    entry_id = d.get("entry", {}).get("entry_id", json_path.parent.parent.name)
    pdb_id = normalize_pdb_id(d.get("entry", {}).get("pdb_id", ""))
    chain_id = normalize_chain_id(d.get("chain", {}).get("chain_id", json_path.stem))

    uniprot_mapping, primary_acc = map_chain_json_to_uniprot(d, sifts_idx)

    uniprot_feature_by_acc = {}
    feature_fetch_notes = []

    for acc in uniprot_mapping.get("uniprot_accessions", []):
        try:
            uj_path = fetch_uniprot_json(acc, cache_dir)
            uj = parse_uniprot_json(uj_path)
            primary = get_uniprot_primary_accession(uj, acc)
            features = extract_uniprot_features(
                uj,
                feature_types=feature_types,
                include_all=include_all_features
            )
            uniprot_feature_by_acc[acc] = features
            feature_fetch_notes.append(f"{acc}:{len(features)}")
        except Exception as e:
            uniprot_feature_by_acc[acc] = []
            feature_fetch_notes.append(f"{acc}:failed:{repr(e)}")

    external_annotations = convert_uniprot_features_to_chain_annotations(
        uniprot_mapping=uniprot_mapping,
        uniprot_feature_by_acc=uniprot_feature_by_acc
    )

    d["uniprot_mapping"] = uniprot_mapping
    d["external_annotations"] = external_annotations

    # Also place results under tracks for downstream/frontend reading
    d.setdefault("tracks", {})
    d["tracks"]["uniprot_features"] = {
        "status": external_annotations["status"],
        "type": "features",
        "source": external_annotations["source"],
        "coordinate_system": "chain_index_1_based",
        "items": external_annotations["items"]
    }

    d.setdefault("status", {})
    d["status"]["uniprot_mapping"] = uniprot_mapping["status"]
    d["status"]["external_annotations"] = external_annotations["status"]

    overall = "done" if uniprot_mapping["status"] == "done" else uniprot_mapping["status"]

    append_log(
        d,
        step="sifts_uniprot_annotations",
        status=overall,
        message=(
            f"SIFTS mapping and UniProt annotations added. "
            f"pdb_id={pdb_id}; chain={chain_id}; "
            f"mapped={uniprot_mapping.get('mapped_count')}; "
            f"coverage={uniprot_mapping.get('coverage')}; "
            f"uniprot={','.join(uniprot_mapping.get('uniprot_accessions', []))}; "
            f"features={';'.join(feature_fetch_notes)}"
        )
    )

    write_json_atomic(json_path, d)

    return {
        "entry_id": entry_id,
        "pdb_id": pdb_id,
        "chain_id": chain_id,
        "json_path": str(json_path),
        "status": uniprot_mapping["status"],
        "sifts_chain_id": uniprot_mapping.get("sifts_chain_id", ""),
        "uniprot_accessions": ";".join(uniprot_mapping.get("uniprot_accessions", [])),
        "mapped_count": uniprot_mapping.get("mapped_count", 0),
        "chain_length": uniprot_mapping.get("chain_length", ""),
        "coverage": uniprot_mapping.get("coverage", 0),
        "identity": uniprot_mapping.get("identity", 0),
        "annotation_status": external_annotations["status"],
        "annotation_count": external_annotations["item_count"],
        "note": uniprot_mapping.get("reason", "")
    }


# ============================================================
# Main workflow
# ============================================================

def write_summary(path: Path, rows: List[Dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "entry_id",
        "pdb_id",
        "chain_id",
        "json_path",
        "status",
        "sifts_chain_id",
        "uniprot_accessions",
        "mapped_count",
        "chain_length",
        "coverage",
        "identity",
        "annotation_status",
        "annotation_count",
        "note"
    ]

    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main():
    p = argparse.ArgumentParser(
        description="Use SIFTS to map PDB residues to UniProt residues, fetch UniProt annotations, and write back to chain JSON."
    )

    p.add_argument("--manifest", default="example_outputs/manifest/manifest_skeleton.csv")
    p.add_argument("--json-root", default="example_outputs/json")
    p.add_argument("--cache-dir", default=".cache/sifts_uniprot")
    p.add_argument("--summary-out", default="example_outputs/qc/sifts_uniprot_annotation_summary.csv")
    p.add_argument("--limit-entries", type=int, default=0)
    p.add_argument("--include-all-uniprot-features", action="store_true")
    p.add_argument(
        "--feature-types",
        default="",
        help="Comma-separated UniProt feature types. Empty means built-in default useful features."
    )

    args = p.parse_args()

    manifest = Path(args.manifest)
    json_root = Path(args.json_root)
    cache_dir = Path(args.cache_dir)
    summary_out = Path(args.summary_out)

    if args.feature_types.strip():
        feature_types = {x.strip() for x in args.feature_types.split(",") if x.strip()}
    else:
        feature_types = DEFAULT_FEATURE_TYPES

    if manifest.exists():
        entries = get_entries_from_manifest(manifest, limit_entries=args.limit_entries)
    else:
        entries = get_entries_from_json_root(json_root, limit_entries=args.limit_entries)

    print("========== SIFTS + UniProt Annotation ==========")
    print(f"Manifest      : {manifest}")
    print(f"JSON root     : {json_root}")
    print(f"Cache dir     : {cache_dir}")
    print(f"Summary out   : {summary_out}")
    print(f"Entries       : {len(entries)}")
    print(f"Include all UniProt features: {args.include_all_uniprot_features}")
    print("================================================")

    summary_rows = []
    sifts_cache = {}

    for ei, entry_id in enumerate(entries, start=1):
        chain_jsons = find_chain_jsons(json_root, entry_id)

        if not chain_jsons:
            summary_rows.append({
                "entry_id": entry_id,
                "pdb_id": "",
                "chain_id": "",
                "json_path": "",
                "status": "failed",
                "sifts_chain_id": "",
                "uniprot_accessions": "",
                "mapped_count": 0,
                "chain_length": "",
                "coverage": 0,
                "identity": 0,
                "annotation_status": "",
                "annotation_count": 0,
                "note": "No chain JSON files found."
            })
            continue

        # Read the first chain JSON to get pdb_id
        first = read_json(chain_jsons[0])
        pdb_id = normalize_pdb_id(first.get("entry", {}).get("pdb_id", ""))

        print(f"[INFO] Entry {ei}/{len(entries)}: {entry_id}, pdb_id={pdb_id}, chains={len(chain_jsons)}")

        if not pdb_id or pdb_id.upper() == "NA":
            for jf in chain_jsons:
                d = read_json(jf)
                summary_rows.append({
                    "entry_id": entry_id,
                    "pdb_id": pdb_id,
                    "chain_id": d.get("chain", {}).get("chain_id", jf.stem),
                    "json_path": str(jf),
                    "status": "failed",
                    "sifts_chain_id": "",
                    "uniprot_accessions": "",
                    "mapped_count": 0,
                    "chain_length": d.get("chain", {}).get("length", ""),
                    "coverage": 0,
                    "identity": 0,
                    "annotation_status": "",
                    "annotation_count": 0,
                    "note": "Missing PDB ID; cannot query SIFTS."
                })
            continue

        try:
            if pdb_id not in sifts_cache:
                xml_gz = fetch_sifts_xml_gz(pdb_id, cache_dir)
                sifts_rows = parse_sifts_xml(xml_gz)
                sifts_cache[pdb_id] = build_sifts_index(sifts_rows)

            sifts_idx = sifts_cache[pdb_id]

            for ci, jf in enumerate(chain_jsons, start=1):
                print(f"  [CHAIN] {ci}/{len(chain_jsons)} {jf.name}")
                try:
                    row = update_chain_json(
                        json_path=jf,
                        sifts_idx=sifts_idx,
                        cache_dir=cache_dir,
                        feature_types=feature_types,
                        include_all_features=args.include_all_uniprot_features
                    )
                except Exception as e:
                    d = read_json(jf)
                    row = {
                        "entry_id": entry_id,
                        "pdb_id": pdb_id,
                        "chain_id": d.get("chain", {}).get("chain_id", jf.stem),
                        "json_path": str(jf),
                        "status": "failed",
                        "sifts_chain_id": "",
                        "uniprot_accessions": "",
                        "mapped_count": 0,
                        "chain_length": d.get("chain", {}).get("length", ""),
                        "coverage": 0,
                        "identity": 0,
                        "annotation_status": "",
                        "annotation_count": 0,
                        "note": repr(e)
                    }
                    print(f"    [ERROR] {repr(e)}")

                summary_rows.append(row)

        except Exception as e:
            print(f"[ERROR] Entry {entry_id} failed: {repr(e)}")
            for jf in chain_jsons:
                d = read_json(jf)
                summary_rows.append({
                    "entry_id": entry_id,
                    "pdb_id": pdb_id,
                    "chain_id": d.get("chain", {}).get("chain_id", jf.stem),
                    "json_path": str(jf),
                    "status": "failed",
                    "sifts_chain_id": "",
                    "uniprot_accessions": "",
                    "mapped_count": 0,
                    "chain_length": d.get("chain", {}).get("length", ""),
                    "coverage": 0,
                    "identity": 0,
                    "annotation_status": "",
                    "annotation_count": 0,
                    "note": repr(e)
                })

    write_summary(summary_out, summary_rows)

    c = Counter(row["status"] for row in summary_rows)
    a = Counter(row["annotation_status"] for row in summary_rows)

    print("========== Finished ==========")
    print(f"Chains processed : {len(summary_rows)}")
    print(f"Mapping status   : {dict(c)}")
    print(f"Annotation status: {dict(a)}")
    print(f"Summary written  : {summary_out}")
    print("==============================")


if __name__ == "__main__":
    main()