#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
11_torcsb_guided_mapping_repair.py

Use emdbMatch .toRCSB files to repair UniProt mapping / external annotations
for chains whose local chain IDs no longer match original RCSB/PDB chain IDs.

Core logic:
    local JSON chain
    -> .toRCSB local_chain -> RCSB chain
    -> SIFTS/CIF gives UniProt candidate accession
    -> current JSON sequence vs UniProt sequence alignment
    -> if quality is good enough, write mapping + annotations back to JSON

Default behavior:
    - dry-run only; use --write to modify JSON
    - read source JSON from --json-root
    - write repaired JSON to --out-json-root
    - default output JSON root is:
      example_outputs/json_backup_before_alignment_fallback
    - do not overwrite existing HIGH / ACCEPTABLE mappings unless --overwrite-good is used
    - write HIGH / ACCEPTABLE by default
    - PARTIAL_USABLE can be written only with:
      --allow-partial-write
      and stricter secondary thresholds:
        coverage >= --partial-write-min-coverage
        identity >= --partial-write-min-identity
        mapped_count >= --partial-write-min-mapped-count
"""

import argparse
import csv
import importlib.util
import json
import os
import re
import shutil
from pathlib import Path
from collections import defaultdict, Counter


def load_module(path, name):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Cannot find module: {path}")
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def write_json_atomic(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = str(path) + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    os.replace(tmp, path)


def norm_seq(s):
    return re.sub(r"[^A-Za-z]", "", str(s or "")).upper().replace("*", "")


def get_chain_sequence(d):
    return norm_seq(d.get("chain", {}).get("sequence", ""))


def get_chain_id_from_json(d, fallback):
    return str(d.get("chain", {}).get("chain_id", fallback)).strip()


def parse_torcsb(path):
    """
    Parse .toRCSB produced by emdbMatch.py.

    Header example:
        >6nm5_1_3D_0

    Meaning:
        local chain = 0
        RCSB/template chain = 3D

    Return:
        local_chain -> {
            rcsb_chain,
            header,
            sequence
        }
    """
    path = Path(path)
    out = {}

    if not path.exists():
        return out

    current = None
    seq_lines = []

    def flush():
        if current is None:
            return
        local_chain = current["local_chain"]
        current["sequence"] = norm_seq("".join(seq_lines))
        out[local_chain] = dict(current)

    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                flush()

                h = line[1:].split()[0]
                parts = h.split("_")

                if len(parts) < 4:
                    current = None
                    seq_lines = []
                    continue

                local_chain = parts[-1]
                rcsb_chain = parts[-2]

                current = {
                    "header": h,
                    "local_chain": local_chain,
                    "rcsb_chain": rcsb_chain,
                }
                seq_lines = []
            else:
                if current is not None:
                    seq_lines.append(line)

    flush()
    return out


def copy_source_json_if_needed(src_json, dst_json):
    src_json = Path(src_json)
    dst_json = Path(dst_json)
    if not dst_json.exists():
        dst_json.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src_json, dst_json)


def json_path_for_chain(json_root, entry_id, chain_id):
    return Path(json_root) / entry_id / "chains" / f"{chain_id}.json"


def select_input_json(src_json, dst_json, args):
    """
    Decide which JSON to read for current old mapping state.

    If --read-from-out-json-if-exists is enabled and target JSON already exists,
    use target JSON as the input state. This is useful when writing into an
    existing test/backup JSON tree.

    Otherwise use source JSON from --json-root.
    """
    src_json = Path(src_json)
    dst_json = Path(dst_json)

    if args.read_from_out_json_if_exists and dst_json.exists():
        return dst_json

    return src_json


def candidate_accessions_from_sifts_chain(sifts_idx, rcsb_chain):
    """
    Try to get UniProt accession candidates from the loaded SIFTS index.

    This function is intentionally tolerant because previous scripts may use
    slightly different keys in sifts_idx rows.
    """
    accs = set()

    rows_by_chain = sifts_idx.get("rows_by_chain", {})
    rows = rows_by_chain.get(rcsb_chain, [])

    for r in rows:
        for k in [
            "uniprot_acc",
            "accession",
            "dbAccessionId",
            "UniProt",
            "SP_PRIMARY",
            "uniprot",
            "unp_acc",
        ]:
            acc = r.get(k)
            if acc and str(acc).strip() not in {".", "?", "None"}:
                accs.add(str(acc).strip())

    return accs


def parse_cif_struct_ref_for_chain(cif_path, pdb_id):
    """
    Parse mmCIF struct_ref / struct_ref_seq and return:
        RCSB chain -> set(UniProt accession)

    This is a backup source when SIFTS rows_by_chain does not expose accession.
    """
    out = defaultdict(set)
    cif_path = Path(cif_path)

    if not cif_path.exists():
        return out

    try:
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    except Exception:
        return out

    def as_list(x):
        if x is None:
            return []
        if isinstance(x, list):
            return x
        return [x]

    try:
        d = MMCIF2Dict(str(cif_path))
    except Exception:
        return out

    ref_ids = as_list(d.get("_struct_ref.id"))
    db_names = as_list(d.get("_struct_ref.db_name"))
    accessions = (
        as_list(d.get("_struct_ref.pdbx_db_accession"))
        or as_list(d.get("_struct_ref.db_code"))
    )

    ref_to_acc = {}

    for rid, dbn, acc in zip(ref_ids, db_names, accessions):
        rid = str(rid)
        dbn = str(dbn or "")
        acc = str(acc or "").strip()

        if not acc or acc in {".", "?"}:
            continue

        if "UNP" in dbn.upper() or "UNIPROT" in dbn.upper():
            ref_to_acc[rid] = acc

    seq_ref_ids = as_list(d.get("_struct_ref_seq.ref_id"))
    strands = as_list(d.get("_struct_ref_seq.pdbx_strand_id"))
    pdb_codes = as_list(d.get("_struct_ref_seq.pdbx_PDB_id_code"))

    for rid, strand, code in zip(seq_ref_ids, strands, pdb_codes):
        rid = str(rid)
        acc = ref_to_acc.get(rid)

        if not acc:
            continue

        code = str(code or "").strip().lower()
        if code not in {"", ".", "?", str(pdb_id).lower()}:
            continue

        for ch in re.split(r"[,; ]+", str(strand or "")):
            ch = ch.strip()
            if ch and ch not in {".", "?"}:
                out[ch].add(acc)

    return out


def get_uniprot_sequence_compat(mod, acc, cache_dir, base):
    """
    Compatible wrapper around previous 07b functions.
    """
    cache_dir = Path(cache_dir)

    if hasattr(mod, "get_uniprot_sequence"):
        return norm_seq(mod.get_uniprot_sequence(acc, cache_dir, base=base))

    if hasattr(mod, "fetch_uniprot_json"):
        uj = mod.fetch_uniprot_json(acc, cache_dir, base=base)
        seq_obj = uj.get("sequence", {})
        if isinstance(seq_obj, dict):
            return norm_seq(seq_obj.get("value", ""))
        return ""

    raise RuntimeError(
        "The map script has neither get_uniprot_sequence nor fetch_uniprot_json"
    )


def alignment_mapping_object_compat(mod, d, acc, uni_seq, thresholds, source_note):
    """
    Compatible call for alignment_mapping_object, because earlier versions may
    use either positional or keyword source_note.
    """
    try:
        return mod.alignment_mapping_object(
            d,
            acc,
            uni_seq,
            thresholds,
            source_note=source_note,
        )
    except TypeError:
        return mod.alignment_mapping_object(
            d,
            acc,
            uni_seq,
            thresholds,
            source_note,
        )


def quality_rank(q):
    order = {
        "HIGH": 5,
        "ACCEPTABLE": 4,
        "PARTIAL_USABLE": 3,
        "PROBLEM": 2,
        "FAILED": 1,
        "": 0,
        None: 0,
    }
    return order.get(q, 0)


def should_try_repair(old_quality, overwrite_good):
    if overwrite_good:
        return True
    return old_quality not in {"HIGH", "ACCEPTABLE"}


def is_write_quality(q, allow_partial, fmap, args):
    """
    Decide whether a new mapping result is allowed to be written.

    HIGH / ACCEPTABLE:
        always writable.

    PARTIAL_USABLE:
        writable only when explicitly allowed and passes stricter secondary
        thresholds. This avoids writing short or low-confidence partial matches.
    """
    if q in {"HIGH", "ACCEPTABLE"}:
        return True, ""

    if q != "PARTIAL_USABLE":
        return False, f"new quality {q} is not allowed for writing"

    if not allow_partial:
        return False, "PARTIAL_USABLE is not allowed without --allow-partial-write"

    try:
        coverage = float(fmap.get("coverage", 0) or 0)
    except Exception:
        coverage = 0.0

    try:
        identity = float(fmap.get("identity", 0) or 0)
    except Exception:
        identity = 0.0

    try:
        mapped_count = int(fmap.get("mapped_count", 0) or 0)
    except Exception:
        mapped_count = 0

    if coverage < args.partial_write_min_coverage:
        return False, (
            f"PARTIAL_USABLE coverage {coverage:.4f} < "
            f"{args.partial_write_min_coverage:.4f}"
        )

    if identity < args.partial_write_min_identity:
        return False, (
            f"PARTIAL_USABLE identity {identity:.4f} < "
            f"{args.partial_write_min_identity:.4f}"
        )

    if mapped_count < args.partial_write_min_mapped_count:
        return False, (
            f"PARTIAL_USABLE mapped_count {mapped_count} < "
            f"{args.partial_write_min_mapped_count}"
        )

    return True, ""


def make_feature_types(mod):
    if hasattr(mod, "DEFAULT_FEATURE_TYPES"):
        return getattr(mod, "DEFAULT_FEATURE_TYPES")

    return {
        "Initiator methionine",
        "Signal peptide",
        "Transit peptide",
        "Propeptide",
        "Chain",
        "Peptide",
        "Topological domain",
        "Transmembrane",
        "Intramembrane",
        "Domain",
        "Repeat",
        "Region",
        "Coiled coil",
        "Motif",
        "Compositional bias",
        "Zinc finger",
        "DNA binding",
        "Nucleotide binding",
        "Calcium binding",
        "Metal binding",
        "Binding site",
        "Active site",
        "Site",
        "Modified residue",
        "Lipidation",
        "Glycosylation",
        "Disulfide bond",
        "Cross-link",
    }


def build_annotations(mod, base, fmap, sifts_idx, cache_dir, feature_types):
    """
    Build UniProt, Pfam, InterPro and merged external annotations using earlier
    07b/07 functions.
    """
    feature_by_acc = {}

    for acc in fmap.get("uniprot_accessions", []):
        try:
            uj = mod.fetch_uniprot_json(acc, Path(cache_dir), base=base)
            feature_by_acc[acc] = mod.extract_uniprot_features(
                uj,
                feature_types,
                include_all=False,
            )
        except Exception:
            feature_by_acc[acc] = []

    uf = mod.convert_uniprot_features(fmap, feature_by_acc)
    pfam = mod.convert_sifts_xrefs(fmap, sifts_idx.get("xref_rows", []), "Pfam")
    interpro = mod.convert_sifts_xrefs(fmap, sifts_idx.get("xref_rows", []), "InterPro")
    external = mod.merge_external(uf, pfam, interpro)

    return uf, pfam, interpro, external


def write_mapping_and_annotations(d, fmap, uf, pfam, interpro, external):
    d["uniprot_mapping"] = fmap

    d["external_annotations"] = external

    d.setdefault("tracks", {})

    d["tracks"]["uniprot_features"] = {
        "status": uf.get("status", ""),
        "type": "features",
        "source": uf.get("source", ""),
        "coordinate_system": "chain_index_1_based",
        "items": uf.get("items", []),
        "item_count": uf.get("item_count", len(uf.get("items", []))),
    }

    d["tracks"]["pfam"] = {
        "status": pfam.get("status", ""),
        "type": "features",
        "source": pfam.get("source", ""),
        "coordinate_system": "chain_index_1_based",
        "items": pfam.get("items", []),
        "item_count": pfam.get("item_count", len(pfam.get("items", []))),
    }

    d["tracks"]["interpro"] = {
        "status": interpro.get("status", ""),
        "type": "features",
        "source": interpro.get("source", ""),
        "coordinate_system": "chain_index_1_based",
        "items": interpro.get("items", []),
        "item_count": interpro.get("item_count", len(interpro.get("items", []))),
    }

    d.setdefault("status", {})
    d["status"]["uniprot_mapping"] = fmap.get("status", "")
    d["status"]["uniprot_mapping_quality"] = fmap.get("quality", "")
    d["status"]["uniprot_mapping_method"] = fmap.get(
        "final_mapping_method",
        fmap.get("method", ""),
    )
    d["status"]["external_annotations"] = external.get("status", "")
    d["status"]["uniprot_features"] = uf.get("status", "")
    d["status"]["pfam"] = pfam.get("status", "")
    d["status"]["interpro"] = interpro.get("status", "")

    d.setdefault("logs", [])
    d["logs"].append(
        {
            "step": "11_torcsb_guided_mapping_repair",
            "status": "done",
            "message": (
                f"quality={fmap.get('quality')}; "
                f"coverage={fmap.get('coverage')}; "
                f"identity={fmap.get('identity')}; "
                f"rcsb_chain={fmap.get('torcsb_guided', {}).get('rcsb_chain', '')}; "
                f"accessions={','.join(fmap.get('uniprot_accessions', []))}"
            ),
        }
    )


def iter_repairable_rows(path, limit=0):
    usable_status = {
        "done_full_toRCSB",
        "partial_toRCSB_usable",
        "partial_toRCSB_usable_full_overlap_returncode_nonzero",
    }

    n = 0
    with open(path, newline="", encoding="utf-8", errors="ignore") as f:
        for r in csv.DictReader(f):
            if r.get("status") not in usable_status:
                continue

            n += 1
            if limit and n > limit:
                break

            yield r


def row_to_cif_path(r):
    for k in ["local_cif", "batch_cif"]:
        p = r.get(k, "")
        if p and Path(p).exists():
            return p
    return ""


def process_chain(
    *,
    entry_row,
    chain_id,
    torcsb_item,
    args,
    mod,
    base,
    sifts_cache,
    feature_types,
):
    entry_id = entry_row["entry_id"]
    pdb_id = entry_row["pdb_id"].lower()

    src_json = json_path_for_chain(args.json_root, entry_id, chain_id)
    dst_json = json_path_for_chain(args.out_json_root, entry_id, chain_id)

    if not src_json.exists() and not dst_json.exists():
        return {
            "entry_id": entry_id,
            "pdb_id": pdb_id,
            "local_chain": chain_id,
            "rcsb_chain": torcsb_item.get("rcsb_chain", ""),
            "json_path": str(src_json),
            "out_json_path": str(dst_json),
            "can_write": "no",
            "write": "no",
            "reason": "missing source json and target json",
        }

    input_json = select_input_json(src_json, dst_json, args)

    if not input_json.exists():
        return {
            "entry_id": entry_id,
            "pdb_id": pdb_id,
            "local_chain": chain_id,
            "rcsb_chain": torcsb_item.get("rcsb_chain", ""),
            "json_path": str(src_json),
            "out_json_path": str(dst_json),
            "can_write": "no",
            "write": "no",
            "reason": "missing input json",
        }

    d = read_json(input_json)

    json_chain_id = get_chain_id_from_json(d, chain_id)
    json_seq = get_chain_sequence(d)
    tor_seq = norm_seq(torcsb_item.get("sequence", ""))

    if args.require_torcsb_seq_match:
        if not tor_seq:
            return {
                "entry_id": entry_id,
                "pdb_id": pdb_id,
                "local_chain": chain_id,
                "rcsb_chain": torcsb_item.get("rcsb_chain", ""),
                "json_path": str(input_json),
                "out_json_path": str(dst_json),
                "can_write": "no",
                "write": "no",
                "reason": "empty toRCSB sequence",
            }

        if json_seq != tor_seq:
            return {
                "entry_id": entry_id,
                "pdb_id": pdb_id,
                "local_chain": chain_id,
                "rcsb_chain": torcsb_item.get("rcsb_chain", ""),
                "json_path": str(input_json),
                "out_json_path": str(dst_json),
                "json_len": len(json_seq),
                "torcsb_len": len(tor_seq),
                "can_write": "no",
                "write": "no",
                "reason": "json sequence and toRCSB sequence mismatch",
            }

    old_map = d.get("uniprot_mapping", {})
    old_quality = old_map.get("quality", "")
    old_status = old_map.get("status", "")
    old_acc = ";".join(old_map.get("uniprot_accessions", []) or [])

    if not should_try_repair(old_quality, args.overwrite_good):
        return {
            "entry_id": entry_id,
            "pdb_id": pdb_id,
            "local_chain": chain_id,
            "rcsb_chain": torcsb_item.get("rcsb_chain", ""),
            "old_quality": old_quality,
            "old_status": old_status,
            "old_accessions": old_acc,
            "can_write": "no",
            "write": "no",
            "reason": "old mapping is already HIGH/ACCEPTABLE",
            "json_path": str(input_json),
            "out_json_path": str(dst_json),
        }

    cache_dir = Path(args.cache_dir)

    sifts_idx = mod.load_sifts_index(
        pdb_id,
        cache_dir,
        base,
        sifts_cache,
    )

    rcsb_chain = torcsb_item["rcsb_chain"]

    candidate_accs = []
    candidate_sources = defaultdict(list)

    for acc in candidate_accessions_from_sifts_chain(sifts_idx, rcsb_chain):
        candidate_accs.append(acc)
        candidate_sources[acc].append(f"toRCSB_sifts_chain:{chain_id}->{rcsb_chain}")

    cif_path = row_to_cif_path(entry_row)
    cif_chain_accs = parse_cif_struct_ref_for_chain(cif_path, pdb_id) if cif_path else {}

    for acc in cif_chain_accs.get(rcsb_chain, []):
        candidate_accs.append(acc)
        candidate_sources[acc].append(f"toRCSB_cif_chain:{chain_id}->{rcsb_chain}")

    seen = set()
    candidates = []
    for acc in candidate_accs:
        if acc and acc not in seen:
            seen.add(acc)
            candidates.append(acc)

    if not candidates:
        return {
            "entry_id": entry_id,
            "pdb_id": pdb_id,
            "local_chain": chain_id,
            "rcsb_chain": rcsb_chain,
            "old_quality": old_quality,
            "old_status": old_status,
            "old_accessions": old_acc,
            "candidate_count": 0,
            "can_write": "no",
            "write": "no",
            "reason": "no UniProt candidate from toRCSB-guided SIFTS/CIF lookup",
            "json_path": str(input_json),
            "out_json_path": str(dst_json),
            "toRCSB_path": entry_row.get("out_toRCSB", ""),
        }

    thresholds = (
        args.high_coverage,
        args.high_identity,
        args.ok_coverage,
        args.ok_identity,
        args.partial_coverage,
    )

    scored = []

    for acc in candidates[: args.max_candidates]:
        try:
            uni_seq = get_uniprot_sequence_compat(mod, acc, cache_dir, base)
        except Exception:
            uni_seq = ""

        if not uni_seq:
            continue

        try:
            fmap = alignment_mapping_object_compat(
                mod,
                d,
                acc,
                uni_seq,
                thresholds,
                source_note="toRCSB_guided_SIFTS_CIF_alignment",
            )
        except Exception as e:
            scored.append(
                {
                    "status": "failed",
                    "quality": "FAILED",
                    "uniprot_accessions": [acc],
                    "coverage": 0.0,
                    "identity": 0.0,
                    "mapped_count": 0,
                    "reason": f"alignment failed: {repr(e)}",
                }
            )
            continue

        fmap["candidate_sources"] = sorted(set(candidate_sources.get(acc, [])))
        scored.append(fmap)

    scored.sort(
        key=lambda x: (
            quality_rank(x.get("quality")),
            float(x.get("identity", 0) or 0),
            float(x.get("coverage", 0) or 0),
            int(x.get("mapped_count", 0) or 0),
        ),
        reverse=True,
    )

    if not scored:
        return {
            "entry_id": entry_id,
            "pdb_id": pdb_id,
            "local_chain": chain_id,
            "rcsb_chain": rcsb_chain,
            "old_quality": old_quality,
            "old_status": old_status,
            "old_accessions": old_acc,
            "candidate_count": len(candidates),
            "can_write": "no",
            "write": "no",
            "reason": "no candidate could be aligned",
            "json_path": str(input_json),
            "out_json_path": str(dst_json),
            "toRCSB_path": entry_row.get("out_toRCSB", ""),
        }

    best = scored[0]
    new_quality = best.get("quality", "")
    new_status = best.get("status", "")
    new_acc = ";".join(best.get("uniprot_accessions", []) or [])
    new_cov = best.get("coverage", "")
    new_id = best.get("identity", "")
    new_mapped = best.get("mapped_count", "")

    best["source"] = "toRCSB-guided UniProt mapping repair"
    best["final_mapping_method"] = "toRCSB_guided_SIFTS_CIF_alignment"
    best["torcsb_guided"] = {
        "local_chain": chain_id,
        "json_chain_id": json_chain_id,
        "rcsb_chain": rcsb_chain,
        "toRCSB_header": torcsb_item.get("header", ""),
        "toRCSB_path": entry_row.get("out_toRCSB", ""),
        "candidate_count": len(candidates),
        "candidate_sources": best.get("candidate_sources", []),
        "input_json": str(input_json),
        "output_json": str(dst_json),
    }

    can_write, reason = is_write_quality(
        new_quality,
        args.allow_partial_write,
        best,
        args,
    )

    if can_write and quality_rank(new_quality) < quality_rank(old_quality):
        can_write = False
        reason = f"new quality {new_quality} is worse than old quality {old_quality}"

    row = {
        "entry_id": entry_id,
        "pdb_id": pdb_id,
        "local_chain": chain_id,
        "rcsb_chain": rcsb_chain,
        "old_status": old_status,
        "old_quality": old_quality,
        "old_accessions": old_acc,
        "new_status": new_status,
        "new_quality": new_quality,
        "new_accessions": new_acc,
        "coverage": new_cov,
        "identity": new_id,
        "mapped_count": new_mapped,
        "candidate_count": len(candidates),
        "candidate_sources": "|".join(best.get("candidate_sources", [])),
        "can_write": "yes" if can_write else "no",
        "write": "no",
        "reason": reason,
        "json_path": str(input_json),
        "src_json_path": str(src_json),
        "out_json_path": str(dst_json),
        "toRCSB_path": entry_row.get("out_toRCSB", ""),
        "toRCSB_header": torcsb_item.get("header", ""),
    }

    if can_write and args.write:
        target_json = src_json if Path(args.out_json_root) == Path(args.json_root) else dst_json

        if target_json != src_json:
            copy_source_json_if_needed(input_json, target_json)
            d = read_json(target_json)

        uf, pfam, interpro, external = build_annotations(
            mod,
            base,
            best,
            sifts_idx,
            cache_dir,
            feature_types,
        )

        write_mapping_and_annotations(d, best, uf, pfam, interpro, external)
        write_json_atomic(target_json, d)

        row["write"] = "yes"
        row["out_json_path"] = str(target_json)

    return row


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--repairable-csv",
        default="example_outputs/qc/emdbmatch_repairable_all.csv",
    )

    ap.add_argument(
        "--json-root",
        default="example_outputs/json",
    )

    ap.add_argument(
        "--out-json-root",
        default="example_outputs/json_backup_before_alignment_fallback",
    )

    ap.add_argument(
        "--base-script",
        default="scripts/07_sifts_uniprot_annotations.py",
    )

    ap.add_argument(
        "--map-script",
        default="scripts/07b_sifts_alignment_fallback.py",
    )

    ap.add_argument(
        "--cache-dir",
        default=".cache/sifts_uniprot",
    )

    ap.add_argument(
        "--summary-out",
        default="example_outputs/qc/torcsb_guided_repair_dryrun.csv",
    )

    ap.add_argument("--limit-entries", type=int, default=0)
    ap.add_argument("--limit-chains", type=int, default=0)

    ap.add_argument("--write", action="store_true")
    ap.add_argument("--overwrite-good", action="store_true")
    ap.add_argument("--allow-partial-write", action="store_true")

    ap.add_argument(
        "--read-from-out-json-if-exists",
        action="store_true",
        default=True,
    )

    ap.add_argument(
        "--no-read-from-out-json-if-exists",
        dest="read_from_out_json_if_exists",
        action="store_false",
    )

    ap.add_argument(
        "--no-require-torcsb-seq-match",
        dest="require_torcsb_seq_match",
        action="store_false",
    )
    ap.set_defaults(require_torcsb_seq_match=True)

    ap.add_argument("--max-candidates", type=int, default=80)

    ap.add_argument("--high-coverage", type=float, default=0.95)
    ap.add_argument("--high-identity", type=float, default=0.98)
    ap.add_argument("--ok-coverage", type=float, default=0.80)
    ap.add_argument("--ok-identity", type=float, default=0.90)
    ap.add_argument("--partial-coverage", type=float, default=0.50)

    ap.add_argument("--partial-write-min-coverage", type=float, default=0.70)
    ap.add_argument("--partial-write-min-identity", type=float, default=0.95)
    ap.add_argument("--partial-write-min-mapped-count", type=int, default=50)

    args = ap.parse_args()

    mod = load_module(args.map_script, "map07b")
    base = mod.load_base_module(args.base_script)

    feature_types = make_feature_types(mod)
    sifts_cache = {}

    rows = []
    counters = Counter()

    entry_seen = 0
    chain_seen = 0

    for entry_row in iter_repairable_rows(args.repairable_csv, limit=args.limit_entries):
        entry_seen += 1

        entry_id = entry_row["entry_id"]
        torcsb_path = entry_row.get("out_toRCSB", "")

        tor_map = parse_torcsb(torcsb_path)
        if not tor_map:
            counters["entry_no_torcsb"] += 1
            continue

        json_chain_dir = Path(args.json_root) / entry_id / "chains"
        out_json_chain_dir = Path(args.out_json_root) / entry_id / "chains"

        if json_chain_dir.is_dir():
            json_chains = sorted(
                [p.stem for p in json_chain_dir.glob("*.json")],
                key=lambda x: (len(x), x),
            )
        elif out_json_chain_dir.is_dir():
            json_chains = sorted(
                [p.stem for p in out_json_chain_dir.glob("*.json")],
                key=lambda x: (len(x), x),
            )
        else:
            counters["entry_missing_json_chain_dir"] += 1
            continue

        for chain_id in json_chains:
            if args.limit_chains and chain_seen >= args.limit_chains:
                break

            if chain_id not in tor_map:
                counters["chain_not_covered_by_torcsb"] += 1
                continue

            chain_seen += 1

            try:
                row = process_chain(
                    entry_row=entry_row,
                    chain_id=chain_id,
                    torcsb_item=tor_map[chain_id],
                    args=args,
                    mod=mod,
                    base=base,
                    sifts_cache=sifts_cache,
                    feature_types=feature_types,
                )
            except Exception as e:
                row = {
                    "entry_id": entry_id,
                    "pdb_id": entry_row.get("pdb_id", ""),
                    "local_chain": chain_id,
                    "rcsb_chain": tor_map.get(chain_id, {}).get("rcsb_chain", ""),
                    "can_write": "no",
                    "write": "no",
                    "reason": f"exception: {repr(e)}",
                    "json_path": str(json_path_for_chain(args.json_root, entry_id, chain_id)),
                    "out_json_path": str(json_path_for_chain(args.out_json_root, entry_id, chain_id)),
                    "toRCSB_path": torcsb_path,
                }

            rows.append(row)

            counters["chains_processed"] += 1
            counters[f"can_write::{row.get('can_write', '')}"] += 1
            counters[f"write::{row.get('write', '')}"] += 1
            counters[f"new_quality::{row.get('new_quality', '')}"] += 1
            counters[f"reason::{row.get('reason', '')}"] += 1

        if args.limit_chains and chain_seen >= args.limit_chains:
            break

    out = Path(args.summary_out)
    out.parent.mkdir(parents=True, exist_ok=True)

    fields = []
    for r in rows:
        for k in r.keys():
            if k not in fields:
                fields.append(k)

    with open(out, "w", newline="", encoding="utf-8") as f:
        if fields:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)

    print("entries seen:", entry_seen)
    print("chains seen:", chain_seen)
    print("write mode:", args.write)
    print("json root:", args.json_root)
    print("out json root:", args.out_json_root)
    print("summary:", out)
    print("counters:")
    for k, v in counters.most_common(80):
        print(k, v)


if __name__ == "__main__":
    main()