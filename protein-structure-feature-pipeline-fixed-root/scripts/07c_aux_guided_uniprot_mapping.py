#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import importlib.util
import json
import os
import re
from pathlib import Path
from collections import defaultdict, Counter


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def write_json_atomic(path, data):
    tmp = str(path) + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    os.replace(tmp, path)


def norm_seq(s):
    return re.sub(r"[^A-Za-z]", "", str(s or "")).upper().replace("*", "")


def parse_fasta_like(path):
    records = []
    h = None
    seq = []

    try:
        with open(path, encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if h is not None:
                        records.append((h, norm_seq("".join(seq))))
                    h = line[1:].strip()
                    seq = []
                else:
                    seq.append(line)
        if h is not None:
            records.append((h, norm_seq("".join(seq))))
    except Exception:
        return []

    return [(h, s) for h, s in records if h]


def load_aux_index(path):
    idx = defaultdict(list)
    p = Path(path)
    if not p.exists():
        return idx

    with open(p, newline="", encoding="utf-8", errors="ignore") as f:
        for r in csv.DictReader(f):
            pdb = r.get("pdb_id", "").strip().lower()
            if pdb:
                idx[pdb].append(r)

    return idx


def aux_files_for_pdb(aux_index, pdb_id):
    pdb_id = str(pdb_id or "").strip().lower()
    files = []
    for r in aux_index.get(pdb_id, []):
        p = Path(r.get("path", ""))
        if p.exists():
            files.append((r.get("suffix", ""), p))
    return files


def parse_torcsb(files, pdb_id):
    """
    Parse emdbMatch.py .toRCSB output.

    Header written by emdbMatch.py:
        >{t2c[t]}_{n}

    Common example:
        >6gyp_2_A_B

    Interpreted as:
        pdb_id = 6gyp
        entity = 2
        rcsb/template chain = A
        local chain = B
    """
    out = {}

    for suffix, p in files:
        if suffix != "toRCSB" and not str(p).endswith(".toRCSB"):
            continue

        for h, seq in parse_fasta_like(p):
            token = h.split()[0]
            parts = token.split("_")

            if len(parts) >= 4 and parts[0].lower() == str(pdb_id).lower():
                local_chain = parts[-1]
                rcsb_chain = parts[-2]
                entity_id = parts[-3]

                out[local_chain] = {
                    "rcsb_chain": rcsb_chain,
                    "entity_id": entity_id,
                    "header": h,
                    "sequence": seq,
                    "file": str(p),
                }

    return out


def parse_cif_struct_ref(files, pdb_id):
    """
    Parse mmCIF struct_ref / struct_ref_seq:
        RCSB chain -> UniProt accession
    """
    out = defaultdict(set)

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

    for suffix, p in files:
        if suffix != "cif" and not str(p).endswith(".cif"):
            continue

        try:
            d = MMCIF2Dict(str(p))
        except Exception:
            continue

        ref_ids = as_list(d.get("_struct_ref.id"))
        db_names = as_list(d.get("_struct_ref.db_name"))
        accessions = as_list(d.get("_struct_ref.pdbx_db_accession")) or as_list(d.get("_struct_ref.db_code"))

        ref_to_acc = {}

        for rid, dbn, acc in zip(ref_ids, db_names, accessions):
            dbn = str(dbn or "")
            acc = str(acc or "").strip()
            if not acc or acc in {".", "?"}:
                continue
            if "UNP" in dbn.upper() or "UNIPROT" in dbn.upper():
                ref_to_acc[str(rid)] = acc

        seq_ref_ids = as_list(d.get("_struct_ref_seq.ref_id"))
        strands = as_list(d.get("_struct_ref_seq.pdbx_strand_id"))
        pdb_codes = as_list(d.get("_struct_ref_seq.pdbx_PDB_id_code"))

        for rid, strand, code in zip(seq_ref_ids, strands, pdb_codes):
            rid = str(rid)
            acc = ref_to_acc.get(rid)
            if not acc:
                continue

            if code and str(code).lower() not in {".", "?", str(pdb_id).lower()}:
                continue

            for ch in re.split(r"[,; ]+", str(strand or "")):
                ch = ch.strip()
                if ch and ch not in {".", "?"}:
                    out[ch].add(acc)

    return out


def accessions_from_sifts_chain(sifts_idx, rcsb_chain):
    accs = set()

    rows_by_chain = sifts_idx.get("rows_by_chain", {})
    for r in rows_by_chain.get(rcsb_chain, []):
        acc = r.get("uniprot_acc") or r.get("accession") or r.get("dbAccessionId")
        if acc:
            accs.add(acc)

    return accs


def make_unreliable_tracks(d, reason):
    d["external_annotations"] = {
        "status": "unreliable_mapping",
        "source": "external annotations skipped because final UniProt mapping is unreliable",
        "coordinate_system": "chain_index_1_based",
        "items": [],
        "item_count": 0,
        "counts": {
            "uniprot_features": 0,
            "pfam": 0,
            "interpro": 0,
        },
        "reason": reason,
    }

    d.setdefault("tracks", {})
    for name in ["uniprot_features", "pfam", "interpro"]:
        d["tracks"][name] = {
            "status": "unreliable_mapping",
            "type": "features",
            "source": "not generated because final UniProt mapping is unreliable",
            "coordinate_system": "chain_index_1_based",
            "items": [],
            "item_count": 0,
        }

    d.setdefault("status", {})
    d["status"]["external_annotations"] = "unreliable_mapping"
    d["status"]["uniprot_features"] = "unreliable_mapping"
    d["status"]["pfam"] = "unreliable_mapping"
    d["status"]["interpro"] = "unreliable_mapping"


def build_annotations(mod, base, fmap, sifts_idx, cache_dir, feature_types):
    feature_by_acc = {}

    for acc in fmap.get("uniprot_accessions", []):
        try:
            uj = mod.fetch_uniprot_json(acc, cache_dir, base=base)
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


def iter_manifest_entries(manifest):
    seen = set()

    with open(manifest, newline="", encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f)

        for r in reader:
            entry_id = r.get("entry_id", "").strip()
            if not entry_id or entry_id in seen:
                continue
            seen.add(entry_id)

            has_structure = r.get("has_structure", "").lower()
            skeleton_status = r.get("skeleton_status", "").lower()

            if has_structure and has_structure not in {"yes", "true", "1"}:
                continue
            if skeleton_status and skeleton_status != "done":
                continue

            yield r


def run_one_chain(jf, args, mod, base, aux_index, sifts_cache, feature_types):
    d = read_json(jf)

    entry = d.get("entry", {})
    chain = d.get("chain", {})

    entry_id = entry.get("entry_id", jf.parent.parent.name)
    pdb_id = str(entry.get("pdb_id", "")).strip().lower()
    local_chain = str(chain.get("chain_id", jf.stem)).strip()
    chain_seq = norm_seq(chain.get("sequence", ""))

    if not pdb_id:
        return {
            "entry_id": entry_id,
            "pdb_id": "",
            "chain_id": local_chain,
            "json_path": str(jf),
            "status": "failed",
            "quality": "FAILED",
            "method": "aux_guided_mapping",
            "reason": "missing pdb_id",
        }

    files = aux_files_for_pdb(aux_index, pdb_id)
    torcsb = parse_torcsb(files, pdb_id)
    cif_chain_accs = parse_cif_struct_ref(files, pdb_id)

    rcsb_chain = ""
    torcsb_file = ""

    if local_chain in torcsb:
        rcsb_chain = torcsb[local_chain]["rcsb_chain"]
        torcsb_file = torcsb[local_chain]["file"]

    sifts_idx = mod.load_sifts_index(
        pdb_id,
        Path(args.cache_dir),
        base,
        sifts_cache,
    )

    candidate_accs = []
    candidate_sources = defaultdict(list)

    # 1. Most trusted: local -> RCSB chain from .toRCSB, then SIFTS accession from that RCSB chain
    if rcsb_chain:
        for acc in accessions_from_sifts_chain(sifts_idx, rcsb_chain):
            candidate_accs.append(acc)
            candidate_sources[acc].append(f"toRCSB_sifts_chain:{local_chain}->{rcsb_chain}")

        for acc in cif_chain_accs.get(rcsb_chain, []):
            candidate_accs.append(acc)
            candidate_sources[acc].append(f"toRCSB_cif_chain:{local_chain}->{rcsb_chain}")

    # 2. CIF accession for local chain directly
    for acc in cif_chain_accs.get(local_chain, []):
        candidate_accs.append(acc)
        candidate_sources[acc].append(f"cif_struct_ref_chain:{local_chain}")

    # 3. If no .toRCSB, fall back to SIFTS local chain accession
    if not candidate_accs:
        for acc in accessions_from_sifts_chain(sifts_idx, local_chain):
            candidate_accs.append(acc)
            candidate_sources[acc].append(f"sifts_local_chain:{local_chain}")

    # 4. Broad fallback: all accessions from SIFTS for this PDB
    if args.allow_all_sifts_fallback:
        for acc in sifts_idx.get("uniprot_accessions", []):
            candidate_accs.append(acc)
            candidate_sources[acc].append("sifts_all_accessions")

    seen = set()
    candidates = []
    for acc in candidate_accs:
        if acc and acc not in seen:
            seen.add(acc)
            candidates.append(acc)

    if not candidates:
        fmap = {
            "status": "failed",
            "quality": "FAILED",
            "method": "aux_guided_mapping",
            "final_mapping_method": "aux_guided_mapping",
            "source": "toRCSB/CIF/SIFTS",
            "uniprot_accessions": [],
            "mapped_count": 0,
            "chain_length": len(chain_seq),
            "coverage": 0.0,
            "identity": 0.0,
            "reason": "no candidate accession from toRCSB/CIF/SIFTS",
        }

        d["uniprot_mapping"] = fmap
        make_unreliable_tracks(d, fmap["reason"])

        if args.write:
            write_json_atomic(jf, d)

        return {
            "entry_id": entry_id,
            "pdb_id": pdb_id,
            "chain_id": local_chain,
            "json_path": str(jf),
            "status": fmap["status"],
            "quality": fmap["quality"],
            "method": fmap["final_mapping_method"],
            "rcsb_chain": rcsb_chain,
            "torcsb_file": torcsb_file,
            "candidate_count": 0,
            "coverage": 0.0,
            "identity": 0.0,
            "reason": fmap["reason"],
        }

    thresholds = (
        args.high_coverage,
        args.high_identity,
        args.ok_coverage,
        args.ok_identity,
        args.partial_coverage,
    )

    scored = []

    for acc in candidates[:args.max_candidates]:
        try:
            uni_seq = mod.get_uniprot_sequence(acc, Path(args.cache_dir), base=base)
        except Exception:
            try:
                uj = mod.fetch_uniprot_json(acc, Path(args.cache_dir), base=base)
                seq_obj = uj.get("sequence", {})
                if isinstance(seq_obj, dict):
                    uni_seq = norm_seq(seq_obj.get("value", ""))
                else:
                    uni_seq = ""
            except Exception:
                uni_seq = ""

        if not uni_seq:
            continue

        try:
            fmap = mod.alignment_mapping_object(
                d,
                acc,
                uni_seq,
                thresholds,
                source_note="aux_guided_toRCSB_CIF_SIFTS_alignment",
            )
        except TypeError:
            fmap = mod.alignment_mapping_object(
                d,
                acc,
                uni_seq,
                thresholds,
                "aux_guided_toRCSB_CIF_SIFTS_alignment",
            )

        fmap["candidate_sources"] = sorted(set(candidate_sources.get(acc, [])))
        scored.append(fmap)

    scored.sort(
        key=lambda x: (
            1 if x.get("quality") in {"HIGH", "ACCEPTABLE", "PARTIAL_USABLE"} else 0,
            float(x.get("identity", 0) or 0),
            float(x.get("coverage", 0) or 0),
            int(x.get("mapped_count", 0) or 0),
        ),
        reverse=True,
    )

    if not scored:
        fmap = {
            "status": "failed",
            "quality": "FAILED",
            "method": "aux_guided_mapping",
            "final_mapping_method": "aux_guided_mapping",
            "source": "toRCSB/CIF/SIFTS",
            "uniprot_accessions": [],
            "mapped_count": 0,
            "chain_length": len(chain_seq),
            "coverage": 0.0,
            "identity": 0.0,
            "reason": "candidate accessions found but no alignment could be scored",
        }

        d["uniprot_mapping"] = fmap
        make_unreliable_tracks(d, fmap["reason"])

        if args.write:
            write_json_atomic(jf, d)

        return {
            "entry_id": entry_id,
            "pdb_id": pdb_id,
            "chain_id": local_chain,
            "json_path": str(jf),
            "status": fmap["status"],
            "quality": fmap["quality"],
            "method": fmap["final_mapping_method"],
            "rcsb_chain": rcsb_chain,
            "torcsb_file": torcsb_file,
            "candidate_count": len(candidates),
            "coverage": 0.0,
            "identity": 0.0,
            "reason": fmap["reason"],
        }

    best = scored[0]
    trusted = best.get("quality") in {"HIGH", "ACCEPTABLE", "PARTIAL_USABLE"}

    best["source"] = "Aux-guided mapping using toRCSB/CIF/SIFTS-constrained UniProt candidates"
    best["final_mapping_method"] = "aux_guided_toRCSB_CIF_SIFTS_alignment"
    best["aux_guided"] = {
        "rcsb_chain_from_toRCSB": rcsb_chain,
        "torcsb_file": torcsb_file,
        "candidate_sources": best.get("candidate_sources", []),
        "candidate_count": len(candidates),
        "aux_file_count": len(files),
    }

    d["uniprot_mapping"] = best

    if trusted:
        uf, pfam, interpro, external = build_annotations(
            mod,
            base,
            best,
            sifts_idx,
            Path(args.cache_dir),
            feature_types,
        )

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
        d["status"]["uniprot_mapping"] = best.get("status", "")
        d["status"]["uniprot_mapping_quality"] = best.get("quality", "")
        d["status"]["uniprot_mapping_method"] = best.get("final_mapping_method", "")
        d["status"]["external_annotations"] = external.get("status", "")
        d["status"]["uniprot_features"] = uf.get("status", "")
        d["status"]["pfam"] = pfam.get("status", "")
        d["status"]["interpro"] = interpro.get("status", "")

    else:
        make_unreliable_tracks(d, best.get("reason", "aux-guided mapping is not reliable"))

    d.setdefault("logs", [])
    d["logs"].append({
        "step": "07c_aux_guided_uniprot_mapping",
        "status": "done" if trusted else "unreliable_mapping",
        "message": (
            f"quality={best.get('quality')}; "
            f"coverage={best.get('coverage')}; "
            f"identity={best.get('identity')}; "
            f"rcsb_chain={rcsb_chain}; "
            f"candidates={len(candidates)}"
        ),
    })

    if args.write:
        write_json_atomic(jf, d)

    return {
        "entry_id": entry_id,
        "pdb_id": pdb_id,
        "chain_id": local_chain,
        "json_path": str(jf),
        "status": best.get("status", ""),
        "quality": best.get("quality", ""),
        "method": best.get("final_mapping_method", ""),
        "rcsb_chain": rcsb_chain,
        "torcsb_file": torcsb_file,
        "candidate_count": len(candidates),
        "candidate_sources": "|".join(best.get("candidate_sources", [])),
        "uniprot_accessions": ";".join(best.get("uniprot_accessions", [])),
        "mapped_count": best.get("mapped_count", ""),
        "chain_length": best.get("chain_length", ""),
        "coverage": best.get("coverage", ""),
        "identity": best.get("identity", ""),
        "external_status": d.get("external_annotations", {}).get("status", ""),
        "pfam_status": d.get("tracks", {}).get("pfam", {}).get("status", ""),
        "interpro_status": d.get("tracks", {}).get("interpro", {}).get("status", ""),
        "reason": best.get("reason", ""),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json-root", default="example_outputs/json_backup_before_alignment_fallback")
    ap.add_argument("--manifest", default="example_outputs/manifest/manifest_skeleton.csv")
    ap.add_argument("--base-script", default="scripts/07_sifts_uniprot_annotations.py")
    ap.add_argument("--map-script", default="scripts/07b_sifts_alignment_fallback.py")
    ap.add_argument("--cache-dir", default=".cache/sifts_uniprot")
    ap.add_argument("--aux-index", default="example_outputs/qc/batch_aux_file_index.csv")
    ap.add_argument("--summary-out", default="example_outputs/qc/07c_aux_guided_mapping_summary.csv")
    ap.add_argument("--limit-entries", type=int, default=0)
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--allow-all-sifts-fallback", action="store_true")
    ap.add_argument("--max-candidates", type=int, default=80)

    ap.add_argument("--high-coverage", type=float, default=0.95)
    ap.add_argument("--high-identity", type=float, default=0.98)
    ap.add_argument("--ok-coverage", type=float, default=0.80)
    ap.add_argument("--ok-identity", type=float, default=0.90)
    ap.add_argument("--partial-coverage", type=float, default=0.50)

    args = ap.parse_args()

    mod = load_module(args.map_script, "map07b")
    base = mod.load_base_module(args.base_script)

    feature_types = getattr(mod, "DEFAULT_FEATURE_TYPES", None)
    if feature_types is None:
        feature_types = {
            "Initiator methionine", "Signal peptide", "Transit peptide", "Propeptide", "Chain", "Peptide",
            "Domain", "Repeat", "Region", "Coiled coil", "Motif", "Zinc finger",
            "Active site", "Binding site", "Metal binding", "Site",
            "Modified residue", "Lipidation", "Glycosylation", "Disulfide bond", "Cross-link",
            "Nucleotide binding", "DNA binding", "Calcium binding",
        }

    aux_index = load_aux_index(args.aux_index)
    sifts_cache = {}

    rows = []
    entry_count = 0
    chain_count = 0

    for r in iter_manifest_entries(args.manifest):
        entry_id = r.get("entry_id", "").strip()
        chain_dir = Path(args.json_root) / entry_id / "chains"

        if not chain_dir.is_dir():
            continue

        entry_count += 1
        if args.limit_entries and entry_count > args.limit_entries:
            break

        for jf in sorted(chain_dir.glob("*.json")):
            chain_count += 1
            print(f"[{chain_count}] {entry_id} {jf.name}")

            try:
                out = run_one_chain(
                    jf,
                    args,
                    mod,
                    base,
                    aux_index,
                    sifts_cache,
                    feature_types,
                )
            except Exception as e:
                out = {
                    "entry_id": entry_id,
                    "pdb_id": "",
                    "chain_id": jf.stem,
                    "json_path": str(jf),
                    "status": "error",
                    "quality": "FAILED",
                    "method": "aux_guided_toRCSB_CIF_SIFTS_alignment",
                    "reason": repr(e),
                }

            rows.append(out)

    out = Path(args.summary_out)
    out.parent.mkdir(parents=True, exist_ok=True)

    fields = []
    for r in rows:
        for k in r.keys():
            if k not in fields:
                fields.append(k)

    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    c = Counter(r.get("quality", "") for r in rows)
    m = Counter(r.get("method", "") for r in rows)
    ex = Counter(r.get("external_status", "") for r in rows)

    print("entries processed:", entry_count)
    print("chains processed:", len(rows))
    print("write:", args.write)
    print("quality:", dict(c))
    print("method:", dict(m))
    print("external:", dict(ex))
    print("summary:", out)


if __name__ == "__main__":
    main()