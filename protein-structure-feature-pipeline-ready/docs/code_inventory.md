# Code Inventory

| File | Role | Suggested stage |
|---|---|---|
| `00_build_manifest.py` | Match FASTA files to manifest entries | manifest |
| `01_build_skeleton.py` | Build entry/chain JSON skeletons from structure files | skeleton |
| `02_compute_local_tracks.py` | Compute DSSP/RSA/hydropathy/clash/bond/cis-peptide tracks | local features |
| `03_compute_disorder_qscore.py` | Add IUPred2A disorder and ChimeraX/QScore features | disorder / validation |
| `04_reparse_qscore_raw.py` | Reparse saved QScore raw text and write back to JSON | validation repair |
| `05_qc_sequence_features.py` | QC local sequence and structure tracks | feature QC |
| `07_sifts_uniprot_annotations.py` | Direct SIFTS and UniProt annotation mapping | annotation |
| `07b_sifts_alignment_fallback.py` | Alignment fallback for weak SIFTS residue mapping | annotation fallback |
| `07c_aux_guided_uniprot_mapping.py` | Auxiliary-file-guided UniProt mapping | annotation repair |
| `08_build_annotation_groups.py` | Build annotation group tracks | annotation groups |
| `09_repair_problem_chains_aux.py` | Repair problematic chains using auxiliary files | repair |
| `10_qc_all_chain_json.py` | Final all-chain QC reports | final QC |
| `11_torcsb_guided_mapping_repair.py` | Repair mapping through `.toRCSB` chain correspondence | repair |
| `emdbMatch_relaxed.py` | Generate/relax chain correspondence from CryoNet/emdbMatch-style inputs | helper |
| `run_emdbmatch_nosource_modify.py` | Run emdbMatch for selected problem entries and classify output usability | helper |
| diagnostic scripts | Inspect batch/web files and individual cases | diagnostics |
