# Data Schema

## Manifest fields

Typical manifest fields include:

- `entry_id`
- `cfdb_id`
- `emdb_id`
- `pdb_id`
- `structure_path`
- `has_structure`
- `fasta_path`
- `has_fasta`
- `json_dir`
- stage-specific status columns

## Chain JSON concepts

Each chain JSON stores:

- entry identifiers
- chain identifier
- sequence
- sequence length and sequence hash
- residue mapping between chain index and PDB residue numbering
- local feature tracks
- UniProt/SIFTS mapping
- external annotations
- annotation groups
- QC logs and stage statuses

## QC outputs

Final QC tables summarize chain-level local features, annotation features, missing tracks, unreliable mappings, and repair status.
