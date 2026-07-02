# QC Rules

The pipeline uses conservative QC checks to identify feature and annotation problems.

## Local feature QC

- Track length should match chain sequence length when residue-level values are expected.
- Segments should be rebuildable from residue-level labels.
- Numeric tracks should contain finite values when marked as successful.
- Optional tools should be reported as `skipped` rather than silently missing.

## Annotation QC

- Direct SIFTS mapping is accepted when coverage and identity pass high-confidence thresholds.
- Alignment fallback is used when direct residue mapping is incomplete or inconsistent.
- Partial mappings are retained with explicit status instead of being treated as clean mappings.
- Repair scripts should write summaries and avoid silent overwrites unless `--write` is supplied.

## Final status categories

- `OK`: all required tracks or mappings are usable.
- `WARNING`: usable with incomplete or optional information.
- `PROBLEM`: requires fallback or repair.
- `FAILED`: no usable mapping or invalid chain record.
