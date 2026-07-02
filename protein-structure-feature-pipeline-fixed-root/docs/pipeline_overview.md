# Pipeline Overview

This workflow is organized around chain-level JSON records. Each entry is first registered in a manifest, then parsed into entry and chain skeletons. Downstream modules add local structural tracks, disorder/Q-score features, UniProt/SIFTS annotations, annotation groups, repair logs, and final QC summaries.

## Stages

1. **Manifest construction**: match structure entries with available FASTA files and record processing status.
2. **Skeleton construction**: parse PDB/mmCIF files and create entry-level and chain-level JSON files.
3. **Local track calculation**: compute structure-derived local tracks, including secondary structure, RSA/buried residues, hydropathy, approximate clashes, bond outliers, and cis-peptide flags.
4. **Disorder and Q-score integration**: run or parse optional external features, including IUPred2A disorder and ChimeraX/QScore outputs, and map residue-level values back to chain coordinates. These features are intentionally separated from Stage 02.
5. **Sequence feature QC**: validate whether local tracks have expected length, coordinate consistency, and usable status.
6. **SIFTS/UniProt annotation**: map PDB residues and chains to UniProt coordinates and convert UniProt features to chain coordinates.
7. **Fallback and repair**: use sequence alignment or auxiliary files to repair weak or failed mappings.
8. **Annotation grouping**: aggregate Pfam, InterPro, UniProt feature, and SIFTS-derived tracks.
9. **Final QC**: produce chain-level summary tables and problem-chain reports.
