# GitHub Update Commands

This archive is prepared as a flattened repository layout. Its contents should be copied directly into the root of the existing GitHub repository, not inside another nested project folder.

Recommended update workflow:

```bash
git clone https://github.com/Andy-ZY-Wong/protein-structure-feature-pipeline.git
cd protein-structure-feature-pipeline

# Unzip the fixed archive into a temporary folder first.
mkdir -p ../psfp_fixed
unzip /path/to/protein-structure-feature-pipeline-fixed-root.zip -d ../psfp_fixed

# Replace the repository content while preserving .git.
rsync -av --delete --exclude='.git' ../psfp_fixed/ ./

# Commit all additions, moves, edits, and deletions.
git status
git add -A
git commit -m "Clean repository layout and clarify pipeline stages"
git push
```

Main changes in this fixed version:

1. Flattened repository layout so that `README.md`, `scripts/`, `docs/`, `example_data/`, and `example_outputs/` live at the repository root.
2. Removed generated `__pycache__/` and `*.pyc` files.
3. Clarified Stage 02 vs Stage 03 responsibility:
   - Stage 02 now handles local structure-derived tracks only.
   - Stage 03 handles optional disorder and Q-score features.
4. Moved optional diagnostic and helper scripts from `scripts/` to `tools/`.
