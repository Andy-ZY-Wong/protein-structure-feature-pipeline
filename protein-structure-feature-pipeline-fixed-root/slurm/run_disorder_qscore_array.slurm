#!/bin/bash
#SBATCH --job-name=dq_all
#SBATCH --output=logs/dq_all_%A_%a.out
#SBATCH --error=logs/dq_all_%A_%a.err
#SBATCH --array=0-5%2
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=72:00:00

set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT:-$PWD}
PYTHON=${PYTHON:-python}
BATCH_ID=${SLURM_ARRAY_TASK_ID:-0}
IUPRED_CMD=${IUPRED_CMD:-""}
CHIMERAX_CMD=${CHIMERAX_CMD:-"chimerax"}

mkdir -p "${PROJECT_ROOT}/logs" "${PROJECT_ROOT}/example_outputs/qc/qscore_raw_all" "${PROJECT_ROOT}/example_outputs/manifest"

${PYTHON} "${PROJECT_ROOT}/scripts/03_compute_disorder_qscore.py"   --manifest "${PROJECT_ROOT}/example_outputs/manifest/manifest_skeleton.csv"   --json-root "${PROJECT_ROOT}/example_outputs/json"   --manifest-out "${PROJECT_ROOT}/example_outputs/manifest/manifest_disorder_qscore_batch${BATCH_ID}.csv"   --batch-size 3000   --batch-index "${BATCH_ID}"   --iupred-cmd "${IUPRED_CMD}"   --map-table "${PROJECT_ROOT}/example_data/raw/maps/map_paths.csv"   --chimerax-cmd "${CHIMERAX_CMD}"   --keep-qscore-raw "${PROJECT_ROOT}/example_outputs/qc/qscore_raw_all"   --qscore-timeout 3600   --iupred-timeout 300
