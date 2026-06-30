#!/bin/bash
#SBATCH --job-name=seq_local
#SBATCH --output=logs/seq_local_%A_%a.out
#SBATCH --error=logs/seq_local_%A_%a.err
#SBATCH --array=0-5%2
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=48:00:00

set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT:-$PWD}
PYTHON=${PYTHON:-python}
BATCH_ID=${SLURM_ARRAY_TASK_ID:-0}

mkdir -p "${PROJECT_ROOT}/logs" "${PROJECT_ROOT}/example_outputs/manifest"

${PYTHON} "${PROJECT_ROOT}/scripts/02_compute_local_tracks.py"   --manifest "${PROJECT_ROOT}/example_outputs/manifest/manifest_skeleton.csv"   --json-root "${PROJECT_ROOT}/example_outputs/json"   --manifest-out "${PROJECT_ROOT}/example_outputs/manifest/manifest_local_tracks_batch${BATCH_ID}.csv"   --batch-size 3000   --batch-index "${BATCH_ID}"
