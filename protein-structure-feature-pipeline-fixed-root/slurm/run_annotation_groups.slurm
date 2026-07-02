#!/bin/bash
#SBATCH --job-name=anno08
#SBATCH --output=logs/anno08_%j.out
#SBATCH --error=logs/anno08_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT:-$PWD}
PYTHON=${PYTHON:-python}
mkdir -p "${PROJECT_ROOT}/logs" "${PROJECT_ROOT}/example_outputs/qc"

${PYTHON} "${PROJECT_ROOT}/scripts/08_build_annotation_groups.py"   --json-root "${PROJECT_ROOT}/example_outputs/json"   --manifest "${PROJECT_ROOT}/example_outputs/manifest/manifest_skeleton.csv"   --summary-out "${PROJECT_ROOT}/example_outputs/qc/annotation_groups_summary.csv"
