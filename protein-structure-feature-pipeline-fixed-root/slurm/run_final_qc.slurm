#!/bin/bash
#SBATCH --job-name=qc10
#SBATCH --output=logs/qc10_%j.out
#SBATCH --error=logs/qc10_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G

set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT:-$PWD}
PYTHON=${PYTHON:-python}
mkdir -p "${PROJECT_ROOT}/logs" "${PROJECT_ROOT}/example_outputs/qc"

${PYTHON} "${PROJECT_ROOT}/scripts/10_qc_all_chain_json.py"   --json-root "${PROJECT_ROOT}/example_outputs/json"   --manifest "${PROJECT_ROOT}/example_outputs/manifest/manifest_skeleton.csv"   --summary-out "${PROJECT_ROOT}/example_outputs/qc/all_chain_qc_summary.csv"   --problem-out "${PROJECT_ROOT}/example_outputs/qc/all_chain_qc_problem_chains.csv"   --annotation-problem-out "${PROJECT_ROOT}/example_outputs/qc/annotation_problem_chains.csv"   --local-problem-out "${PROJECT_ROOT}/example_outputs/qc/local_parameter_problem_chains.csv"
