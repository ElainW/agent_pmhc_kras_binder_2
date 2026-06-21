#!/bin/bash
# AF2 initial-guess filter (complex pae_interaction/binder_rmsd, or monomer binder pLDDT).
#
# taskset -c 0-15 works around the container's 1024 pids-cgroup limit: nproc reports 256
# cores, and jax/tensorflow size their CPU threadpools to the visible core count, which
# blows past the cgroup cap and crashes with pthread_create() failed (11). See
# pmhc_design/README.md "Environment notes".
#
# Usage: 04_run_af2_initial_guess.sh <pdbdir> <outdir> [complex|monomer]
set -euo pipefail

PDBDIR="${1:?usage: 04_run_af2_initial_guess.sh <pdbdir> <outdir> [complex|monomer]}"
OUTDIR="${2:?usage: 04_run_af2_initial_guess.sh <pdbdir> <outdir> [complex|monomer]}"
MODE="${3:-complex}"

mkdir -p "$OUTDIR"
cd /workspace/dl_binder_design

EXTRA_ARGS=()
if [ "$MODE" = "monomer" ]; then
    EXTRA_ARGS+=(-force_monomer)
fi

taskset -c 0-15 /venv/af2_binder_design/bin/python af2_initial_guess/predict.py \
  -pdbdir "$PDBDIR" \
  -outpdbdir "$OUTDIR" \
  -scorefilename "$OUTDIR/out.sc" \
  "${EXTRA_ARGS[@]}"
