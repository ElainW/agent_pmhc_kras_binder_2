#!/bin/bash
# ProteinMPNN sequence design + 1 FastRelax cycle on RFdiffusion backbones.
# Binder = chain A (first chain, fully redesigned); target = chain B (last chain, fixed).
#
# Usage: 03_run_proteinmpnn_fastrelax.sh <pdbdir> <outpdbdir> [runlist]
set -euo pipefail

PDBDIR="${1:?usage: 03_run_proteinmpnn_fastrelax.sh <pdbdir> <outpdbdir> [runlist]}"
OUTPDBDIR="${2:?usage: 03_run_proteinmpnn_fastrelax.sh <pdbdir> <outpdbdir> [runlist]}"
RUNLIST="${3:-}"

mkdir -p "$OUTPDBDIR"
cd /workspace/dl_binder_design

EXTRA_ARGS=()
if [ -n "$RUNLIST" ]; then
    EXTRA_ARGS+=(-runlist "$RUNLIST")
fi

/venv/dl_binder_design/bin/python mpnn_fr/dl_interface_design.py \
  -pdbdir "$PDBDIR" \
  -outpdbdir "$OUTPDBDIR" \
  "${EXTRA_ARGS[@]}"
