#!/bin/bash
# ProteinMPNN sequence design + 1 FastRelax cycle on RFdiffusion backbones.
# Binder = chain A (first chain, fully redesigned); target = chain B (last chain, fixed).
#
# Usage: 03_run_proteinmpnn_fastrelax.sh <pdbdir> <outpdbdir>
set -euo pipefail

PDBDIR="${1:?usage: 03_run_proteinmpnn_fastrelax.sh <pdbdir> <outpdbdir>}"
OUTPDBDIR="${2:?usage: 03_run_proteinmpnn_fastrelax.sh <pdbdir> <outpdbdir>}"

mkdir -p "$OUTPDBDIR"
cd /workspace/dl_binder_design
/venv/dl_binder_design/bin/python mpnn_fr/dl_interface_design.py \
  -pdbdir "$PDBDIR" \
  -outpdbdir "$OUTPDBDIR"
