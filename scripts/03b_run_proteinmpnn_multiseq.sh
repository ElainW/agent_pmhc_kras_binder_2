#!/bin/bash
# ProteinMPNN sequence design WITHOUT FastRelax, N sequences/backbone (paper's supplement
# specifies 4-32 MPNN sequences/backbone -- see docs/03_design_log.md 2026-06-22 entries).
# dl_interface_design.py disallows seqs_per_struct>1 together with relax_cycles>0, and
# RFdiffusion's own README notes FastRelax didn't yield large in silico improvements anyway
# (it's mainly there for more shots-on-goal) -- so we drop it here in exchange for the
# paper's intended per-backbone sequence multiplicity. Run 07_esmfold_filter.py afterward
# to triage down to the best sequence/backbone before the expensive AF2-complex step.
#
# Binder = chain A (first chain, fully redesigned); target = chain B (last chain, fixed).
#
# Usage: 03b_run_proteinmpnn_multiseq.sh <pdbdir> <outpdbdir> [num_seqs_per_struct] [runlist]
set -euo pipefail

PDBDIR="${1:?usage: 03b_run_proteinmpnn_multiseq.sh <pdbdir> <outpdbdir> [num_seqs_per_struct] [runlist]}"
OUTPDBDIR="${2:?usage: 03b_run_proteinmpnn_multiseq.sh <pdbdir> <outpdbdir> [num_seqs_per_struct] [runlist]}"
NUM_SEQS="${3:-8}"
RUNLIST="${4:-}"

mkdir -p "$OUTPDBDIR"
cd /workspace/dl_binder_design

EXTRA_ARGS=()
if [ -n "$RUNLIST" ]; then
    EXTRA_ARGS+=(-runlist "$RUNLIST")
fi

/venv/dl_binder_design/bin/python mpnn_fr/dl_interface_design.py \
  -pdbdir "$PDBDIR" \
  -outpdbdir "$OUTPDBDIR" \
  -relax_cycles 0 \
  -seqs_per_struct "$NUM_SEQS" \
  "${EXTRA_ARGS[@]}"
