#!/bin/bash
# Generate hotspot-conditioned binder backbones against the prepped 9UV8 chain-B target.
# Hotspots = peptide p4-p6 (B184-186), centered on the G12D mutant Asp(p5)=B185, per
# RFdiffusion's own "Picking Hotspots" guidance (3-6 hotspots; README.md lines 296-299)
# and the ctnnb1-15 single-hotspot analog in paper_output/design_hotspots.csv.
#
# Usage: 02_run_rfdiffusion.sh <output_prefix> <num_designs>
set -euo pipefail

OUTPUT_PREFIX="${1:?usage: 02_run_rfdiffusion.sh <output_prefix> <num_designs>}"
NUM_DESIGNS="${2:?usage: 02_run_rfdiffusion.sh <output_prefix> <num_designs>}"

mkdir -p "$(dirname "$OUTPUT_PREFIX")"

/venv/SE3nv/bin/python /workspace/RFdiffusion/scripts/run_inference.py \
  inference.output_prefix="$OUTPUT_PREFIX" \
  inference.input_pdb=/workspace/designs/prep/9UV8_target_chainB.pdb \
  'contigmap.contigs=[B1-189/0 70-100]' \
  'ppi.hotspot_res=[B184,B185,B186]' \
  inference.num_designs="$NUM_DESIGNS" \
  denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0
