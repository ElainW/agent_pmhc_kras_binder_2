#!/bin/bash
# Generate hotspot-conditioned binder backbones against the prepped 9UV8 chain-B target.
# Default hotspots = peptide p4-p6 (B184-186), centered on the G12D mutant Asp(p5)=B185, per
# RFdiffusion's own "Picking Hotspots" guidance (3-6 hotspots; README.md lines 296-299)
# and the ctnnb1-15 single-hotspot analog in paper_output/design_hotspots.csv.
# inference.write_trajectory=False -- the per-step denoising trajectory isn't used downstream
# and silently consumed ~15GB for round 1's 1300 backbones (94MB of actual final PDBs).
#
# Usage: 02_run_rfdiffusion.sh <output_prefix> <num_designs> [hotspot_res, e.g. 'B184,B185,B186']
set -euo pipefail

OUTPUT_PREFIX="${1:?usage: 02_run_rfdiffusion.sh <output_prefix> <num_designs> [hotspot_res]}"
NUM_DESIGNS="${2:?usage: 02_run_rfdiffusion.sh <output_prefix> <num_designs> [hotspot_res]}"
HOTSPOT_RES="${3:-B184,B185,B186}"

mkdir -p "$(dirname "$OUTPUT_PREFIX")"

/venv/SE3nv/bin/python /workspace/RFdiffusion/scripts/run_inference.py \
  inference.output_prefix="$OUTPUT_PREFIX" \
  inference.input_pdb=/workspace/designs/prep/9UV8_target_chainB.pdb \
  inference.write_trajectory=False \
  'contigmap.contigs=[B1-189/0 70-100]' \
  "ppi.hotspot_res=[$HOTSPOT_RES]" \
  inference.num_designs="$NUM_DESIGNS" \
  denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0
