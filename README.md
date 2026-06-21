# pmhc_design scripts

Scripts and analysis code for the KRAS(G12D) pMHC miniprotein binder design campaign
(target: PDB 9UV8, HLA-A*11:01). See `/workspace/CLAUDE.md` and `/workspace/docs/` for
full project context, conventions, and the per-round design log.

Heavy outputs (RFdiffusion backbones, ProteinMPNN sequences, AF2/AF3 structures) live
under `/workspace/designs/<round>/` and are not committed here; this repo holds the
reusable pipeline code.

## Scripts

- `scripts/make_rfdiff_target.py` — builds the RFdiffusion docking target from `input/9UV8.pdb`:
  HLA-A*11:01 heavy chain residues 1-180 + KRAS(G12D) peptide residues 1-9, renumbered
  continuously (1-189) and relabeled to chain B (matches `contact_filter` convention: peptide
  = residues 181-189, p5/Asp(G12D) = residue 185). Drops beta-2-microglobulin (not part of the
  peptide-binding groove).
- `scripts/prepare_af3_jobs.py` — given a CSV of candidate binder sequences that passed the
  AF2 + ProteinMPNN-specificity filters, generates AlphaFold Server (AF3) batch-upload JSON
  files for manual submission: one on-target (9UV8, peptide VVGADGVGK) and one off-target
  (8I5E, peptide VVGAGGVGK) job per candidate, binder chain folded with no MSA, batched to the
  AF3 Server's 30-job/day quota. Drop the resulting structures back under
  `designs/<round>/af3_jobs/<batch>/results/` for downstream scoring.
