# Design Campaign Output — KRAS(G12D) pMHC Miniprotein Binders

**Target:** HLA-A\*11:01 / KRAS G12D 9-mer `VVGADGVGK` (PDB 9UV8)
**Off-target:** HLA-A\*11:01 / KRAS WT 9-mer `VVGAGGVGK` (PDB 8I5E)
**Objective:** α-helical miniprotein binders (≤120 aa) specific for p5 Asp (G12D mutation)

This directory (`/workspace/designs/`) holds all intermediate and final outputs from three
rounds of the design campaign. `/workspace/pmhc_design/designs/` holds copies of the
aggregated summary tables.

---

## Pipeline Overview

```
01_make_rfdiff_target.py       → prep/9UV8_target_chainB.pdb   (target PDB, single chain B)
02_run_rfdiffusion.sh          → round*/backbones/              (RFdiffusion backbone PDBs)
02b_run_partial_diffusion.py   → round*/backbones/              (partial-diffusion backbone PDBs)
05b_min_distance_filter.py     → round*/min_dist_*.csv          (backbone pre-filter)
03b_run_proteinmpnn_multiseq.sh → round*/mpnn_*/               (ProteinMPNN sequences, 8/backbone)
07_esmfold_filter.py           → round*/esmfold_*.csv / *.txt  (ESMFold monomer triage)
04_run_af2_initial_guess.sh    → round*/af2_complex_*/         (AF2 complex PDBs + scores)
08_mpnn_specificity.py         → round*/mpnn_spec_*.csv        (MPNN log-prob specificity)
09_make_peptide_only_ipsae_input.py → ipsae_peptide_only/      (binder+peptide PDBs for ipSAE)
10_run_pmhc_fold_3chain.py     → ipsae_af2_out/                (pMHC-fold 3-chain AF2 output)
prepare_af3_jobs.py            → af3_jobs_negdelta/            (AF3 Server batch JSON jobs)
15_compute_binder_metrics.py   → binder_metrics_13.tsv         (structural metrics)
16_compile_ipsae_gt0_stats.py  → final_designs/ipsae_gt0_stats/ (merged stats for 6 ipSAE>0 designs)
aggregate_summary.py           → final_designs/top_candidate_summary_stats.tsv
11_rosetta_reference_metrics.py → validated_reference_metrics.csv (calibration cutoffs)
```

---

## Scripts (`/workspace/pmhc_design/scripts/`)

### Preparation

**`01_make_rfdiff_target.py`**
Extracts chains A (MHC, residues 1–180) and C (peptide, residues 1–9) from `input/9UV8.pdb`,
merges them into a single chain B renumbered 1–189, and writes
`designs/prep/9UV8_target_chainB.pdb`. This is the fixed target fed to all downstream steps.
The peptide occupies residues 181–189 (p5 Asp = residue 185).

---

### Backbone Generation

**`02_run_rfdiffusion.sh`** `<output_prefix> <num_designs> [hotspot_res]`
Full de novo backbone generation via RFdiffusion. Hotspots default to `B184,B185,B186`
(peptide p4–p6, centered on the G12D Asp). Binder length sampled 70–100 aa; trajectories
suppressed (`write_trajectory=False`) to save disk.
*Outputs:* numbered `<prefix>_N.pdb` backbone files + a `trb` info file per design.

**`02b_run_partial_diffusion.py`** `--seeds_csv --backbone_dir --out_prefix [--variants_per_seed N] [--partial_T T] [--hotspot_res ...]`
Round-2/3 backbone diversification via partial diffusion. Reads a CSV of seed backbone tags,
noises/denoises each binder chain at `partial_T` steps (default 18/50) while holding the
target fixed. Binder length is fixed to match the seed. Per-seed `partial_T` override
supported via a `partial_T` column in the seeds CSV.
*Outputs:* `<out_prefix>_<seed>_<variant>.pdb` backbone files.

---

### Sequence Design

**`03_run_proteinmpnn_fastrelax.sh`** `<pdbdir> <outpdbdir> [runlist]`
ProteinMPNN sequence design with 1 FastRelax cycle (1 sequence per backbone).
Uses `dl_binder_design/mpnn_fr/dl_interface_design.py`. Binder = chain A (redesigned);
target = chain B (fixed).
*Outputs:* `<tag>_dldesign_0_cycle1.pdb` files in `<outpdbdir>`.

**`03b_run_proteinmpnn_multiseq.sh`** `<pdbdir> <outpdbdir> [num_seqs=8] [runlist]`
ProteinMPNN sequence design without FastRelax, generating multiple sequences per backbone
(default 8, matching the paper's supplement). Used for rounds 2–3 to expand sequence
diversity before the ESMFold triage step.
*Outputs:* `<tag>_dldesign_0.pdb` through `<tag>_dldesign_N.pdb` in `<outpdbdir>`.

---

### Filters

**`05_contact_filter.py`** `--pdbdir <dir> --out_csv <path>`
PyRosetta ContactMolecularSurface filter on designed/AF2-predicted structures (after MPNN,
when real side chains exist). Reports `peptide_cms` (binder vs all 9 peptide residues) and
`p5_cms` (binder vs Asp 185 only).
*Output:* CSV with columns `tag, peptide_cms, p5_cms`.

**`05b_min_distance_filter.py`** `--pdbdir <dir> --out_csv <path>`
Backbone-stage pre-filter (before MPNN). Computes minimum Cα–Cα distance (`min_ca_dist`)
and minimum backbone-atom distance (`min_bb_dist`) between the binder chain A and p5
residue 185. Does not require side chains.
*Output:* CSV with columns `tag, min_ca_dist, min_bb_dist`.

**`06_select_contact_pass.py`** `--csv <contact_csv> --out_runlist <txt> [--min_peptide_cms 10] [--min_p5_cms 1] [--mpnn_suffix ...]`
Converts `05_contact_filter.py` output to a runlist for `dl_interface_design.py`. Default
gates: `peptide_cms > 10` AND `p5_cms > 1` (generous, meant to kill clear non-contacts only).
*Output:* newline-delimited tag list (`.txt`).

**`07_esmfold_filter.py`** `--pdbdir <dir> --out_csv <csv> --out_runlist <txt>`
ESMFold (HuggingFace `facebook/esmfold_v1`) monomer foldability screen for all MPNN
sequences per backbone. Selects the single highest-pLDDT sequence per backbone and writes
its tag to a runlist for the AF2-complex step, keeping AF2 call volume at ~1/backbone.
*Outputs:*
- `*.csv` — pLDDT and Cα RMSD for every sequence
- `*.txt` — best-per-backbone tag runlist

**`08_mpnn_specificity.py`** `--pdbdir <dir> --out_csv <csv> --resnum 185 --chain B --on_target_aa D`
Local ProteinMPNN log-probability specificity screen. For each complex PDB (chain A = binder,
chain B = target), scores `log P(Asp|context)` and `log P(Ala|context)` at residue 185B.
`mpnn_spec_score = log_p(D) − log_p(A)`: higher = binder more specifically "reads" the
G12D mutation vs a neutral Ala. No AF3/off-target structure needed.
*Output:* CSV with columns `design, position, aa, log_p_on, log_p_ala, mpnn_spec_score`.

---

### ipSAE / pMHC-fold Confirmation

**`09_make_peptide_only_ipsae_input.py`** `<pdb_in> <json_in> <out_prefix>`
Prepares binder+peptide-only structures for the ipSAE (inter-protein SAE) metric.
Extracts chain A (binder) + peptide residues (181–189 of chain B), relabels the peptide
to chain C renumbered 1–9, and slices the corresponding PAE/pLDDT from the AF2 JSON.
*Outputs:* `<out_prefix>.pdb` and `<out_prefix>.json`.

**`10_run_pmhc_fold_3chain.py`**
Runs the pMHC-fold 3-chain AF2 scoring model from `paper/pMHCI_binder_design/`. Combines
the binder PDB (chain A) with the full target (MHC + peptide, chain B), rebuilds a
template PDB with continuously renumbered residues to avoid overlap, and calls
`predict_utils` to produce AF2 scores in the pMHC-fold framework.
*Outputs:* per-design AF2 score files in the `ipsae_af2_out/` directory.

---

### Rosetta / Structural Analysis

**`11_rosetta_reference_metrics.py`** `--inputs <glob> --out <csv> [--crystal <pdb>] [--relax]`
Calibrates Rosetta interface cutoffs from the 19 experimentally validated designs in
`paper_output/af3_nomsa/`. Runs `InterfaceAnalyzerMover` on FastRelaxed structures to
compute ddG, shape complementarity (sc), packstat, buried unsatisfied H-bonds (buns), and
a CMS proxy at p5. The resulting 10th-percentile/median values define the scoring-plan
cutoffs rather than literature guesses.
*Output:* `designs/validated_reference_metrics.csv` (and `.tsv` copy in `pmhc_design/designs/`).

**`12_check_p5_contact.py`** `<pdb_or_cif> [--two_chain]`
Checks binder-to-p5 minimum heavy-atom distance in an independently-folded AF3 structure
(chain layout: A=binder, B=MHC, C=peptide, D=b2m; p5 = chain C residue 5). Used to verify
that the p5 contact is robust in untemplated predictions, not a pose-bias artifact from the
templated AF2/pMHC-fold pipeline.

**`13_find_basic_to_asp.py`** `<pdb> [tag] [--two_chain]`
Single-structure scan for basic residues (Arg/Lys/His) on the binder chain A that are
within H-bond/salt-bridge range of the Asp p5 carboxylate oxygens (OD1, OD2). Prints any
hits to stdout. Used for manual inspection of candidate structures.

**`14_batch_find_basic_to_asp.py`**
Batch version of `13_find_basic_to_asp.py` hardcoded to
`designs/round3/spec_scan_winners/*.pdb`. Prints the closest basic–Asp pair per design and
a summary table of designs with at least one hit within 4 Å.

**`15_compute_binder_metrics.py`** `--input_dir <dir> --output_tsv <tsv>`
Computes per-design structural metrics from post-filter pMHC-fold PDB files (chain layout:
A=MHC groove, B=b2m, C=peptide, D=binder). Metrics:
- Binder Cα radius of gyration (`rg`)
- Number of helix segments, helix fraction (DSSP)
- Binder–peptide orientation angle (SVD principal axes)
- Lateral and vertical COM displacement from peptide COM (along MHC A principal axis)
*Output:* TSV; see `designs/binder_metrics_13.tsv`.

**`16_compile_ipsae_gt0_stats.py`**
Merges `af3_stats_13` (13-candidate AF3 run) with `af3_outward_asp_stats` (outward-Asp AF3
run for dldesign_6 and _7 variants) for the 6 designs with ipSAE > 0.1. Classifies p5 Asp
orientation (inward/outward) by chi1 dihedral (|chi1| < 120° → inward). Writes merged stats
to `/workspace/final_designs/ipsae_gt0_stats/`.

**`final_designs/aggregate_summary.py`**
Joins all per-candidate data sources — `final_designs/candidates_with_sequences.csv`,
`final_designs/ipsae_gt0_stats/af3_design_stats.tsv`,
`final_designs/ipsae_gt0_stats/fastrelax_af3_scores.tsv`, and
`final_designs/ipsae_gt0_stats/redesign_audit.tsv` — into a single 156-column TSV for the
6 ipSAE > 0 designs. Separately labels Asp inward and outward CMS scores and polar contacts
(`cms_p5_inward`, `cms_p5_outward`, `n_contacts_p5_polar_inward/outward`,
`contacts_polar_p5_inward/outward`). Requires pandas (run with `/venv/proteinmpnn/bin/python3`).
*Output:* `final_designs/top_candidate_summary_stats.tsv`.

---

### AF3 Job Preparation

**`prepare_af3_jobs.py`** `--candidates <csv> --out_dir <dir> [--batch_size 30]`
Generates AlphaFold Server (AF3) batch-upload JSON files for the G12D specificity screen.
For each candidate, emits two jobs: `<id>_on` (binder + 9UV8 on-target) and `<id>_off`
(binder + 8I5E off-target, WT peptide `VVGAGGVGK`). Uses the `alphafoldserver` v1 schema
(`proteinChain`). MHC heavy chain gets a real MSA (server-fetched); binder and b2m are
no-MSA; peptide uses `useStructureTemplate:False`.
*Output:* `af3_batch_NNN.json` files (≤30 jobs each) in `<out_dir>`.

---

## Output Directory Structure

### `prep/`
| File | Description |
|------|-------------|
| `9UV8_target_chainB.pdb` | Prepared design target: MHC+peptide merged to chain B, residues 1–189 (peptide = 181–189, p5 Asp = 185) |
| `make_rfdiff_target.py` | Copy of `01_make_rfdiff_target.py` used to generate the above |
| `pmhc_fold_weights/mixed_mhc_pae_run6_af_mhc_params_20640.pkl` | pMHC-fold AF2 weights |

---

### `round1/` — Full De Novo RFdiffusion

Generated ~1300 backbones (de novo, hotspots p4–p6) across two batches:
- `r1_*`: 300 backbones, 1 MPNN+FastRelax sequence each → `mpnn/`, AF2-complex in `af2_complex/`
- `r1b_*`: 1000 backbones, 8 MPNN sequences each (multiseq) → `mpnn_r1b/`

| Path | Description |
|------|-------------|
| `backbones/` | Raw RFdiffusion PDB backbones |
| `contact_filter_backbones.csv` / `contact_filter_r1b.csv` | CMS pre-filter scores (05_contact_filter.py) |
| `contact_pass_all.txt` / `contact_pass_r1_mpnn.txt` / `contact_pass_r1b.txt` | Runlists of backbones passing CMS gate |
| `mpnn/` | ProteinMPNN+FastRelax sequences for r1 batch (1/backbone) |
| `mpnn_r1b/` | ProteinMPNN multiseq (8/backbone) for r1b batch |
| `esmfold_r1b.csv` | ESMFold pLDDT + Cα RMSD for all r1b sequences |
| `esmfold_best_per_backbone_r1b.txt` | Best-per-backbone runlist from ESMFold triage |
| `af2_complex/` | AF2 initial-guess complex predictions for r1 batch (288 entries) |
| `af2_complex_r1b/` | AF2 initial-guess complex predictions for r1b ESMFold winners |
| `af2_complex_recheck32/` | AF2 re-run on 32 top specificity candidates |
| `af2_monomer_r1b/` | AF2 monomer predictions for r1b shortlist |
| `mpnn_spec_shortlist.csv` / `mpnn_spec_recheck_all32.csv` | MPNN specificity scores |
| `shortlist_pdbs/` | PDB copies of the r1b top candidates |
| `shortlist_r1b_tags.txt` | Tags of round-1b shortlisted designs |
| `pilot/`, `pilot_af2complex/`, `pilot_af2mono/`, `pilot_mpnn/` | Early pilot runs (pre-campaign) |
| `recheck_all_seqs/` | Full sequence re-evaluation batch |

**Round 1 outcome:** 5 designs passed `pae_interaction < 10`:
`r1_184`, `r1b_273_5`, `r1b_401_7`, `r1b_403_2`, `r1b_870_2`.
None passed the specificity gate (strong AF2 binding but near-zero `mpnn_spec_score`; all
contact was to the MHC framework, not the peptide). → Triggered round 2 with peptide-only
hotspot conditioning and partial diffusion.

---

### `round2/` — Partial Diffusion Diversification

Partial diffusion (T=18) on the round-1 top-5 seeds, 200 variants/seed = 1000 backbones.
Also a separate r1b_403 re-run with tighter hotspots.

| Path | Description |
|------|-------------|
| `backbones/` | Partial-diffusion backbone PDBs |
| `min_dist_r2.csv` | Backbone-stage p5 distance filter (05b_min_distance_filter.py) |
| `mindist_pass_r2.txt` / `contact_pass_r2.txt` | Runlists from backbone pre-filters |
| `mpnn_r2_full/` | ProteinMPNN multiseq (8/backbone) for full r2 set |
| `mpnn_r2_403/` | ProteinMPNN multiseq for r2 r1b_403 re-run |
| `mpnn_random8/` / `mpnn_top15/` / `mpnn_top21/` | MPNN sequences for sub-selections |
| `esmfold_random8.csv` / `esmfold_top15.csv` / `esmfold_top21.csv` | ESMFold screens |
| `esmfold_best_random8.txt` / `esmfold_best_top15.txt` / `esmfold_best_top21.txt` | Best-per-backbone runlists |
| `esm_r2_403.csv` / `esm_r2_403_winners.txt` | ESMFold for r1b_403 re-run |
| `af2_complex_random8/` (8) | AF2 complex for random 8 seeds |
| `af2_complex_top15/` (15) | AF2 complex for top-15 seeds |
| `af2_complex_top21/` (21) | AF2 complex for top-21 seeds |
| `af2_r2_403/` | AF2 complex for r1b_403 sub-batch |
| `mpnn_spec_r2.csv` / `mpnn_spec_r2_all.csv` | MPNN specificity scores |
| `mpnn_spec_screen/` / `mpnn_spec_screen_all/` | Complex PDBs used for MPNN spec scoring |
| `seed_candidates_top100.csv` | Top 100 seeds ranked for partial diffusion |
| `seed_candidates_top5_validated.csv` | Final top-5 seeds chosen for round-2 partial diffusion |
| `top_specificity_backbones.txt` / `top21_specificity_backbones.txt` | Tag lists for spec winners |
| `random_filtered_backbones.txt` | Random subset runlist |

---

### `round3/` — Expanded Partial Diffusion + MPNN Spec Screen

~1000 new backbones (partial diffusion from round-2 winners), full 8-seq MPNN, ESMFold
triage, and MPNN specificity screen before AF2-complex. The first campaign pass to run
the specificity screen upstream of AF2.

| Path | Description |
|------|-------------|
| `backbones/` | ~1000 partial-diffusion backbones (1001 files) |
| `esmfold_all1000.csv` | ESMFold scores for all ~8000 MPNN sequences (8/backbone) |
| `esmfold_all1000.log` | Run log |
| `esmfold_best_all1000.txt` | Best-per-backbone runlist (1000 tags) |
| `esmfold_top20.csv` | ESMFold scores for the top-20 specificity sub-batch |
| `esmfold_best_top20.txt` | Best-per-backbone runlist for top-20 sub-batch |
| `mpnn_all1000/` | ProteinMPNN multiseq (8/backbone) for all 1000 backbones |
| `mpnn_spec_screen/` / `mpnn_spec_screen2/` | Complex PDBs for MPNN spec scoring (all 994) |
| `mpnn_spec_r3.csv` (666 rows) | MPNN spec scores, first batch |
| `mpnn_spec_r3_part2.csv` | MPNN spec scores, second batch |
| `mpnn_spec_r3_combined.csv` | Merged spec scores |
| `mpnn_spec_winners994.csv` (994 rows) | Full spec scores for all round-3 sequences |
| `mpnn_spec_top8_all64.csv` | Spec scores for expanded 64-seq sub-batch of top 8 backbones |
| `af2_complex_all994/` (994 entries) | AF2 complex for all round-3 ESMFold winners |
| `af2_complex_top20/` | AF2 complex for top-20 specificity sub-batch |
| `af2_complex_top8_all64/` | AF2 complex for top-8 backbone expanded sequences |

---

### ipSAE Candidate Pipeline

After AF2-complex filtering, 13–19 candidates were taken through the pMHC-fold / ipSAE
confirmation pipeline.

| Path | Description |
|------|-------------|
| `ipsae_candidates/` | Binder+target complex PDBs for the ~19 AF3-candidate designs |
| `ipsae_peptide_only/` | Binder+peptide-only PDBs + AF2 JSON (for ipSAE; 09_make_peptide_only_ipsae_input.py output) |
| `ipsae_af2_out/` | pMHC-fold 3-chain AF2 predictions + contact contact analysis (10_run_pmhc_fold_3chain.py output) |
| `ipsae_af2_out/out.sc` | Score table for all ipSAE runs |

---

### AF3 Specificity Screen

**`af3_jobs_negdelta/`** — AF3 Server job files for the on-target / off-target specificity screen (13 candidates × 2 = 26 jobs).

| Path | Description |
|------|-------------|
| `af3_batch_001.json` | AF3 Server upload JSON (alphafoldserver v1 format) |
| `folds_2026_06_24_21_02.zip` | Downloaded AF3 Server results archive |
| `extracted/` | Unzipped AF3 results, one subdirectory per fold job (5 models + summary JSONs per job) |

**`af3_ipsae_13/`** — ipSAE metric files extracted from the 13 on-target AF3 folds.
One `.txt` file per design with the ipSAE score table (`ipsae_binder_peptide`, chain-pair
scores, contact counts).

---

### AF3 Full Statistics (`af3_stats_13/`, `af3_outward_asp_stats/`)

These directories store the post-FastRelax Rosetta analysis of the 13-candidate AF3 run
and the two "outward Asp" variant runs (`r3_r1b_870_87_dldesign_6` and `_7`), computed
by `15_compute_binder_metrics.py` and the Liu et al. pipeline's BindCraft analysis module.

| Path | Contents |
|------|---------|
| `af3_stats_13/af3_design_stats.tsv` | Full per-design stats table: iptm, ipSAE, CMS by peptide position, contact counts, H-bonds, salt bridges, Rosetta ddG/sc/packstat/buns, charge audit, sequence |
| `af3_stats_13/fastrelax_af3_scores.tsv` | Rosetta FastRelax scores + InterfaceAnalyzer metrics on relaxed AF3 structures |
| `af3_stats_13/redesign_audit.tsv` | Charge, Cys, Met redesign flags per design |
| `af3_stats_13/relaxed_pdbs/` | FastRelaxed PDB structures (4-chain: A=binder, B=MHC, C=peptide, D=b2m) |
| `af3_stats_13/contacts/` | Per-design contact TSVs (all, H-bond, hydrophobic, polar) |
| `af3_stats_13/targets.csv` | Target CSV used for this stats run |
| `af3_outward_asp_stats/` | Same layout as `af3_stats_13/`, covering the two outward-Asp AF3 re-folds |
| `af3_outward_asp/` | Raw AF3 output (CIF + JSON) for the outward-Asp jobs |

**`af3_outward_asp/`** — Raw AF3 outputs for `r3_r1b_870_87_dldesign_6_on` and
`r3_r1b_870_87_dldesign_7_on`, which were re-folded with an outward-Asp Arg-gated hotspot
variant. Contains CIF models and full-data JSONs (1 model each, single AF3 Server seed).

---

### Binder Metrics (`binder_metrics_input/`)

Relaxed PDB structures (FastRelax on AF3 models) for the 13 shortlisted candidates, used
as input to `15_compute_binder_metrics.py`. Chain layout: A=MHC groove, B=b2m, C=peptide,
D=binder (4-chain pMHC-fold convention).

---

### Final Designs (`/workspace/final_designs/`)

Top-of-funnel outputs for the 6 candidates with ipSAE > 0 (excluding
`r3_r1b_870_64_dldesign_0` at ipSAE = 0.012). This is the authoritative per-candidate
data for downstream analysis and reporting.

| Path | Description |
|------|-------------|
| `candidates_with_sequences.csv` | Master candidate table with binder sequences appended; covers all 21 evaluated designs with columns for every metric through `ipsae_af3_binder_pep`. Source of truth for sequence, round, and gate-pass status. |
| `ipsae_gt0_stats/af3_design_stats.tsv` | Full per-design structural stats for the 6 ipSAE>0 designs: iptm, ipSAE (observed and outward-Asp re-run), CMS per peptide position (p1–p9), inward/outward Asp chi1 classification, p5 contact counts (all/polar/hbond/hydrophobic/salt-bridge) split by Asp conformation, secondary structure, Rosetta interface metrics, BindCraft filter flags |
| `ipsae_gt0_stats/fastrelax_af3_scores.tsv` | Rosetta FastRelax scores + full InterfaceAnalyzer energy-term breakdown (fa_atr, fa_rep, fa_sol, fa_elec, hbond terms, fa_dun, etc.) for the 6 designs |
| `ipsae_gt0_stats/redesign_audit.tsv` | Charge, Cys, Met redesign flags; interface-fixed residue lists; pass/fail per design |
| `top_candidate_summary_stats.tsv` | **Aggregated summary (156 cols × 6 rows).** All columns from the three `ipsae_gt0_stats/` TSVs merged with off-target AF3 scores, MPNN spec scores, and pMHC-fold Δ from `candidates_with_sequences.csv`. Key column groups: identity, AF3 ipTM/ipSAE (observed + outward-Asp re-run), CMS per position, **inward/outward Asp CMS and polar contacts**, Rosetta interface + energy terms, BindCraft filters, charge/Cys audit. Generated by `aggregate_summary.py`. |
| `aggregate_summary.py` | Script that produces `top_candidate_summary_stats.tsv`; run with `/venv/proteinmpnn/bin/python3` |

**Note on r3_r1b_273_28_dldesign_0:** this design carries the `_control` suffix in its
relaxed PDB filename (`r3_r1b_273_28_dldesign_0_control_relaxed.pdb`) and is flagged
`Y_but_wrong_direction` because the p5 Asp adopts an outward conformation that was
predicted but not targeted during the binder generation run. It is included in the 6
ipSAE > 0 designs as a structural reference.

---

### Validated Reference (`validated_af3_staged/`, `validated_stats_out/`)

**`validated_af3_staged/`** — AF3 Server outputs for the 19 Liu et al. validated binders,
used to calibrate scoring thresholds. Each subdirectory (e.g. `ctnnb1-15/`) holds:
- `fold_<id>_model_0.cif` — top AF3 model
- `fold_<id>_full_data_0.json` — full confidence data (PAE, pLDDT)
- `fold_<id>_summary_confidences_0.json` — iptm, ptm, ranking score

**`validated_stats_out/`** — BindCraft pipeline stats for the 19 validated designs:
`af3_design_stats.tsv`, `fastrelax_af3_scores.tsv`, `redesign_audit.tsv`, `contacts/`,
`relaxed_pdbs/`. Same schema as `af3_stats_13/`. Used as the reference distribution by
`11_rosetta_reference_metrics.py`.

---

## Aggregated Summary Files (`/workspace/designs/` and `/workspace/pmhc_design/designs/`)

| File | Script that produced it | Contents |
|------|------------------------|----------|
| `candidates.csv` | manual / pipeline | Master candidate table. One row per design, columns: `candidate_id`, `round`, `seed_lineage`, `binder_len`, `net_charge`, `n_cys`, `n_met`, `af2_pae_interaction`, `pmhcfold_delta_on_minus_off`, `af3_on_iptm`, `af3_on_binder_pep_iptm`, `af3_off_iptm`, `af3_off_binder_pep_iptm`, `passes_af3_gates`, `real_cms_p5`, `their_ipsae_binder_pep_on/off`, Rosetta metrics, `mpnn_spec_score`, `ipsae_af3_binder_pep`, ESMFold + structural metrics |
| `final_designs/candidates_with_sequences.csv` | manual / pipeline | `candidates.csv` with `binder_sequence` column added; 21 evaluated designs; authoritative source for sequences + gate flags |
| `funnel_summary.csv` | manual | Round-by-round funnel: stage, n_in, n_out, pass_rate, gate, notes |
| `final_designs/top_candidate_summary_stats.tsv` | `final_designs/aggregate_summary.py` | 156-column × 6-row aggregated summary for the 6 ipSAE>0 designs; includes inward/outward Asp CMS and polar contact columns — primary file for reporting |
| `af3_design_stats.tsv` | `16_compile_ipsae_gt0_stats.py` | Merged full per-design stats for the 6 ipSAE>0 designs (from `af3_stats_13` + `af3_outward_asp_stats`), with inward/outward p5 Asp labeling; also at `final_designs/ipsae_gt0_stats/af3_design_stats.tsv` |
| `af3_negdelta_candidates.csv` | `prepare_af3_jobs.py` input | Candidate sequences submitted to AF3 Server for the on/off specificity screen |
| `af3_negdelta_results.csv` | `16_compile_ipsae_gt0_stats.py` | AF3 Server on/off results: iptm, binder_pep_iptm, delta values, gate pass flag |
| `binder_metrics_13.tsv` | `15_compute_binder_metrics.py` | Structural metrics for 13 candidates: Rg, n_helix_segs, helix_fraction, orientation_angle, lateral_dist, vertical_dist |
| `binder_metrics_13.log` | — | Run log for `15_compute_binder_metrics.py` |
| `real_cms_13_ontarget.tsv` | BindCraft pipeline | CMS per peptide position (p1–p9) for 13 on-target AF3 relaxed structures |
| `their_ipsae_13_onoff.tsv` | BindCraft pipeline | ipSAE (binder–peptide) on and off-target for the 13-candidate set |
| `tier1_af3_stats.tsv` | BindCraft pipeline | Full BindCraft stats table for the tier-1 candidates (same schema as `af3_design_stats.tsv`) |
| `round3_rosetta_metrics.csv` | `11_rosetta_reference_metrics.py` | Rosetta interface metrics for the round-3 shortlist |
| `round3_rosetta_metrics.log` | — | Run log |
| `validated_design_spec_scores.csv` | `08_mpnn_specificity.py` | MPNN specificity scores for Liu et al. validated designs (reference distribution) |
| `validated_reference_metrics.csv` | `11_rosetta_reference_metrics.py` | Rosetta ddG/sc/packstat/buns/dSASA reference values from the 19 validated hits |
| `ipsae_af2_out.log` | — | Run log for the pMHC-fold 3-chain AF2 run (10_run_pmhc_fold_3chain.py) |

Files with the same name in `/workspace/pmhc_design/designs/` are copies or earlier
versions of the above (produced before the working directory was restructured to
`/workspace/designs/`). The `af3_stats_13/` and `af3_outward_asp_stats/` subdirectories
in `pmhc_design/designs/` are earlier snapshots; `/workspace/designs/` versions are current.

---

## Design Naming Convention

`<round>_<seed_backbone>_<variant>_dldesign_<N>`

- `r1_184` — round 1, backbone 184 (de novo)
- `r1b_403` — round 1b (multiseq batch), backbone 403
- `r3_r1b_870_87` — round 3 partial diffusion, seeded from `r1b_870`, variant 87
- `_dldesign_N` — ProteinMPNN sequence index (0–7)
- `_on` / `_off` — on-target (9UV8 G12D) / off-target (8I5E WT) AF3 fold suffix

## Key Metrics Glossary

| Metric | Definition | Target threshold |
|--------|-----------|-----------------|
| `af2_pae_interaction` | AF2 initial-guess PAE across the binder–target interface | < 10 |
| `af3_iptm` | AF3 Server interface TM-score for the full complex | ≥ 0.90 |
| `af3_on_binder_pep_iptm` | AF3 chain-pair ipTM between binder and peptide chains | ≥ 0.77 |
| `ipsae_af3_binder_pep` | ipSAE metric (inter-protein SAE, binder–peptide pair) | > 0 (best ≥ 0.93) |
| `mpnn_spec_score` | `log P(D185) − log P(A185)` from ProteinMPNN; reads p5 Asp vs Ala | > 0 preferred |
| `pmhcfold_delta_on_minus_off` | pMHC-fold AF2 Δ(on−off) interaction score | > 0 preferred |
| `real_cms_p5` | Rosetta ContactMolecularSurface binder vs p5 Asp, on AF3 structure | > 0 required |
| `rosetta_ddG` | Rosetta InterfaceAnalyzer `dG_separated` on FastRelaxed structure | calibrate vs validated set |
| `rosetta_sc` | Rosetta shape complementarity at the interface | calibrate vs validated set |
