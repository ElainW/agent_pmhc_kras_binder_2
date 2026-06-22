# pmhc_design scripts

Scripts and analysis code for the KRAS(G12D) pMHC miniprotein binder design campaign
(target: PDB 9UV8, HLA-A*11:01). See `/workspace/CLAUDE.md` and `/workspace/docs/` for
full project context, conventions, and the per-round design log.

Heavy outputs (RFdiffusion backbones, ProteinMPNN sequences, AF2/AF3 structures) live
under `/workspace/designs/<round>/` and are not committed here; this repo holds the
reusable pipeline code.

## Scripts

Numbered scripts mirror the order they're run in; each is also the audit trail for the
corresponding commands actually executed during the campaign.

- `scripts/00_setup_dl_binder_design_env.sh` — one-time environment fixes for
  `dl_binder_design` (see "Environment notes" below); idempotent, safe to re-run.
- `scripts/01_make_rfdiff_target.py` — builds the RFdiffusion docking target from `input/9UV8.pdb`:
  HLA-A*11:01 heavy chain residues 1-180 + KRAS(G12D) peptide residues 1-9, renumbered
  continuously (1-189) and relabeled to chain B (matches `contact_filter` convention: peptide
  = residues 181-189, p5/Asp(G12D) = residue 185). Drops beta-2-microglobulin (not part of the
  peptide-binding groove).
- `scripts/02_run_rfdiffusion.sh <output_prefix> <num_designs>` — hotspot-conditioned binder
  backbone generation (hotspots B184-186 / peptide p4-p6, centered on Asp(p5)=B185; see
  `docs/03_design_log.md` round-1 entries for the rationale).
- `scripts/03_run_proteinmpnn_fastrelax.sh <pdbdir> <outpdbdir> [runlist]` — ProteinMPNN sequence
  design + 1 FastRelax cycle, 1 sequence/backbone (binder = chain A, redesigned; target = chain B,
  fixed). Used for the r1 (300-backbone) pilot batch.
- `scripts/03b_run_proteinmpnn_multiseq.sh <pdbdir> <outpdbdir> [num_seqs_per_struct] [runlist]` —
  ProteinMPNN **without** FastRelax, N sequences/backbone (default 8, paper's supplement range is
  4-32 — see "MPNN multiplicity gap" in `docs/03_design_log.md`). `dl_interface_design.py`
  disallows `seqs_per_struct>1` with FastRelax on, and FastRelax wasn't found to matter much
  in silico per RFdiffusion's own README, so this trades it for the paper's intended per-backbone
  sampling depth. Follow with `07_esmfold_filter.py` to triage to 1 sequence/backbone before AF2.
- `scripts/04_run_af2_initial_guess.sh <pdbdir> <outdir> [complex|monomer] [runlist]` — AF2
  initial-guess filter: complex mode scores `pae_interaction`/`binder_rmsd`, monomer mode scores
  binder-alone `pLDDT`.
- `scripts/05_contact_filter.py --pdbdir <dir> --out_csv <path>` — cheap PyRosetta
  `ContactMolecularSurface` pre-filter on raw RFdiffusion backbones (binder vs. peptide
  181B-189B, and vs. Asp/p5=185B alone) — ~1-2s/design vs AF2's ~200s. Run before MPNN to
  avoid spending MPNN+AF2 compute on backbones that don't reach the peptide.
- `scripts/06_select_contact_pass.py --csv <05 output> --out_runlist <path>` — turns the
  contact-filter scores into an AF2/MPNN `-runlist` of passing backbone tags.
- `scripts/07_esmfold_filter.py --pdbdir <mpnn_multiseq_out> --out_csv <path> --out_runlist <path>` —
  folds every ProteinMPNN sequence (binder chain alone) with ESMFold (`facebook/esmfold_v1`,
  HuggingFace `transformers`, no MSA) and keeps only the best-pLDDT sequence per backbone,
  written as an AF2-complex `-runlist`. This is what makes the N-sequences/backbone multiplicity
  affordable — AF2-complex call volume stays at ~1/backbone regardless of N.
- `scripts/prepare_af3_jobs.py` — given a CSV of candidate binder sequences that passed the
  AF2 + ProteinMPNN-specificity filters, generates AlphaFold Server (AF3) batch-upload JSON
  files for manual submission: one on-target (9UV8, peptide VVGADGVGK) and one off-target
  (8I5E, peptide VVGAGGVGK) job per candidate, binder chain folded with no MSA, batched to the
  AF3 Server's 30-job/day quota. Drop the resulting structures back under
  `designs/<round>/af3_jobs/<batch>/results/` for downstream scoring.

## Environment notes (this Vast.ai instance)

- **CPU affinity:** the container's pids cgroup caps at 1024 (`/sys/fs/cgroup/pids.max`)
  but `nproc` reports 256 visible cores. JAX/TensorFlow size their CPU threadpools to
  the visible core count and blow past the pids cap, crashing with
  `pthread_create() failed (11)`. Always wrap jax/tensorflow processes (AF2
  initial-guess, pMHC-fold) with `taskset -c 0-15` (or similar) to cap visible cores.
- **AF2 weights:** `dl_binder_design/af2_initial_guess/predict.py` expects
  `model_weights/params/params_model_1_ptm.npz` (not present in the base image — download
  from `https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar`,
  extract just `params_model_1_ptm.npz`, ~5.5GB tar/~370MB needed file).
- **ProteinMPNN dependency:** `dl_binder_design/mpnn_fr/dl_interface_design.py` expects a
  `ProteinMPNN` checkout inside `mpnn_fr/` — symlinked to the project's own
  `/workspace/ProteinMPNN` rather than re-cloning.
- **`dl_binder_design` venv was missing `torch`** (needed by `dl_interface_design.py`
  alongside `pyrosetta`) — installed `torch` (cu128 build) into that venv; coexists fine
  with its existing `jax`/`pyrosetta`.
- **Blackwell (RTX 5090, sm_120) + vendored AlphaFold/jax compatibility:** the
  `af2_binder_design` venv's original `jax==0.4.28` ships a `ptxas` that can't target
  sm_120 (`cannot be compiled to future architecture`). Upgraded to latest `jax[cuda12]`
  (0.10.2) to get Blackwell kernel support, which in turn requires patching the vendored
  2022-era AlphaFold code for jax API removals (`jax.tree_map` family, `jax.util.wraps`,
  `jax.lib.xla_bridge`, `jnp.clip(a_min=/a_max=)`). See
  `patches/dl_binder_design_blackwell_jax_compat.patch` (applied directly to the
  `/workspace/dl_binder_design` checkout) — re-apply with
  `git -C /workspace/dl_binder_design apply patches/dl_binder_design_blackwell_jax_compat.patch`
  if that checkout is ever reset.
- **ESMFold** (`07_esmfold_filter.py`) installed via `pip install transformers accelerate` into
  the `proteinmpnn` venv (already had Blackwell-compatible `torch`) — uses
  `transformers.models.esm.modeling_esmfold.EsmForProteinFolding` (`facebook/esmfold_v1`), not
  Meta's original `fair-esm`+`openfold` (avoids a painful custom-CUDA-kernel build). Needs
  `taskset -c 0-15` for the same pids-cgroup reason as AF2. HF token cached at
  `/workspace/.hf_home/token` (per-instance `HF_HOME`, not in any repo) for faster downloads.
  CPU inference works but is impractically slow (~minutes/sequence) — GPU is the point of using
  ESMFold here at all.
