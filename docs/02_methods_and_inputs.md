# Inputs, Tools & Method (from Liu et al.)

> Covers the provided inputs, the cloned pipeline/tools, and the approach + parameters extracted from Liu et al., *Science* 2025. **Update this file as the paper/repo are read more closely** (e.g. refined cutoffs, additional repo conventions). Scientific rationale is frozen in [`01_scientific_background.md`](01_scientific_background.md); per-design results go in [`03_design_log.md`](03_design_log.md).

## Structural Input

**PDB: 9UV8**

The PDB is provided as-is. No cleaned or curated target structure is supplied. The 9UV8 construct is the peptide–MHC class I complex: the **HLA-A\*11:01 heavy-chain ectodomain (α1–α3)**, **β2-microglobulin**, and the **KRAS(G12D) 9-mer peptide (VVGADGVGK)** in the binding groove. It does **not** include the transmembrane/cytoplasmic regions of the heavy chain. Note that the deposited chains carry expression/biotinylation artifacts (e.g. a C-terminal AviTag–linker on β2-microglobulin) that are not part of the native complex.

The intended target surface is the solvent-exposed face of the bent 9-mer, with the mutant Asp12 (p5) at the apex — **not** the conserved MHC framework alone. A binder that contacts only the MHC helices will not be mutation- or peptide-specific.

**Correction (round 2, 2026-06-23):** earlier phrasing here ("together with the flanking HLA α1/α2 helices") was read as license to use HLA helix residues as RFdiffusion **hotspots** alongside the peptide — that's wrong and has been reverted. Hotspots must be **peptide-only**: conditioning RFdiffusion on HLA framework residues biases the docked interface toward the (easier, larger) conserved helix surface rather than the peptide, which is exactly the specificity failure mode this project screens for (round 1's candidates: strong AF2 binding, weak `mpnn_spec_score` vs. the `ctnnb1-15` positive control — see `docs/03_design_log.md`). A real binder's *final* footprint may still incidentally reach the flanking helices once docked on the peptide; that's a consequence of geometry, not something to steer toward via hotspots.

Participants are free to decide how to work with the structure — including which chains to use, whether to strip tags/β2-microglobulin, or how to preprocess or repack the groove. These decisions are considered part of the challenge.

## Sequences (9UV8 Construct)

**Target peptide — KRAS(G12D) 9-mer (chain C):**
```
VVGADGVGK
```

**HLA-A\*11:01 heavy chain, α1–α3 ectodomain (chain A):**
```
MGSHSMRYFYTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQTDRVDLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAAHAAEQQRAYLEGRCVEWLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWE
```

**β2-microglobulin (chain B; native sequence):**
```
MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM
```
*(In 9UV8 the deposited β2m chain continues `...GSGGSGAGLNDIFEAQKIEWHE` — a Gly/Ser linker plus AviTag for biotinylation, not part of native β2m.)*

Reference sequences: KRAS — NCBI RefSeq **NP_004976.1** (G12D substitution); HLA-A\*11:01 — IPD-IMGT/HLA; β2-microglobulin — NCBI RefSeq **NP_004039.1**.

## Design Constraints

- Each sequence must be **≤ 120 amino acids**
- Designs should make their predicted interface with the **peptide** (the flanking HLA helices are not a hotspot target — see correction above)

## Scoring

Designs are filtered and ranked on the metrics below. **You do not need to run all of these metrics** — use the subset that best discriminates promising designs at each stage of the funnel.

- **Folding confidence** (e.g. AF2/AF3 pLDDT, pTM)
- **RMSD** between the designed backbone and the predicted (refolded) binder after alignment
- **Net charge** of the binder
- **Secondary structure** composition (% residues in α-helix / β-sheet / loop)
- **CMS (contact molecular surface) scores for the 5th residue of the KRAS(G12D) peptide** (the mutant Asp12, peptide position 5 in the 9-mer)
- **Number of hydrophobic and polar interactions** between binder and peptide at the interface
- **ipSAE** between the peptide and the binder
- **Energetics after FastRelax** (Rosetta interface energy / ddG)
- **Shape complementarity (sc) and packstat** between binder and target
- **Unsaturated hydrogen bonds at the interface**
- **Buried unsaturated hydrogen bonds within the complex**
- **Number of lysines and methionines at the binding interface** (developability / oxidation liabilities)
- **Fraction of exposed hydrophobic residues**
- **Presence of poly-alanine stretches** (a flag for low-information / degenerate designs)

**Specificity is the primary objective, not raw interface score.** ~85% of the pMHC surface is conserved MHC, and the only therapeutically meaningful contact is with the **mutant peptide (Asp12, peptide position 5)**. A high absolute ipSAE or interface ddG against the whole complex will reward binders that simply clamp the **conserved HLA-A\*11:01 framework** — high-scoring but specificity-blind and clinically useless.

Designs are therefore **counter-screened against wild-type KRAS / HLA-A\*11:01** (peptide VVGAGGVGK on the same allele) and judged on the **differential**: a design must show predicted-interface preference and direct contact to the **G12D Asp** over WT, not merely a high absolute score. The **CMS at peptide position 5** is the key per-residue handle for this. Designs that score well only via MHC-framework contacts are deprioritized regardless of absolute ipSAE.

Experimental validation of binding and specificity is downstream and out of scope for the current computational run.

Reference: [ipSAE on GitHub](https://github.com/DunbrackLab/IPSAE)

## Resources & Inputs Provided

Everything needed to start is in `/workspace`. The agent is a **computational protein design scientist** whose job is to design miniprotein binders to the KRAS(G12D) pMHC — **learning from the experimentally validated designs** in the reference paper rather than starting from scratch.

**Target inputs (`input/`):**

| File | Role |
|------|------|
| `input/9UV8.pdb` | **Target** — KRAS(G12D) 9-mer (VVGADGVGK) / HLA-A\*11:01. Chains A (HLA heavy), B (β2m), C (peptide). The complex to design against. |
| `input/8I5E.pdb` | **Wild-type off-target** — wild-type KRAS peptide on the same MHC, for specificity counter-screening. Chains H / L / P. |
| `input/pmhc_fold_scaffold/*.pdb` | **HLA-A\*11:01 pMHC template PDBs** for `paper/pMHCI_binder_design/pMHC_fold/pmhc_fold.py` (`6jtp_prep`, `6o4y_prep`, `8i5c_prep`, `8i5d_prep`, `8i5e_prep`). These are the **already-prepared output** of the template search — so it does **not** need to be re-run. (The source table `pMHC_fold/alignments_all_alleles_vs_pdb_June29_2024.csv` has been **removed to save space**.) |

**Reference approach & validated-design materials (context — `paper/`, `paper_output/`):**

| File | Role |
|------|------|
| `paper/pMHCI_binder_design/` | **Cloned code repo** for Liu et al. (the paper's official pipeline). Modules: `pMHC_fold/` (specificity folding), `contact_filter/` (CMS), `mpnn_spec_filter/` (ProteinMPNN specificity), `scaffolds/` (recyclable scaffold library), `NGS_analysis/`, `software/`. |
| `paper/science.adv0185.pdf` | Main paper — Liu et al., *Science* 2025, "Design of high-specificity binders for peptide-MHC-I complexes." |
| `paper/science.adv0185_sm.pdf` | Supplement — method details, parameters, and filter cutoffs. |
| `paper_output/science.adv0185_data_s1.xlsx` | Sequences of the **19 experimentally validated** binders. |
| `paper_output/af3_nomsa/` | AF3 (no-MSA binder) predicted structures of the 18 validated binders — **provided by the user**, not from Liu et al.'s released materials. Paired with `designs/af3_design_stats.tsv` (user-provided, comprehensive per-design metrics including real CMS per peptide position, ipSAE, DSSP, hotspot positions, Rosetta IAM output). |
| `paper_output/af3_nomsa/binder_pMHC_full_noMSA_MHC.json` | **Real AF3 Server submission example** — confirms the actual input schema is `dialect: "alphafoldserver"`, `version: 1`, chains as `"proteinChain"` (not `"protein"`/`"alphafold3"` v2, which is the *output* `data.json` schema — `prepare_af3_jobs.py` used the wrong one until 2026-06-24). This example also strips MSA from the MHC chain; **we deliberately don't** — keeping a real MSA for MHC avoids false-negative binding predictions (user direction). |
| `paper_output/design_hotspots.csv` | Per-target binding-hotspot analysis (Claude-generated from the paper) — which peptide positions each validated binder reads. |
| `paper_output/mmdb_9O5S.pdb` | An **experimentally determined** structure of a validated binder–pMHC complex (ground-truth check vs. the AF3 model). |

## Reference Approach — Liu et al., *Science* 2025

The pipeline is **structure-based design, not sequence-based**: RFdiffusion generates novel binder *backbones* docked against the target structure, and ProteinMPNN designs sequences onto those backbones — we do not start from any existing binder sequence. The workflow follows Liu et al.: **RFdiffusion** (hotspot-conditioned backbone generation) → **ProteinMPNN** (sequence design) → **ESMFold** (foldability triage) → **AF2 initial-guess** (binding check) → **pMHC-fold** (on/off-target specificity) → **AF3 (no-MSA)** (independent confirmation) → **Rosetta InterfaceAnalyzer** (interface energy, shape complementarity, buried unsatisfied H-bonds; run after AF3, *not* from Liu et al.'s pipeline — added at user direction for calibrated interface metrics beyond Liu et al.'s CMS). The same machinery is applied here to 9UV8, with the KRAS-G12D-specific twist that the discriminating residue is the **mutant Asp at peptide position 5**.

**AF3 (no-MSA) structures of the 18 validated designs** (`paper_output/af3_nomsa/`): these were provided by the user and serve as the reference set. They are not from Liu et al.'s released materials.

## Lessons from the Experimentally Validated Designs

These are the empirical patterns from the 18 validated hits (`data_s1.xlsx`, `af3_nomsa/`, `af3_design_stats.tsv`) that should shape our KRAS designs. Note: each validated design targets a **different peptide on a different allele** — "p5" means **position 5 of that design's own target peptide**, which is a distinct residue and mutation in every case. The common principle is that the *neoepitope mutation*, whatever it is and wherever it falls in the peptide, must be contacted directly.

- **Format:** every validated binder is a **small α-helical miniprotein, 59–106 aa** (within the 120-aa cap). Several reuse a privileged "mage-513" scaffold; others are de novo. Helical bundles, not β or loop-heavy folds.
- **Peptide-centric docking:** each binder reads only **1–4 upward-facing peptide residues**, and **the neoepitope mutation itself is the primary specificity determinant** — e.g. `ctnnb1-15` (A\*03:01, TTAPFLSGK) is driven by **F5**, the S5F mutation (cation-π/π-π/hydrophobic). This is the **direct analog of our target**: same A3-supertype family, 9-mer with C-terminal Lys anchor, single mutated residue at **p5** as the specificity handle. For KRAS that residue is **Asp (p5)** — a charged side chain, so favor a complementary basic/H-bonding pocket (cf. `phox2b-5` R6-bidentate, `wt1-5` R1-bidentate via binder Glu).
- **AF3 (no-MSA) confirmation thresholds** (observed across all 19 hits): overall **iptm 0.92–0.97**, **ptm 0.94–0.97**, with **binder–peptide chain-pair iptm ≈ 0.77–0.91**. Treat **iptm ≥ 0.90** + strong binder–peptide chain-pair iptm as the bar for the final shortlist.
- **Specificity is enforced explicitly**, not assumed — see the specificity filters in the pipeline below.
- **Expect multiple rounds.** High-quality binders rarely come from a single pass — plan on several design rounds (generate → filter → inspect survivors → adjust hotspots/scaffolds/cutoffs → regenerate). Per-round yield through the AF2 → CMS → specificity funnel is low, and the paper re-tuned cutoffs **per target and per round**; the final shortlist is the survivor of repeated narrowing, not one generation.

## Translating the validated designs into evaluation criteria (KRAS G12D)

Quantitative patterns from the 18–19 validated binders (`data_s1.xlsx`, `af3_nomsa/`, `design_hotspots.csv`, crystal `mmdb_9O5S.pdb`, and Rosetta metrics in `designs/validated_reference_metrics.csv`), turned into pass/flag criteria for KRAS(G12D) candidates. Numbers are guidance from a mixed-allele set (only `ctnnb1-15`/`phox2b` are close A3-supertype/9-mer analogs), not hard law.

**`designs/af3_design_stats.tsv`** (user-provided, 2026-06-25) is now the most comprehensive validated-design reference: real PyRosetta `ContactMolecularSurface` per peptide position (not the lightweight atom-count proxy below), DSSP secondary structure, multi-position `hotspot_positions` per design (most hits read 2-4 peptide positions, not only p5), salt-bridge detection, and full InterfaceAnalyzer output. Prefer it over `validated_reference_metrics.csv` going forward.

**Sequence criteria**
- **Length:** 59–106 aa (median ~95); compact A3/C-allele analogs are 59–83 aa. Target **~60–100 aa** (cap 120).
- **Net charge: negative.** All hits are −3 to −11 (median −7). **Flag any net-positive design.** (KRAS p5 is Asp⁻, so the binder still needs a *local* basic/H-bond donor at the p5 contact while net-acidic overall.)
- **No cysteines.** 0 across all hits — treat a Cys as a near-automatic reject.
- **Methionine ≤ 2; no poly-Ala run > 3.** Flag Met-rich (oxidation) or poly-Ala (degenerate) designs.
- **Helical composition:** ~35–49% hydrophobic, helix-favoring; expect a helical bundle, not β/loop-heavy.

**Topology & interface criteria**
- **Single-chain α-helical miniprotein** spanning the groove, contacting both MHC α-helices *and* the peptide (9O5S confirms a 102-aa helical binder). Don't penalize MHC contact per se — specificity comes only from the peptide contacts.
- **Must contact the neoepitope mutation.** Each of the 18 validated designs contacts the mutation-bearing residue of *its own* target peptide (the hotspot positions vary — e.g. p5 in ctnnb1/mage/sars/hiv-9, p4-p7 in gp100-3, p1/p8 in hiv-10/wt1). For KRAS(G12D), **that residue is p5 (Asp12)**. **A candidate that does not make a direct contact to Asp(p5) is rejected regardless of other scores.** (Caveat: direct contact here means within real CMS/heavy-atom range — the AF3 structural results from this campaign show that backbone H-bonds to the p5 carbonyl are real but do not by themselves prove G12D specificity; ideally a binder Arg/Lys/His side chain H-bonds or salt-bridges the Asp carboxylate (OD1/OD2), though no current candidate achieves this in the AF3 structure — see `docs/03_design_log.md` 2026-06-25.)
- **Interface size:** validated binders bury large interfaces — ~7–14 binder residues at the peptide interface and **ΔSASA 2100–3300 Å²** (median ~2900). Contact **2–5 peptide positions including p5 and ≥1 flanking** (p4 and/or p6–p8). **Flag interfaces < ~1800 Å²** as too small.
- **Provide a complementary partner for Asp p5:** position an Arg/Lys/His or backbone amide to H-bond/salt-bridge the Asp carboxylate — mirroring `wt1-5` (binder Glu↔peptide Arg1) and `phox2b-5` (bidentate to Arg6). Highest-value feature for mutation specificity.

**Rosetta interface cutoffs — calibrated, not guessed** (18 validated AF3 designs, FastRelax + InterfaceAnalyzer; `designs/validated_reference_metrics.csv`, via `scripts/rosetta_reference_metrics.py`, PyRosetta in the **`dl_binder_design`** env):
- `rosetta_ddG` (dG_separated): validated −57 to −88 REU (median −71). **Gate: < −50 REU**; target ≤ −70.
- `shape_complementarity`: 0.51–0.71 (median 0.62). **Gate: ≥ 0.50**; target ≥ 0.62.
- `packstat`: 0.56–0.71 (median 0.63). Target ≥ 0.60.
- `cms_Asp_p5` **(proxy)**: `scripts/rosetta_reference_metrics.py`'s `cms_p5_proxy` is a lightweight stand-in (binder heavy-atom-to-p5 contact count, not real CMS) — every hit makes **≥ 2 such contacts** (median ~2-4). **Gate: ≥ 2 (target ≥ 4).** Don't compare this number to `cms_p5` in `af3_design_stats.tsv` below — different metric, different scale.
- `cms_p5` **(real PyRosetta ContactMolecularSurface, Å²)**, from `af3_design_stats.tsv`: validated range **10.7–89.7, median ~47.9**. Use this for any future comparison against real CMS computed on our candidates (requires running the actual `ContactMolecularSurface` filter, not the proxy).
- `buried_unsat_hbonds`: median 12, up to 22 — **diagnostic only**, do not hard-gate low (large pMHC interfaces carry some unsat).
- NOTE: these are protocol-specific — score KRAS candidates with the **same** FastRelax + InterfaceAnalyzer pipeline for the cutoffs to transfer.
- **ipSAE caveat**: `af3_design_stats.tsv`'s `ipsae_binder_peptide` column (validated range 0.86-0.97) uses a different, simplified d0 formula (global contact-count-based, not the per-residue Schaeffer & Dunbrack `d0res`) than the official `DunbrackLab/IPSAE` `ipsae.py` used elsewhere in this project (validated range 0.42-0.87 at pae/dist=10/10). Confirmed by directly importing and running both on identical AF3 PAE data — not a cutoff difference, a genuinely different statistic. Don't mix the two when comparing candidates.

**Confidence criteria (AF3 no-MSA confirmation)**
- **iptm ≥ 0.90** (validated 0.92–0.97), **ptm ≥ ~0.94**, **binder–peptide chain-pair iptm ≥ 0.77** (validated 0.77–0.91).

**Specificity (the objective)**
- Counter-screen G12D vs WT (`8I5E`): require a positive G12D-minus-WT differential (CMS@p5 and Δ pAE_interaction). A binder that scores equally on WT is reading the MHC/backbone, not the mutation.

## Design Pipeline (parameters from the paper supplement)

> Scripts/configs will live under this repo and be committed as built. **Compute:** local **NVIDIA RTX 5090** (single GPU); a pre-provisioned **conda environment** on the GPU host carries the model stack (standard Python packages may be `pip`/`conda` installed as needed).

> **When applying any model or filter (RFdiffusion, ProteinMPNN, pMHC-fold, contact/specificity filters, scaffold docking, etc.), follow the cloned repo's own instructions** in `paper/pMHCI_binder_design/README.md` and the per-module folders — use its documented commands, scripts (e.g. `align_chainB.py` for docking scaffolds onto the target), input conventions (renumber the pMHC from residue 1, relabel MHC+peptide to chain B), and recommended settings rather than improvising invocations. Note the repo's **2026 update: AF3 is now recommended for folding pMHCs** — prefer it for the structure-prediction/confirmation steps.

1. **RFdiffusion** — generate binder backbones docked against the full truncated 9UV8 target (HLA α1/α2 + peptide, as the structural context RFdiffusion docks against), but **hotspots restricted to peptide residues only** (never HLA helix residues — see correction above), conditioned toward high contact at the exposed Asp(p5) apex. Use the recyclable scaffold library with **partial diffusion `partial_t` 12–25** (of 50 denoising steps) after docking scaffolds to the new target. **1,000–20,000 backbones** per RFdiffusion cycle; iterate partial diffusion until enough designs pass.
2. **ProteinMPNN** — sequence design on the output backbones.
3. **AF2 initial-guess** — primary binding filter: **pAE_interaction (general cutoff < 5)**, binder pLDDT, and binder RMSD to the design.
4. **Monomer foldability** — binder-alone pLDDT (foldability of the designed chain in isolation). Uses **ESMFold**, not AF2's `-force_monomer` mode: ESMFold (no MSA, ~1-2s/seq) is a faithful substitute for this specific check and is already computed as a byproduct of `07_esmfold_filter.py`'s per-backbone sequence triage (vs. AF2 monomer's ~180s/design for the same number) — see `docs/03_design_log.md` round-2 entries. AF2's initial-guess templating trick (step 3) has no ESMFold equivalent and is *not* substituted.
5. **Specificity screening (the crux):**
   - **pMHC-fold** each design against the MHC loaded with the **on-target (G12D)** peptide and **2–3 off-target peptides in the same allele** (wild-type KRAS via `8I5E`, plus other G12 mutants). Compute **peptide-only ↔ binder pAE_interaction**, then a **Δ pAE_interaction = on-target − off-target** (more negative = more specific) to filter.
   - **ProteinMPNN log-probability / alanine-scan**: score the likelihood of each peptide residue (esp. Asp p5) when the binder is present; filter on the on-target log-prob or the on−off log-prob difference. **Rule (from paper author Bingxu Liu, 2026-06-23): never rank/select *backbones* with `mpnn_spec_score` — only use it to rank *sequences within the same backbone*.** Cross-backbone comparisons are confounded by backbone-specific baseline/scale effects unrelated to true specificity; verified empirically in round 2 (see `docs/03_design_log.md`).
6. **AF3 (no-MSA) confirmation** — fold the shortlist (binder + pMHC, no MSA for the binder; real MSA for MHC to avoid false-negative binding) and keep those meeting the validated-design thresholds (iptm ≥ 0.90; strong binder–peptide chain-pair iptm ≥ 0.77). Run on/off-target (G12D and WT peptide) for both absolute confidence and specificity delta. Run on the **AlphaFold3 web server (≤ 30 jobs/day** — stage accordingly). This orthogonal check matters: AF2/pMHC-fold are pose-templated (biased toward the designed coordinates) and AF3 without templating is the only fully independent structural prediction in the pipeline.
7. **Rosetta InterfaceAnalyzer** (*user-directed, after AF3; not in Liu et al. pipeline*) — FastRelax the AF3 on-target structure (`scripts/11_rosetta_reference_metrics.py --relax`), then run InterfaceAnalyzer for ddG, shape complementarity, packstat, buried-unsat H-bonds, and **real CMS per peptide position** (`ContactMolecularSurface`, Å²-scale; calibrated against the 18 validated designs via `designs/af3_design_stats.tsv`, range 10.7–89.7 at p5, median 47.9). Also check polar contacts and H-bonds to p5 side chain (OD1/OD2) to verify the binder is reading the Asp mutation rather than only the invariant backbone.

**This is iterative.** Steps 1–6 are expected to run over **multiple design rounds**: inspect what survives each round, re-tune hotspots/scaffolds/cutoffs (cutoffs are target- and round-specific), and regenerate until enough designs clear the funnel. A single pass will not produce high-quality binders. **Across rounds, it is fine to use the evaluation metrics in the [Scoring](#scoring) section to guide and improve the designs** — e.g. steer hotspots/sequence design toward better CMS at Asp(p5), stronger shape complementarity/ddG, cleaner secondary structure, fewer buried unsatisfied H-bonds, and a larger on- vs off-target specificity margin — using them as optimization targets between rounds, not just as a final ranking.

**Output:** a ranked shortlist of **< 50 prioritized designs** (≤ 120 aa each) with the supporting metric table and AF3-confirmed predictions.

### Scoring plan for the run (decided)

The Scoring section lists everything that *can* be computed; the metrics actually used during the run, and how, are fixed here. Each stage has a small number of **hard gates** (fail → drop, don't score downstream), with the rest kept as **diagnostics** that rank survivors or flag liabilities. Cutoffs are starting points and may be re-tuned per round (see the iterative note).

| Stage | Metric | Use | Starting cutoff |
|---|---|---|---|
| Backbone (RFdiffusion) | peptide-contact filter (repo `contact_filter`) | **gate** | adequate contact over peptide incl. p5 |
| AF2 initial-guess | `af2_pae_interaction` | **gate** | **< 5** |
| AF2 initial-guess | `binder_rmsd` (design→pred) | **gate** | **< 2 Å** |
| ESMFold (substitutes AF2 monomer) | `af2_monomer_pLDDT` (column name kept; value is ESMFold's) | **gate** | **≥ 80** |
| Rosetta/contact | `cms_Asp_p5` (CMS on mutant Asp) | **gate** | ≥ 2 binder contacts to p5 (target ≥ 4) |
| Rosetta/contact | `rosetta_ddG` (post-FastRelax) | **gate** | **< −50 REU (validated median −71)** |
| Rosetta/contact | `shape_complementarity` | diagnostic/rank | ≥ 0.50 (target 0.62) |
| Rosetta/contact | `buried_unsat_hbonds` | diagnostic/rank | fewer better |
| Specificity | `delta_pae_specificity` (G12D − WT/other) | **gate** | **< 0** (more negative better) |
| Specificity | `mpnn_spec_score` (Ala-scan at p5) | **within-backbone rank only — never a cross-backbone gate** | pick the best sequence per backbone, not which backbones advance |
| AF3 (no-MSA) | `af3_iptm` | **gate** | **≥ 0.90** |
| AF3 (no-MSA) | `af3_binder_pep_iptm` | **gate** | **≥ 0.77** |
| Final ranking | `ipSAE` (peptide↔binder) | **rank only** | compute on AF3-confirmed set; primary leaderboard sort |
| Developability | net charge, frac exposed hydrophobic, interface Lys/Met, polyAla, SS % | diagnostic/flag | applied to shortlist, not early gates |

**Rationale:** the hard gates are the cheap, high-signal filters the paper relied on (AF2 `pae_interaction`, CMS on the mutant residue, the specificity Δ) ordered cheapest-first so most designs die before expensive steps. `ipSAE` is **not** an early gate (it correlates with the AF3/pae metrics and is costly); it is computed only on the AF3-confirmed set and used as the final leaderboard sort. Developability metrics never filter early — they only break ties and flag liabilities on the shortlist.

**If too few designs pass, loosen the gates.** These cutoffs are not hard requirements — if a round yields too few survivors to be useful, it is fine to **relax the gates somewhat** (e.g. accept `af2_pae_interaction` slightly above 5, or `af3_iptm` modestly below 0.90) to keep a working set of the **best-available** candidates, then **carry the top-ranking ones into the next round and improve them** (partial diffusion / resampling around them, sequence redesign, hotspot adjustment) until they clear the original thresholds. Prefer advancing and refining the current best over discarding everything and stalling; tighten the gates back up as design quality improves across rounds.

## Alternative methods (optional)

The Liu et al. funnel above is the primary approach, but two recent binder-design methods are worth keeping in reserve — especially for the **specificity** problem:

- **PXDesign** (ByteDance; [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.08.15.670450v2), [server/repo](https://github.com/bytedance/PXDesign)). Same generate→sequence→filter shape as ours, but native to the **Protenix** (AF3-style) model: a diffusion backbone generator (PXDesign-d) → ProteinMPNN → filtering by **ipTM with Protenix + AF2-initial-guess** and **Foldseek** diversity clustering. Reports ~17–82% nanomolar hit rates, beating RFdiffusion/AlphaProteo on benchmarks. Turnkey and high-throughput; a drop-in alternative to the RFdiffusion→AF2 stages.
- **Mosaic** (escalante-bio; [repo](https://github.com/escalante-bio/mosaic)). A **gradient-based composite-objective** framework: the binder sequence is a continuous distribution optimized *through* multiple predictors at once (AF2 / Boltz / Protenix + ProteinMPNN + pLMs + handcrafted terms) summed into one differentiable loss, then discretized. Won the Adaptyv Nipah de novo track and produced the top TREM2 binder (1.11 nM). **Most relevant twist for us:** the **G12D-vs-WT specificity margin** (and CMS@p5) could be folded *directly into the optimization objective* rather than enforced only as a downstream filter — a route to bake mutation-specificity into design from the start.

These are not set up in this repo; use the Liu et al. pipeline by default and consider these only if it stalls on specificity or yield.

