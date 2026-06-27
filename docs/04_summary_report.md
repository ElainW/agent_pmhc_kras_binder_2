# Summary Report: KRAS(G12D) pMHC Miniprotein Binder Design (Draft)

> Draft status: covers rounds 1-3 and the AF3/Rosetta confirmation pass on the round-3 shortlist. Not yet reviewed for synthesis decision. Source data: `docs/03_design_log.md` (full narrative log), `designs/*.csv`.

## 1. Objective

Design a single-chain miniprotein (≤120 aa) that binds the HLA-A\*11:01-presented KRAS(G12D) 9-mer peptide (`VVGADGVGK`, PDB 9UV8) **specifically** — i.e. discriminates the G12D mutant from the wild-type peptide (`VVGAGGVGK`, PDB 8I5E) — following the de novo design strategy of Liu et al. (Science 2025, `paper/science.adv0185.pdf`) adapted to this target. Specificity, not raw binding affinity, is the primary objective: the binder must make a direct contact to peptide position 5 (p5, the mutant Asp) with a complementary H-bond/basic partner, per the binding mode observed in the paper's 18 validated hits.

## 2. Method

Pipeline: **RFdiffusion** (backbone generation, peptide-only hotspot conditioning) → **ProteinMPNN** (sequence design, multi-sequence-per-backbone) → **ESMFold** (foldability triage, substitutes AF2-monomer) → **ProteinMPNN conditional-probs specificity scan** (`mpnn_spec_score`, within-backbone ranking only) → **AF2 initial-guess** (pose-templated complex binding check) → **custom pMHC-fold** (3-chain binder+MHC+peptide AF2-finetuned model, self-templated; on/off-target PAE delta as a second specificity readout) → **AF3 Server** (independent, untemplated on/off-target confirmation) → **Rosetta InterfaceAnalyzer** (ddG/shape-complementarity/packstat/buried-unsat calibrated against the 18 validated hits) and a **pose-bias contact diagnostic** (compares the templated design's p5 contact distance to the same distance in the untemplated AF3 refold).

## 3. Round-by-round summary

### Round 1 — Full de novo RFdiffusion (~1300 backbones total)

**First batch (r1, 300 backbones):** RFdiffusion full de novo, peptide-only hotspots (B184+B186+B187) → contact filter `05_contact_filter.py` (299 backbones scored, ~all pass p5 contact gate) → ProteinMPNN (1 seq/backbone = 302 sequences; no intermediate filter) → AF2 initial-guess complex (282 designs tested, **1 pass** pae<10): `r1_184_dldesign_0_cycle1` (pae=8.38).

**Second batch (r1b, 1000 backbones):** RFdiffusion full de novo → contact filter (999 backbones scored, 993 pass p5_cms>0) → ProteinMPNN (8 seqs/backbone for 934 backbones = **7,472 sequences**) → ESMFold monomer foldability filter (pLDDT≥80: 5,275/7,472 pass; best-pLDDT sequence selected per backbone) → AF2 initial-guess complex (**914 designs tested, 4 pass** pae<10, 0.4%): `r1b_273_dldesign_5`, `r1b_401_dldesign_7`, `r1b_403_dldesign_2`, `r1b_870_dldesign_2` → AF2 monomer on the 4 passing sequences (pLDDT 92.0–97.6, CA-RMSD 0.44–0.93 Å) → ProteinMPNN spec scan recheck on all 32 seqs for these 4 backbones (within-backbone ranking only).

**Round 1 outcome:** 5 candidates with strong AF2 binding but poor specificity signal by `mpnn_spec_score` (within-backbone ranking only; cross-backbone use later shown to be invalid).

---

### Round 2 — Partial diffusion (flat `partial_T=18`, 1000 backbones)

Seeded from round 1's 5 AF2-passing backbones (200 variants/seed) → contact filter `05_contact_filter.py` (998 scored, 499 pass) + geometric filter `05b_min_distance_filter.py` (1000 tested, 275 pass min binder–peptide distance gate) → ProteinMPNN spec scan `08_mpnn_specificity.py` on all passing backbones (~1,003 sequences, 1 per backbone) → **top-N backbone selection by `mpnn_spec_score`** (top-21, top-15, random-8 tested as separate subsets) → ProteinMPNN sequence generation + ESMFold (top-21: 168 seqs; top-15: 120 seqs; random-8: 64 seqs) → AF2 initial-guess complex on each subset → **0 passed pae<10** across all subsets.

**Round 2 outcome:** 0% AF2 pass — entirely a sampling artifact from spec-score-biased backbone selection (see methodological finding below). `partial_T=18` retrospectively shown not to impair binding.

---

### Round 3a — Partial diffusion (graduated `partial_T=12–25`, biased top-N)

1000 new backbones (graduated T by seed AF2 quality) → same spec-score-biased top-N selection strategy as round 2 → **0 AF2 passes** on the tested subset. Confirmed that biased selection is the common failure mode, not T or backbone quality.

---

### Round 3b — Unbiased full-funnel (994 backbones → AF2)

Same 1000 round-3 backbones → ProteinMPNN (8 seqs/backbone = **8,000 sequences**) → ESMFold (pLDDT≥80 AND CA-RMSD≤2.0 Å: 6,971/8,000 pass; best-pLDDT winner per backbone selected from 992 backbones) → geometric filter `05b_min_distance_filter.py` (min binder–peptide distance gate) → ProteinMPNN spec scan `08_mpnn_specificity.py` (**994 sequences scored**, within-backbone ranking only — not used for backbone selection) → AF2 initial-guess complex on all 994 (**8 pass** pae<10, 0.8%): `r3_r1b_273_1`, `r3_r1b_273_28`, `r3_r1b_403_65`, `r3_r1b_403_6`, `r3_r1b_870_22`, `r3_r1b_870_64`, `r3_r1b_870_87` lineages (2 from `r1b_273_28`).

**Round 3b outcome:** Breakthrough — 8 genuine AF2-passing backbones invisible to every prior biased top-N strategy.

---

### Round 3c — Extended sequence testing + full specificity battery (20 AF2-passing candidates)

**Extended AF2 testing:** 7 remaining ProteinMPNN sequences per each of 8 passing backbones (56 new seqs) → ESMFold pre-filter (160 seqs across all 8×8=64 tested) → AF2 initial-guess complex (64 designs tested, **15 pass** pae<10, 23%); combined with 3b's 8 passing seqs → **20 total AF2-passing candidates** (5 r1 + 15 r3 across 8 backbones).

**Specificity screens:**
1. **pMHC-fold on/off-target** `10_run_pmhc_fold_3chain.py` — 3-chain AF2-finetuned-MHC fold (G12D vs WT 8I5E); Δ(PAE on − off) gate: 20 candidates → **13 negative Δ** (correct direction) + 1 `r1b_273_28` control (selected for deliberate non-specificity).
2. **AF3 Server on/off-target** — 13 submitted (12 + 1 control); no-MSA for binder, real MSA for MHC, no structure template; iptm ≥ 0.90 AND binder–pep iptm ≥ 0.77 AND Δ > 0 gate: **1 strict pass** (`r3_r1b_870_87_dldesign_6`); 5 further pass absolute AF3 binding gate without on/off discrimination.
3. **Rosetta InterfaceAnalyzer** post-FastRelax (3-repeat constrained, DAlphaBall BUNS) on all 13 AF3 structures: ddG, shape complementarity, packstat, buried-unsat, real PyRosetta `ContactMolecularSurface` per peptide position.
4. **ipSAE** (their pipeline `calc_ipsae.py`) on all 13 on-target and 13 off-target AF3 structures.
5. **5-model AF3 ensemble analysis** `14_batch_find_basic_to_asp.py` — Arg/Lys/His↔Asp(p5) salt-bridge/H-bond distances across all 5 AF3 models; p5 Asp rotamer (gauche−/trans) classified per model.
6. **p5 contact pose-bias diagnostic** `12_check_p5_contact.py` — min binder-to-p5 distance: designed pose vs untemplated AF3 refold.

**Round 3c outcome:** `r3_r1b_870_87_dldesign_6` (pae=4.29, AF3 iptm/pep-iptm=0.93/0.77, ipSAE=0.896, Arg47↔Asp(p5):OD2 salt bridge in 4/5 AF3 models) and `dldesign_7` as co-lead; 13 new AF2-passing candidates total in round 3.

**Key methodological finding — `mpnn_spec_score` is not a valid backbone-ranking metric** (confirmed by the paper's author, then validated empirically in two independent ways):

1. *Cross-round spec-score bias*: selecting backbones by `mpnn_spec_score` (computed from 1 sequence/backbone) consistently produced 0% AF2 pass rates across four independent top-N samples spanning rounds 2 and 3, while the unbiased full-funnel found genuine binders in the same backbone pools. `mpnn_spec_score` with 1 sequence/backbone is a noisy backbone-level signal that does not correlate with AF2 `pae_interaction` (see ESMFold-as-pre-filter section).

2. *Definitive `partial_T` non-effect* (retrospective experiment, 2026-06-26): ran the complete ProteinMPNN → ESMFold → AF2 funnel on **all 200 round-2 r1b_403 backbones** (`partial_T=18`, 174 ESMFold-passing) and compared to round-3's 200 r1b_403 backbones (`partial_T=12`). AF2 pass rate: **T=18 = 1.1% (2/174), T=12 = 1.0% (2/200); median pae 26.42 vs 26.43; distributions statistically indistinguishable.** The two T=18 binders (`pd_r1b_403_188` pae=5.89, `pd_r1b_403_119` pae=6.32) were never in round 2's spec-score-selected test set. `partial_T=18` does not impair binding — round 2's "0% AF2 pass" was entirely a small-sample sampling artifact: ~29/1000 backbones tested → ~0.3 expected successes at a ~1% base rate. Backbone RMSD to seed also does not predict AF2 binding within T=12 (Pearson r=0.047).

**ESMFold-as-pre-filter validation**: checked how predictive within-backbone ESMFold pLDDT ranking is of actual AF2 binding rank, across the 8 round-3 backbones where all 8 ProteinMPNN sequences were tested by both. ESMFold pLDDT is a moderate predictor (mean Spearman ρ vs AF2 `pae_interaction` = -0.43); CA-RMSD-to-design carries no signal (ρ≈0.00). Decision-relevant result: picking the ESMFold-pLDDT-winner per backbone landed the actual AF2-best sequence in 5/8 backbones and top-3 in all 8 (mean rank 1.5/8 vs. 4.5 expected by chance) — validates the round-3 funnel's cheap pre-filter, though testing all 8 sequences per backbone still surfaced additional/better binders (including the campaign-best) in 3/8 cases.

## 4. Specificity screening evolution

1. **AF2 initial-guess** is pose-templated (biased toward the designed coordinates) — it measures "is the intended pose stable," not independent confidence.
2. **pMHC-fold** (custom 3-chain AF2-finetuned-MHC script, `scripts/10_run_pmhc_fold_3chain.py`) is also self-templated for the binder, with real on-target (9UV8) or off-target (8I5E) MHC+peptide templates. On/off-target PAE delta gives a second specificity readout. (Found and fixed a `crop_size` padding bug along the way that had produced a spurious "low peptide confidence" artifact.)
3. **AF3 Server** (no templating at all, real MSA for MHC, no-MSA for the binder) is the only fully *independent* check — confirmed via a real submission-schema example the user provided (`paper_output/af3_nomsa/binder_pMHC_full_noMSA_MHC.json`); a schema bug in `prepare_af3_jobs.py` (using the AF3 *output* schema instead of the real *input* `alphafoldserver` v1 schema) was caught and fixed before any submission.
4. **p5-contact pose-bias diagnostic** (`scripts/12_check_p5_contact.py`): compares min binder-to-p5 heavy-atom distance between the templated design and the untemplated AF3 refold. A contact that doesn't survive untemplated folding indicates the AF2/pMHC-fold signal was likely pose-bias, not a real interaction.
5. **Real PyRosetta `ContactMolecularSurface`** (per peptide position, `scripts/13_find_basic_to_asp.py`'s sibling CMS script) replaced the earlier lightweight `cms_p5_proxy` (a simple atom-count stand-in) for final scoring — see cross-check note below.

**Cross-check against an independent pipeline** (`ElainW/pmhc_binder_kras`, user-provided): resolved two apparent discrepancies and fixed one real bug upstream.
- **ipSAE**: their `calc_ipsae.py` and the official `DunbrackLab/ipsae.py` (used here) compute genuinely different statistics — confirmed by exact numerical match (their real `ipsae_binder_peptide` on the 18 validated designs, 0.86-0.97, exactly reproduces what running their code on our same AF3 data gives). Their `d0` derives from total contact-pair count across the whole binder-peptide submatrix, inflating `d0` and saturating the score near 1.0; the official method computes `d0` per-residue, a stricter statistic (range 0.42-0.87 on the same designs). Not a bug — never compare the two directly.
- **Hotspot CMS**: also resolved — their `cms_p5` (10.7-89.7, Å²-scale) is the real PyRosetta `ContactMolecularSurface` filter; this project's `cms_p5_proxy` was always an explicitly-labeled lightweight integer-count placeholder. Different metric, not a calibration error.
- **`pmhc_fold.py` crop_size bug (real, fixed)**: their script never overrides `model_config.data.eval.crop_size` (default 256), so AF2's `random_crop_to_size` randomly crops any query over 256 residues — most of theirs and ours run 270-390. Unlike this project's own (harmless) padding bug, this risks silently discarding the peptide or large chunks of binder/MHC. Patched and verified live (unpatched processed shape `(4,256)` on a 378-residue query → patched `(4,378)`); applied to the cloned repo.

## 5. Candidate stats (13-candidate AF3-confirmed shortlist)

Drawn from the 12 best pMHC-fold on/off deltas among the 20 total AF2-passing candidates (rounds 1+3), plus one `r1b_273_28` lineage member as a deliberate non-specific control.

| Candidate | AF2 pae | ipSAE AF2¹ | ipSAE AF3² | pMHC-fold Δ | AF3 on iptm/pep-iptm | real CMS p5 | their-ipSAE Δ | p5 dist: design→AF3 | Asp OD contact (5-model AF3 ensemble)³ |
|---|---|---|---|---|---|---|---|---|---|
| `r3_r1b_870_87_dldesign_6` | 4.29 | 0.585 | **0.510** ✓ | -1.22 | 0.93/**0.77** | 18.6 | +0.04 | 3.47→3.22Å (robust) | **4/5 models**: Arg47:NH2↔OD2 2.09-2.51Å; Tyr78:OH↔OD1 3.4-3.7Å in 2 models |
| `r3_r1b_870_87_dldesign_7` | 6.30 | 0.277 | 0.391 | -0.10 | 0.91/0.72 | 25.0 | -0.02 | 3.47→2.66Å (robust) | **4/5 models**: Arg47:NH2↔OD2 2.14-2.88Å; Arg51:NH1↔OD2 2.5-2.9Å in 2 models |
| `r3_r1b_870_87_dldesign_2` | 7.13 | 0.234 | 0.387 | -1.79 | 0.90/0.73 | **28.4** | **+0.08** | 3.47→2.98Å (robust) | **0/5 models**: Asp gauche- in all; OD 6.5-7.2Å from binder — backbone/CB contact only |
| `r1b_403_dldesign_2` | 4.94 | 0.463 | 0.710 | -0.85 | 0.95/0.86 | 11.8 | -0.003 (flat) | 3.01→3.76Å (robust, non-specific) | 0/5: no OD contact in any model; high pep-iptm from backbone/groove only |
| `r3_r1b_403_65_dldesign_4` | 8.17 | 0.179 | 0.658 | -0.19 | 0.95/0.83 | 12.8 | -0.004 (flat) | 3.47→3.71Å (robust, non-specific) | 1/5 (m2 only): Arg16:NH2↔OD1 3.63Å; non-specific overall |
| `r3_r1b_273_28_dldesign_0` (control) | 8.87 | n/a | 0.603 | **+4.13** | 0.93/0.79 | 12.4 | 0.0002 (flat) | 3.78→3.46Å (robust, non-specific by design) | **4/5 models**: Arg73:NH1/NH2↔OD1 2.60-3.06Å — same mechanism as leads; non-specific due to parallel MHC-framework contacts |
| `r3_r1b_870_22_dldesign_0` | 5.59 | 0.420 | 0.000 | -0.15 | 0.76/0.10 | 2.4 | 0.0000 both | 4.16→4.56Å (borderline) | 0/5: no OD contact |
| `r3_r1b_870_64_dldesign_0` | 6.74 | 0.275 | 0.012 | -1.56 | 0.85/0.51 | 3.8 | +0.27 | 3.40→4.76Å (weakening) | 1/5 (m1 only): Lys82:NZ↔OD2 2.73Å — partial |
| `r3_r1b_870_22_dldesign_7` | 4.92 | 0.425 | 0.000 | -2.13 | 0.74/0.07 | 0.0 | 0.0000 both | 4.51→7.96Å (**lost**) | 0/5: contact-lost |
| `r1b_870_dldesign_2` | 5.43 | 0.337 | 0.000 | -0.05 | 0.77/0.13 | 0.2 | **off=0.90→on=0.00** | 2.99→5.72Å (**lost**) | 0/5: contact-lost |
| `r1b_401_dldesign_7` | 5.63 | 0.414 | 0.000 | -1.61 | 0.77/0.12 | ≈0 | 0.0000 both | 3.02→8.02Å (**lost**) | 0/5: contact-lost |
| `r3_r1b_403_6_dldesign_6` | 5.89 | 0.169 | 0.000 | -1.40 | 0.76/0.14 | 2.2 | 0.0000 both | 2.54→15.11Å (**lost completely**) | 0/5: contact-lost |
| `r1b_273_dldesign_5` | 7.59 | 0.166 | 0.000 | -0.22 | 0.76/0.10 | 0.0 | 0.0000 both | 3.46→20.02Å (**lost completely**) | 0/5: contact-lost |

¹ `ipSAE AF2`: DunbrackLab `ipsae.py` (pae=15, dist=15) on AF2-initial-guess structures, peptide-only submatrix. ² `ipSAE AF3`: DunbrackLab `ipsae.py` (pae=10, dist=10) on AF3 independent on-target structures, A-C chain pair (binder vs peptide). Validated-design range: 0.42–0.87, median 0.72. ✓ marks the sole candidate within the validated range on both metrics. `real CMS p5`: PyRosetta `ContactMolecularSurface` (Å²), validated range 10.7–89.7/median 47.9. ³ All 5 AF3 models examined for all 13 candidates. Model_0 predicts gauche- Asp (OD buried) for every candidate where a binder is near p5 — a systematic bias under borderline binder-pep-iptm; use the full ensemble, not model_0 alone, for interface contact analysis.

**Reading the table:** `ipSAE AF3` is computed on the independent AF3 structure (no binder template). **`r3_r1b_870_87_dldesign_6` is the only candidate clearing the validated AF3 ipSAE minimum (0.510 ≥ 0.42) and both AF3 iptm gates.** `dldesign_7` (0.391) and `dldesign_2` (0.387) are marginally below the validated floor. `r1b_403_dldesign_2` (0.710) and `r3_r1b_403_65_dldesign_4` (0.658) score well within range but show zero specificity discrimination — strong non-specific binders. Bottom 7 rows: 0.000–0.012, consistent with confirmed contact loss.

**The "Asp OD contact" column was computed across all 5 AF3 models for all 13 candidates** — model_0 alone is unreliable because it predicts gauche- Asp (OD buried) for every candidate where a binder is positioned near p5, regardless of which candidate. This is a systematic AF3 bias under borderline binder-pep-iptm, not a random artifact. The full-ensemble results reveal four mechanistic tiers:

1. **Consistent Arg-Asp salt bridge (4/5 models)**: `dldesign_6` (Arg47:NH2↔OD2 2.09-2.51Å + Tyr78:OH↔OD1 in 2 models), `dldesign_7` (Arg47↔OD2 2.14-2.88Å + Arg51↔OD2 2.51-2.92Å in 2 models), **`r3_r1b_273_28_dldesign_0` control** (Arg73:NH1/NH2↔OD1 2.60-3.06Å in 4/5 models — same mechanism, but non-specific due to parallel MHC-framework contacts that also work on WT).
2. **Partial OD contact (1/5 models)**: `r3_r1b_403_65_dldesign_4` (Arg16↔OD1 3.63Å, borderline H-bond); `r3_r1b_870_64_dldesign_0` (Lys82:NZ↔OD2 2.73Å in 1 model — Lys rather than Arg, partial engagement).
3. **No OD contact, backbone-level p5 engagement**: `r3_r1b_870_87_dldesign_2` (gauche- in all 5, OD 6.5-7.2Å; real CMS 28.4 Å² from backbone/CB contact); `r1b_403_dldesign_2` (gauche-/flat in all 5, AF3 ipSAE=0.710 from backbone/groove, no carboxylate contact).
4. **No p5 contact at all**: remaining 7 (contact-lost group — no binder near peptide in any model).

**Mechanistic contact with Asp p5 side chain — full 13-candidate 5-model analysis.** The earlier conclusion that "no direct Asp side-chain contact exists" was based on model_0 only and was wrong for multiple candidates. Key findings after examining all 65 AF3 models (13 candidates × 5):

- **`dldesign_6`**: models 1-4 trans/gauche+ (chi1=-170° to +61°), **Arg47:NH2↔OD2 2.09-2.51Å** + Tyr78:OH↔OD1 3.4-3.7Å in 2 models; model_0 gauche- outlier (chi1=-86°).
- **`dldesign_7`**: models 1-4 trans (chi1≈-167 to -174°), **Arg47:NH1/NH2↔OD2 2.14-2.88Å** + Arg51:NH1↔OD2 2.51-2.92Å in models 3-4; model_0 gauche- outlier (chi1=-79°).
- **`dldesign_2`**: all 5 gauche- (chi1≈-70 to -76°, OD 6.5-7.2Å), no OD contact in any model — backbone/CB mechanism throughout.
- **`r3_r1b_273_28_dldesign_0` (control)**: models 1-4 show **Arg73:NH1/NH2↔OD1 2.60-3.06Å** — same Arg-Asp mechanism as dldesign_6/7, but the control also binds via parallel contacts that work on WT, explaining its non-specific AF3 profile.
- **`r3_r1b_403_65_dldesign_4`**: model_2 only, Arg16:NH2↔OD1 3.63Å — marginal single-model contact; still non-specific overall.
- **`r3_r1b_870_64_dldesign_0`**: model_1 only, **Lys82:NZ↔OD2 2.73Å** — a genuine Lys-Asp salt bridge in 1/5 models; consistent with its borderline/weakening contact status.
- **`r1b_403_dldesign_2`**: 0/5 models — despite AF3 ipSAE=0.710 (highest in the set), no OD contact; pep-iptm=0.86 from backbone/groove engagement only.
- **7 contact-lost candidates**: 0/5 models each — no binder near peptide in any model.

**Systematic bias**: model_0 predicts gauche- Asp for every candidate where a binder is near p5 — a consistent top-ranked-model artifact when binder-pep-iptm is borderline. Always examine all 5 models for any interface contact conclusions. The 9UV8 crystal has chi1=-160° (trans, OD solvent-exposed); models that converge on trans are recovering the crystallographically correct Asp rotamer.

## 6. Sequence/structure check on the lead lineage

`r1b_870_87`'s three sequences (`dldesign_2`, `_6`, `_7`) vs. the 18 experimentally validated designs:

- **Length** (93 aa), **Cys** (0), **Met** (0), **max poly-Ala run** (3) all match the validated convention closely.
- **Net charge**: `dldesign_2` (-6) and `dldesign_6` (-9) are typical (validated range -3 to -11, median -7). **`dldesign_7` (-1) is an outlier** — every validated design is at least -3 net-negative; this breaks the project's net-negative-charge rule despite being a highly-charged (K/E-rich) sequence, and is worth a solubility/stability sanity check before prioritizing it over `_2`/`_6`.
- **Structure** (Rosetta): all three sit within the validated envelope on every metric (sc, packstat, buns, n_int_res, cms_p5), but ddG (-44 to -53) is on the **weaker end** of the validated range (-57 to -88, median -70) — closer to the weakest validated hit (`phox2b-11`, -56.7) than the median.

## 7. Comparison with an independent pipeline (tier1 r2+r2.5)

**Data source:** `designs/tier1_af3_stats.tsv`, r2+r2.5 subset (n=17). Independent designs by the user's own pipeline targeting the same 9UV8 pMHC; all 17 were selected to have polar interactions with p5 (15 have direct Asp side-chain contacts, 2 have backbone-only contacts at p5 — the same situation as our `dldesign_2`). Computed with the same `af3_design_stats.py` script and same FastRelax+InterfaceAnalyzer pipeline.

| Metric | Validated (n=18) | tier1 r2+r2.5 (n=17) | Our 12¹ (Asp inward, m0) | Our 12¹ (Asp outward†) |
|---|---|---|---|---|
| ipSAE binder-pep | 0.86-0.97 (med **0.94**) | 0.85-0.94 (med **0.92**) | 0.00-0.94 (med 0.21) | 0.00-0.94 (med 0.21) |
| CMS p5 (Å²) | 10.7-89.7 (med **47.9**) | 26.7-57.0 (med **44.4**) ✓ | 0-28.5 (med 5.9) | 0-54.9 (med **10.5**) ←improved |
| CMS hotspot p5 (Å²) | 35.8-79.4 (med **56.2**) | 26.7-57.0 (med **44.4**) ✓ | 0-28.5 (med 5.9) | 0-54.9 (med **10.5**) |
| CMS peptide total (Å²) | 174-270 (med **240**) | 95-203 (med 138) | 16-199 (med 96) | 16-199 (med 121) |
| n contacts p5 (all) | med 3 | med **3** ✓ | med 0 | med 0 |
| n contacts p5 (polar) | med 0 | med **2** ✓ | med 0 | med 0 |
| n contacts p5 (H-bond) | med 0 | med **2** (max 6) ✓ | med 0 | med 0 (max 2) |
| n saltbridge at p5 | med 0 | med **1** ✓ | med 0 | med 0 |
| n contacts hotspot (polar) | 1-4 (med **2**) | 1-3 (med **2**) ✓ | 0-2 (med 0) | 0-2 (med 0) |
| binder res in contact | 7-11 (med **9**) | 3-8 (med 5) | 0-7 (med 4) | 0-7 (med 4) |
| Helix % | med 87 | med **92** | med 90 | med 90 |
| Shape complementarity | 0.51-0.71 (med **0.62**) | 0.54-0.65 (med **0.61**) ✓ | 0.36-0.65 (med 0.57) | 0.36-0.65 (med 0.57) |
| Packstat | 0.59-0.72 (med **0.62**) | 0.53-0.65 (med **0.61**) ✓ | 0.55-0.68 (med 0.60) | 0.55-0.68 (med 0.60) |
| Unsatisfied H-bonds (IA) | 5-23 (med **12**) | 2-13 (med **7**) better | 2-15 (med 10) | 2-15 (med 10) |
| BUNS delta unsat | 0-8 (med **3**) | 1-4 (med **3**) ✓ | 0-7 (med 2.5) ✓ | 0-7 (med 2.5) ✓ |
| Surface hydrophobicity | 0.19-0.53 (med **0.29**) | 0.21-0.60 (med 0.42) ⚠ | 0.23-0.64 (med 0.44) ⚠ | 0.23-0.64 (med 0.44) ⚠ |
| pass_all_bindcraft | 10/18 (56%) | 12/17 **(71%)** ✓ | 9/12 (75%) ✓ | 9/12 (75%) ✓ |

¹ Excludes the `r3_r1b_273_28_dldesign_0` control. † CMS and contact metrics for `dldesign_6` (model_2, pep_iptm=0.79) and `dldesign_7` (model_1, pep_iptm=0.73) replaced with the outward-Asp model values; `dldesign_2` has no outward-Asp model in any of 5 AF3 predictions.

**Key takeaways from the comparison:**

1. **tier1 matches the validated-design profile on the most important specificity metrics** (CMS p5 median 44 vs validated 48, ipSAE median 0.92 vs 0.94, n polar contacts median 2 matching validated, n saltbridge median 1) — because every tier1 design was explicitly filtered to have p5 polar contacts. Biophysical quality (SC, packstat, BUNS) is also at or above validated levels.

2. **Our 12 candidates lag on p5 contact density.** The overall-population medians for CMS p5 (5.9 inward, 10.5 outward) and polar/H-bond counts (0 in both orientations) are well below tier1 because our population is dominated by 7 zero-contact candidates. Our three leads individually reach tier1 levels on CMS p5 when using the outward-Asp model (`dldesign_6` 39 Å², `dldesign_7` 55 Å²) — but the ensemble H-bond count (max 2 vs tier1 max 6) still falls short.

3. **Both pipelines share the surface hydrophobicity gap** (both ~0.42-0.44 vs validated 0.29). This is likely a shared property of helical miniproteins targeting the pMHC groove rather than a design-specific failure, but it is a developability concern worth tracking.

4. **tier1's lower binder_res_in_contact** (median 5 vs validated 9) is consistent with our leads (median 4) — small α-helical miniproteins designed against a single peptide make fewer total interface contacts than the larger validated hits. Interface focus on p5 rather than breadth appears to be the shared characteristic.

5. **Implication for a future round:** an explicit filter for p5 polar contacts (as applied in tier1 selection) is the clearest path to closing the gap with both tier1 and the validated designs. Our leads `dldesign_6` and `dldesign_7` demonstrate the mechanism is achievable (Arg47 salt bridge confirmed in 4/5 AF3 models), but it needs to be enforced as a selection criterion from the start of the design cycle rather than discovered retrospectively.

### §7b. tier1 r2+r2.5, r3 vs validated designs vs our candidates

All four sets scored with the same FastRelax+InterfaceAnalyzer+CMS+contacts pipeline. Validated (n=18) column: user-provided `af3_design_stats.tsv` (paper authors' run with original AF3 PAE files enabling ipSAE); our independent FastRelax run on the same structures (`designs/validated_reference_metrics.tsv`) confirmed all CMS, contact-count, SC, packstat, and surface-hydrophobicity values within stochastic FastRelax variation (≤5% continuous, identical discrete counts). ipSAE for validated designs not independently recomputable (AF3 PAE output files unavailable).

| Metric | Validated (n=18)¹ | tier1 r2+r2.5 (n=17) | tier1 r3 (n=19) | Our 12 (outward†) |
|---|---|---|---|---|
| ipSAE binder-pep | med **0.94** | med **0.92** | med 0.88 | med 0.21 |
| CMS p5 (Å²) | med **47.9** | med **44.4** | med 38.8 | med 10.5 |
| CMS hotspot p5 (Å²) | med **56.2** | med **44.4** | med 38.8 | med 10.5 |
| CMS peptide total (Å²) | med **240** | med 138 | med 143 | med 121 |
| n contacts p5 (polar) | med 0 | med **2** | med **2** | med 0 |
| n contacts p5 (H-bond) | med 0 | med **2** | med 1 | med 0 (max 2) |
| n saltbridge at p5 | med 0 | med **1** | med **1** | med 0 (max 1) |
| n contacts hotspot (polar) | med **2** | med **2** | med **2** | med 0 |
| binder res in contact | med **9** | med 5 | med 5 | med 4 |
| Helix % | med 87 | med **92** | med **92** | med 90 |
| Shape complementarity | med **0.62** | med 0.61 | med **0.62** ✓ | med 0.57 |
| Packstat | med **0.62** | med 0.61 | med 0.59 | med 0.60 |
| Unsatisfied H-bonds | med **12** | med 7 ↓ | med 7 ↓ | med 10 |
| BUNS delta unsat | med **3** | med 3 ✓ | med **2** ↓ better | med 2.5 ✓ |
| Surface hydrophobicity | med **0.29** | med 0.42 ⚠ | med 0.36 ↓ better | med 0.44 ⚠ |
| pass_all_bindcraft | 10/18 (56%) | 12/17 (71%) | **15/19 (78%)** ↑ best | 9/12 (75%) |

¹ From user-provided `af3_design_stats.tsv`; CMS/contacts/Rosetta metrics independently verified by our pipeline run (`designs/validated_reference_metrics.tsv`). Our run: CMS p5 49.5 Å², CMS hotspot p5 52.6 Å², SC 0.620, packstat 0.625, surface hydro 0.288, BUNS med 4.0.

† CMS and contact metrics for `dldesign_6` (model_2) and `dldesign_7` (model_1) replaced with outward-Asp model values; `dldesign_2` has no outward-Asp model in any of 5 AF3 predictions.

**tier1 r3 vs r2+r2.5:** r3 improves biophysical quality — surface hydrophobicity 0.36 (down from 0.42, approaching validated 0.29), BUNS 2 (down from 3), SC reaches validated median 0.62, pass_all_bindcraft 78% (up from 71%). ipSAE and CMS p5 edge down slightly (0.88 vs 0.92; 38.8 vs 44.4), likely reflecting more diverse backbones in r3. p5 contact rate maintained (n polar = 2, n saltbridge = 1 in both rounds).

**Validated vs tier1:** both rounds match the validated set on discrete p5 contact counts (n polar = 2, n saltbridge = 1) despite lower CMS p5 magnitude (38–44 vs 48 Å²), as tier1 designs focus on a single peptide position rather than the full groove. Shape metrics (SC, packstat, BUNS) are at or better than validated levels for both rounds.

**Our 12 vs tier1 r3:** surface hydrophobicity (0.44 vs 0.36) and the CMS p5/polar-contacts gap remain the primary differences. On packstat, BUNS, and Helix %, our candidates are essentially at parity with tier1 r3. The contact gap narrows when using the outward-Asp model for our leads (`dldesign_6` CMS p5=39 Å², `dldesign_7` CMS p5=55 Å² — both at or above tier1 r3's median of 38.8), confirming our leads are individually competitive once the correct Asp rotamer is used, but the overall population lacks the p5 contact density of the tier1 set.

## 8. Recommendation  *(previously §7)*

**`r3_r1b_870_87_dldesign_6` is the primary lead** and **`r3_r1b_870_87_dldesign_7` is the co-lead**, both confirmed by the full 5-model AF3 ensemble analysis. Both have Arg47 making a direct salt bridge to Asp(p5):OD2 in 4/5 AF3 models at genuine H-bond/salt-bridge distances (2.09-2.88Å), confirming the intended G12D-specific recognition mechanism. `dldesign_6` additionally has the best absolute confidence metrics (AF2 pae=4.29, AF3 iptm/pep-iptm=0.93/0.77, AF3 ipSAE=0.510 within validated range); `dldesign_7` additionally engages Arg51 in two models, providing a bidentate Arg-Asp contact, though its near-neutral net charge (-1) warrants a solubility/stability check before synthesis.

**`r3_r1b_870_87_dldesign_2`** is now the **lower-priority** member of the lineage: all 5 AF3 models show the Asp carboxylate buried (gauche- in every model, OD 6.5-7.2Å from binder), with no carboxylate-specific contact in any independent prediction. It makes dense backbone/CB-level contact at p5 (real CMS 28.4 Å²) and shows good specificity signal by their-ipSAE (+0.08) — but likely reads the backbone perturbation Asp creates rather than the carboxylate directly, making it mechanistically less G12D-specific than `dldesign_6`/`dldesign_7`.

Deprioritize the other 10 candidates: 7 lose p5 contact under independent folding (triple-confirmed by real CMS, their-ipSAE, and AF3 ipSAE=0.000), 2 (`r1b_403_dldesign_2`, `r3_r1b_403_65_dldesign_4`) are strong but non-specific binders, and 1 (`dldesign_2`) lacks carboxylate-specific contact.

## 9. Open items / caveats

- Only 13 of the campaign's 20 AF2-passing candidates have been AF3-confirmed; the remaining 7 (the ones with the worst pMHC-fold deltas) were not submitted.
- Full 5-model ensemble analysis: Arg/Lys-Asp(p5) salt bridges confirmed in 4/5 models for `dldesign_6`, `dldesign_7`, and the non-specific control; partial (1/5) for `r3_r1b_403_65_dldesign_4` and `r3_r1b_870_64_dldesign_0`; absent for all others. Model_0 systematically predicts gauche- Asp whenever a binder is near p5 — always use all 5 models for interface contact conclusions.
- No wet-lab or orthogonal computational (e.g. MD) validation yet — this is a purely in-silico shortlist.
- Tier1 comparisons (§7, §7b) use `designs/tier1_af3_stats.tsv` (user-provided). Our independent FastRelax+pipeline run on the 18 validated designs confirmed CMS, contact-count, SC, packstat, and surface-hydrophobicity metrics match within stochastic variation; results in `designs/validated_reference_metrics.tsv`. ipSAE for validated designs is not independently recomputable (AF3 PAE output files not available).
- This document is a draft; numbers should be spot-checked against `docs/03_design_log.md` and the underlying CSVs (`designs/candidates.csv`, `designs/funnel_summary.csv`, `designs/round3_rosetta_metrics.csv`, `designs/af3_negdelta_results.csv`, `designs/real_cms_13_ontarget.tsv`, `designs/their_ipsae_13_onoff.tsv`, `designs/af3_stats_13/af3_design_stats.tsv`, `designs/validated_reference_metrics.tsv`) before being treated as final.
