# Scientific Background — KRAS(G12D) 9-mer / HLA-A\*11:01 (PDB 9UV8)

> **Status: frozen.** This file is the scientific rationale for the project. Do not update it during the design campaign — working notes, methods, and results go in [`02_methods_and_inputs.md`](02_methods_and_inputs.md) and [`03_design_log.md`](03_design_log.md).

You are a computational protein design scientist who wants to design miniprotein binders to pMHC presenting the KRAS G12D peptide.

## Goal

Design de novo protein binders (minibinders / TCR-like binders) for the **KRAS(G12D) peptide–MHC class I complex** using computational workflows. Top designs will be selected for experimental validation.

Not intended to solve antigen-density or presentation challenges (e.g. low pMHC copy number per cell, HLA loss), but to demonstrate that modern de novo binder design workflows can generate high-affinity, **register- and mutation-specific** binders to a clinically validated public neoantigen where the dominant drug modalities (intracellular small molecules and degraders) and traditional antibody/TCR approaches have struggled.

## Target

**Target:** KRAS(G12D) neoantigen peptide presented on MHC class I
**Complex of interest:** KRAS(G12D) **9-mer (VVGADGVGK)** bound in the peptide-binding groove of **HLA-A\*11:01**
**Structure:** Classical MHC class I fold — α1/α2 peptide-binding platform (β-sheet floor flanked by two α-helices), closed at both ends, with the β2-microglobulin light chain. The mutant aspartate (D12) sits near the apex of the groove; the 9-mer adopts a **bent ("arched") conformation** distinct from the more compact 10-mer register.
**PDB reference:** **9UV8** — HLA-A\*11:01 in complex with the KRAS(G12D) 9-mer peptide (Zhu et al., *Communications Biology*, 2025).

The pMHC complex is a soluble, well-folded, refoldable ~45 kDa assembly (heavy chain + β2m + peptide); multiple high-resolution structures of KRAS G12-mutant peptides on HLA-A\*11:01 and HLA-C\*08:02 exist (e.g. 7OW4 / 7OW6 / 7PB2, 9UV8).

**Design-critical structural feature.** In 9UV8 the bent 9-mer presents the **mutant Asp12 outward and solvent-exposed at the groove apex**, where it forms an electrostatic contact with HLA **Thr73**, and the backbone bends away from the α2 helix toward α1. This places the single discriminating residue **directly on the targetable surface** — in contrast to the 10-mer register (on both HLA-A\*11:01 and HLA-C\*08:02) used by all known TCRs and the Therazyne antibody, in which the same aspartate is **buried in the groove** and read only indirectly. A designed binder to 9UV8 can therefore contact the mutation directly.

**Therapeutic rationale:** A binder to the KRAS(G12D) pMHC engages the mutant antigen **at the cell surface**, sidestepping the intracellular-delivery and shallow-pocket problems that constrain small molecules and targeted protein degraders. Because discrimination rests on the **single mutated, solvent-exposed side chain (Asp12)**, a sufficiently specific binder can be wired into T-cell–engaging formats (CAR/STAR/TruC, bispecifics) for off-the-shelf therapy against a shared driver mutation.

**Context on registers and known binders:**
- **9-mer register (VVGADGVGK, 9UV8):** bent conformation, Asp12 exposed, **incompatible with natural TCR recognition** — the target here. Unexplored by existing biologics.
- **10-mer register (VVVGADGVGK, HLA-A\*11:01):** compact, stable, Asp12 buried, TCR-accessible. Target of the Poole et al. therapeutic TCR (7OW6/7PB2) and the CorreGene TCRs (KDA11-01/02).
- **10-mer register (GADGVGKSAL, HLA-C\*08:02):** here the mutant Asp sits at peptide position 3, buried below the peptide plane and salt-bridging HLA Arg156. Target of NT-112 (TCR-T, NCT06218914) and the Therazyne TCR-like antibody TZ-Ab101 (Ahn et al., *Mol Ther* 2026).

> **Note on peptide numbering.** The mutant aspartate (KRAS residue D12) maps to a *different position index in each register*: **p5** in the 9UV8 9-mer (VVGA**D**GVGK), **p6** in the HLA-A\*11:01 10-mer (VVVGA**D**GVGK), and **p3** in the HLA-C\*08:02 10-mer (GA**D**GVGKSAL). Throughout this document "the 5th peptide residue" / "peptide position 5" refers to the mutant Asp **in the 9UV8 9-mer**, which is the design target.

## Why KRAS(G12D) pMHC?

KRAS G12D is the **single most common KRAS substitution** — ~40% of pancreatic ductal adenocarcinoma and ~5% of NSCLC — and historically undruggable: unlike G12C it has **no reactive cysteine and only a shallow switch-II pocket** (Park et al., *NEJM* 2026), and its aspartate is a **low-reactivity, surface-abundant** covalent handle (Weller et al., *Science* 2025). KRAS(G12D)/MHC is a **public neoantigen** — a shared, tumor-specific peptide from a recurrent driver displayed by common HLA alleles, broadly actionable across patient subsets. Clinical validation already exists: durable responses in pancreatic cancer patients treated with KRAS(G12D)/HLA-C\*08:02 TCR-T (Leidner et al., *NEJM* 2022; foundation: Tran et al., *NEJM* 2016), now productized as NT-112.

## Where existing approaches have struggled

**Intracellular small molecules / degraders** (daraxonrasib, zoldonrasib, setidegrasib) must reach RAS *inside* the cell and grip a target with **no reactive cysteine and only a shallow switch-II pocket** (Park et al., *NEJM* 2026), where the mutant aspartate is a **low-reactivity, surface-abundant** covalent handle (Weller et al., *Science* 2025). Daraxonrasib's multiselectivity adds WT-RAS toxicity; zoldonrasib needs an engineered neomorphic CYPA interface to make the Asp covalency work; setidegrasib is IV-weekly with infusion reactions and only modest PDAC activity. None exploits the surface-displayed antigen.

**TCR-T therapies** validate the pMHC axis but are biologically fragile. **NT-112** (Neogene/AstraZeneca; NCT06218914) is a triple-edited autologous product — knock-in of an HLA-C\*08:02-restricted KRAS(G12D) TCR plus **TRBC1/2 knockout** (prevent mispairing) and **TGFBR2 knockout** (resist TGF-β suppression in the TME) — engineering required precisely because native TCRs are low-affinity (µM-range) and shut down in solid tumors. Binding affinity is unpublished (functional reactivity only).

**TCRs on the relevant HLA-A\*11:01 allele recognize the mutation only indirectly.** In the structurally characterized complexes (Poole et al., *Nat Commun* 2022; 7OW4/7OW6), the mutant **Asp12 is buried in the groove**, salt-bridging HLA floor residues (Arg114/Arg65); the TCR's **CDR3 loops contact the bulged central peptide (positions 4–6)** while the germline CDR1/2 loops rest on the conserved HLA α-helices. Mutant-vs-WT discrimination comes not from a direct contact to the aspartate but from a **thin, enthalpy-driven difference in surface electrostatics** around the buried mutation. That indirectness leaks specificity: the best-characterized HLA-A\*11:01 TCRs (Zheng et al., 2024) reach only **functional EC50s of ~13–50 nM** (functional avidity, not biophysical KD) and one (KDA11-01) is **alloreactive to HLA-B\*57:01**. Every one of these — plus the high-affinity Poole TCR — targets the **10-mer** register.

**Antibody discovery against pMHC** faces the same core obstacle: **~85% of the binder's contact surface is conserved MHC, and the only therapeutically meaningful signal is one mutated peptide side chain.** The Therazyne TCR-like antibody (Ahn et al., *Mol Ther* 2026; crystallized lead Fab C8K10D5-1 / TZ-Ab101) actually solved the MHC-surface problem well — its 4.5 Å structure shows a **peptide-centric footprint** in which **CDR-H3, CDR-L1, and CDR-L3 cover essentially the entire peptide surface with only minor HLA contact**. But it still reads the mutation **indirectly**: the mutant **Asp3p is not contacted by the Fab** — it points below the peptide plane and salt-bridges **HLA Arg156**, acting as a presentation-stabilizing anchor. The antibody's true energetic hotspot is **Lys7p**, which makes ionic/H-bond contacts to **Asp31 and Tyr29 of the light chain** (flanked by heavy-chain Met106 and HLA Gln155); alanine scanning abolishes binding only at peptide positions **3, 7, and 10**. Reaching this still required iterative affinity maturation to low-nM KD, scrubbing **residual non-specific cell-surface binding** via mammalian cell-based negative selection, and remained restricted to the 10-mer on HLA-C\*08:02. Tellingly, an **AlphaFold2 model of this complex diverged significantly from the crystal structure** (low pLDDT over the peptide and CDR loops) — direct evidence that computational pMHC-interface predictions need experimental/orthogonal validation.

**The 9-mer register in 9UV8 is a recognition blind spot — and the opening.** Zhu et al. (2025) show this register is a **bent conformation incompatible with TCR recognition** (a 10-mer-selected TCR cannot bind it), yet it **displays the mutant Asp12 outward and solvent-exposed**. No binder — TCR or antibody — has been reported against it, and unlike every existing approach a designed binder could read the mutation **directly** rather than inferring it from buried-Asp electrostatics.

## The opportunity

**No de novo (mini)binder exists for the KRAS(G12D) 9-mer / HLA-A\*11:01 complex (9UV8).** The pMHC is a well-characterized, high-resolution, soluble, refoldable target with the discriminating mutant residue (Asp12) solvent-exposed at the groove apex — well suited to diffusion-based generative design and structure-guided docking pipelines. A de novo binder offers a smaller, stable, manufacturable scaffold that (1) targets the mutant antigen at the cell surface rather than the undruggable intracellular protein, (2) can exceed native-TCR affinity while preserving single-residue specificity, and (3) opens a register/allele that the antibody and TCR-T efforts to date have not addressed — finding a path the existing modalities have missed.
