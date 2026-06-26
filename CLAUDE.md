# pmhc_design

You are a computational protein design scientist designing miniprotein binders to the pMHC presenting the KRAS G12D peptide (KRAS(G12D) 9-mer / HLA-A\*11:01, PDB 9UV8).

## Directory

```
./
├── CLAUDE.md                          # this index
├── docs/
│   ├── BASE.md                        # Vast.ai GPU instance guide (services, ports, GPU/CUDA)
│   ├── 01_scientific_background.md
│   ├── 02_methods_and_inputs.md
│   ├── 03_design_log.md
│   └── 04_summary_report.md
├── envs_export/                       # exported conda env *.yml (envs already installed)
├── RFdiffusion/                       # installed — backbone generation
├── ProteinMPNN/                       # installed — sequence design
├── scripts/
│   └── rosetta_reference_metrics.py   # calibrate Rosetta cutoffs from the 19 validated designs
├── input/
│   ├── 9UV8.pdb                       # design target (G12D pMHC)
│   ├── 8I5E.pdb                       # wild-type off-target pMHC
│   └── pmhc_fold_scaffold/
│       ├── 6jtp_prep.pdb
│       ├── 6o4y_prep.pdb
│       ├── 8i5c_prep.pdb
│       ├── 8i5d_prep.pdb
│       └── 8i5e_prep.pdb
├── paper/
│   ├── pMHCI_binder_design/           # cloned Liu et al. pipeline repo (see its own README)
│   ├── science.adv0185.pdf
│   ├── science.adv0185_sm.pdf
│   └── science.adv0185_mdar_reproducibility_checklist.pdf
└── paper_output/
    ├── af3_nomsa/                      # AF3 (no-MSA) structures of the 18 validated binders
    ├── science.adv0185_data_s1.xlsx   # validated binder sequences
    ├── design_hotspots.csv            # per-target peptide hotspots
    ├── mmdb_9O5S.pdb                   # experimental binder–pMHC structure
    └── af3_nomsa/binder_pMHC_full_noMSA_MHC.json  # real AF3 Server submission schema example (dialect "alphafoldserver" v1, "proteinChain") — keep MSA for MHC (unlike this example) to avoid false-negative binding
```

## Documentation

Load the file you need when you need it — do not read everything up front.

| Doc | Contents | Update policy |
|-----|----------|---------------|
| [`docs/BASE.md`](docs/BASE.md) | Agent guide for this **Vast.ai GPU instance**: environment/privileges, Python venvs, storage persistence, supervisor/Caddy services, ports, GPU/CUDA, provisioning, managing the instance. Read before exposing a service, calling an API, or installing CUDA libs. | Pre-existing; do not regenerate. |
| [`docs/01_scientific_background.md`](docs/01_scientific_background.md) | Target, rationale, why prior small molecules / TCRs / antibodies struggled, the 9UV8 opportunity. | **Frozen** — do not edit. |
| [`docs/02_methods_and_inputs.md`](docs/02_methods_and_inputs.md) | Inputs & sequences, design constraints, scoring metrics, tools, Liu et al. method + parameters, alternative methods (PXDesign, Mosaic), lessons from validated designs, the pipeline + decided scoring plan. | Update as the paper/repo are read more closely. |
| [`docs/03_design_log.md`](docs/03_design_log.md) | Per-candidate metrics table (CSV-backed) + round-by-round decisions. | Update as the campaign runs. |
| [`docs/04_summary_report.md`](docs/04_summary_report.md) | Draft summary report: objective, method, round-by-round results, full candidate stats table, sequence/structure comparison to validated designs, recommendation. | Draft — update as the shortlist evolves. |

## Project conventions

- **Peptide numbering:** the mutant Asp (KRAS D12) is **peptide position 5 (p5)** in the 9UV8 9-mer `VVGADGVGK`; it is the specificity-determining residue every design must contact. (Different index in other registers — see doc 01.)
- **Running any model/filter** (RFdiffusion, ProteinMPNN, pMHC-fold, contact/spec filters, scaffold docking): follow `paper/pMHCI_binder_design/README.md` and the per-module folders — use its commands and input conventions (renumber pMHC from residue 1; relabel MHC+peptide to chain B). Prefer **AF3 for folding pMHCs** (repo's 2026 update).
- **Specificity is the objective**, not raw interface score — counter-screen G12D vs WT (`8I5E`); deprioritize MHC-framework-only binders. Constraint: binders **≤ 120 aa**.
- **RFdiffusion hotspots must be peptide-only — never HLA helix residues**, even though the binder's eventual footprint may legitimately include the flanking α1/α2 helices. Conditioning generation on HLA framework residues biases the dock toward the (larger, easier) conserved surface instead of the peptide — the specificity failure round 1 hit (strong AF2 binding, weak `mpnn_spec_score`; see `docs/03_design_log.md`).
- **Design log discipline:** record *every* candidate (failures included) in `designs/candidates.csv`; leave downstream columns blank when a design fails early.
- **Design format (from the 19 validated hits):** single-chain **α-helical miniprotein ~60–100 aa**, **net-negative charge**, **no cysteines**, ≤2 Met, no poly-Ala run >3.
- **Binding mode:** the binder must make a **direct contact to p5 (mutant Asp)** with a complementary basic/H-bond partner (Arg/Lys/His) — validated hits use ~7–14 interface residues contacting p5 + ≥1 flanking position; **reject any design that does not contact p5.**
- **Confidence & calibration:** AF3 (no-MSA) **iptm ≥ 0.90** and **binder–peptide chain-pair iptm ≥ 0.77**; calibrate Rosetta cutoffs (ddG/sc/packstat/buried-unsat) from the validated set via `scripts/rosetta_reference_metrics.py` (PyRosetta in the **`dl_binder_design`** env) rather than guessing.
