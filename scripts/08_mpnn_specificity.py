#!/usr/bin/env python
"""Local ProteinMPNN log-probability specificity screen -- no AF3/pMHC-fold needed.

Per the supplement's "ProteinMPNN log probability screening" section (separate from the
AF-MHC/pMHC_fold screening): score the conditional probability ProteinMPNN assigns to the
on-target peptide residue (Asp/p5) given the rest of the complex (including the bound
binder) as fixed context, and filter on either the raw log-prob or the difference to the
alanine-mutant log-prob. Both come for free from a single conditional_probs_only pass
(it returns the full 21-AA distribution at the scored position), so no off-target
structure or extra AF3 job is needed for this signal.

For each input complex PDB (chain A = binder, fixed/context; chain B = target+peptide):
  1. parse_multiple_chains.py -> assign_fixed_chains.py (designed_chain=B) ->
     make_fixed_positions_dict.py (--specify_non_fixed, only the scored residue(s) free)
  2. protein_mpnn_run.py --conditional_probs_only 1
  3. read log P(Asp) and log P(Ala) at the scored position from the output .npz

mpnn_spec_score = log_p(on-target aa) - log_p(Ala) -- per the paper's difference filter.
Higher (less negative) = binder more strongly "reads" the on-target residue vs a neutral Ala.

Usage:
    python 08_mpnn_specificity.py --pdbdir <dir of complex pdbs, chain A=binder, chain B=target> \
        --out_csv mpnn_spec.csv --resnum 185 --chain B --on_target_aa D
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

import numpy as np

MPNN_HELPER = "/workspace/ProteinMPNN/helper_scripts"
MPNN_RUN = "/workspace/ProteinMPNN/protein_mpnn_run.py"
ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"


def run(cmd):
    subprocess.run(cmd, check=True)


def score_one(pdb_path, resnum, chain, work_dir):
    os.makedirs(work_dir, exist_ok=True)
    pdbs_dir = os.path.join(work_dir, "pdbs")
    if os.path.exists(pdbs_dir):
        shutil.rmtree(pdbs_dir)
    os.makedirs(pdbs_dir)
    shutil.copy(pdb_path, pdbs_dir)

    parsed = os.path.join(work_dir, "parsed.jsonl")
    assigned = os.path.join(work_dir, "assigned.jsonl")
    fixed = os.path.join(work_dir, "fixed.jsonl")

    run([sys.executable, f"{MPNN_HELPER}/parse_multiple_chains.py",
         f"--input_path={pdbs_dir}", f"--output_path={parsed}"])
    run([sys.executable, f"{MPNN_HELPER}/assign_fixed_chains.py",
         f"--input_path={parsed}", f"--output_path={assigned}", f"--chain_list={chain}"])
    run([sys.executable, f"{MPNN_HELPER}/make_fixed_positions_dict.py",
         f"--input_path={parsed}", f"--output_path={fixed}", f"--chain_list={chain}",
         f"--position_list={resnum}", "--specify_non_fixed"])
    run([sys.executable, MPNN_RUN,
         "--suppress_print", "1",
         "--jsonl_path", parsed, "--chain_id_jsonl", assigned, "--fixed_positions_jsonl", fixed,
         "--out_folder", work_dir, "--num_seq_per_target", "8", "--batch_size", "8",
         "--conditional_probs_only", "1"])

    tag = os.path.basename(pdb_path)[:-4]
    npz = np.load(os.path.join(work_dir, "conditional_probs_only", tag + ".npz"))
    log_p = npz["log_p"]
    if log_p.ndim == 3:
        log_p = log_p.mean(axis=0)
    mask = npz["design_mask"].astype(bool)
    log_p = log_p[mask]
    if log_p.shape[0] == 0:
        raise ValueError(f"scored position not found for {tag}")
    return log_p[0]  # 21-dim log-prob vector at the scored residue


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdbdir", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--resnum", type=int, default=185, help="PDB residue number (chain B is 1-189, no offset)")
    ap.add_argument("--chain", default="B")
    ap.add_argument("--on_target_aa", default="D", help="G12D Asp at p5")
    ap.add_argument("--work_dir", default="/tmp/mpnn_spec_work")
    args = ap.parse_args()

    aa_idx = {aa: i for i, aa in enumerate(ALPHABET)}
    on_idx = aa_idx[args.on_target_aa]
    ala_idx = aa_idx["A"]

    rows = []
    for fname in sorted(os.listdir(args.pdbdir)):
        if not fname.endswith(".pdb"):
            continue
        tag = fname[:-4]
        try:
            log_p = score_one(os.path.join(args.pdbdir, fname), args.resnum, args.chain, args.work_dir)
        except Exception as e:
            print(f"FAILED {tag}: {e}")
            continue
        log_p_on = float(log_p[on_idx])
        log_p_ala = float(log_p[ala_idx])
        mpnn_spec_score = log_p_on - log_p_ala
        rows.append((tag, log_p_on, log_p_ala, mpnn_spec_score))
        print(f"{tag}\tlog_p({args.on_target_aa})={log_p_on:.3f}\tlog_p(A)={log_p_ala:.3f}\tspec_score={mpnn_spec_score:.3f}")

    with open(args.out_csv, "w") as f:
        f.write("tag,log_p_on_target,log_p_ala,mpnn_spec_score\n")
        for tag, lp_on, lp_ala, spec in rows:
            f.write(f"{tag},{lp_on},{lp_ala},{spec}\n")
    print(f"wrote {len(rows)} rows to {args.out_csv}")


if __name__ == "__main__":
    main()
