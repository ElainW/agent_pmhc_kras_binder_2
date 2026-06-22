#!/usr/bin/env python
"""Fast monomer-foldability triage for ProteinMPNN-designed binder sequences, using
ESMFold (facebook/esmfold_v1 via HuggingFace transformers) instead of AF2 monomer --
no MSA, single forward pass, ~10s vs AF2's ~180-190s per sequence.

Motivation: matching the paper's 4-32 MPNN sequences/backbone (see docs/03_design_log.md)
multiplies the number of structures that would otherwise need the expensive AF2-complex
pae_interaction gate. ESMFold screens binder-alone foldability for *all* sequences per
backbone cheaply, and only the single best-folding sequence per backbone (by mean pLDDT)
is forwarded to AF2-complex -- so AF2 call volume stays flat at ~1/backbone while MPNN
still samples the full range.

Expects input PDBs named "<backbone_tag>_dldesign_<i>.pdb" (dl_interface_design.py's
-relax_cycles 0 naming) with the binder as chain A.

Usage:
    python 07_esmfold_filter.py --pdbdir <mpnn_out_dir> --out_csv esmfold_scores.csv \
        --out_runlist best_per_backbone.txt [--mpnn_suffix _dldesign]
"""
import argparse
import csv
import glob
import os
import re
import time

import torch
from transformers import AutoTokenizer
from transformers.models.esm.modeling_esmfold import EsmForProteinFolding

AA3to1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def binder_chain_a_sequence(pdb_path):
    seq = []
    seen = set()
    for line in open(pdb_path):
        if line.startswith("ATOM") and line[21] == "A" and line[12:16].strip() == "CA":
            resnum = line[22:27]
            if resnum in seen:
                continue
            seen.add(resnum)
            seq.append(AA3to1.get(line[17:20].strip(), "X"))
    return "".join(seq)


def backbone_tag(mpnn_tag, mpnn_suffix):
    return re.sub(re.escape(mpnn_suffix) + r"_\d+$", "", mpnn_tag)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdbdir", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--out_runlist", required=True)
    ap.add_argument("--mpnn_suffix", default="_dldesign")
    ap.add_argument("--min_plddt", type=float, default=70.0, help="below this, drop the backbone entirely")
    args = ap.parse_args()

    tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1")
    model = model.cuda().eval()

    rows = []
    pdbs = sorted(glob.glob(os.path.join(args.pdbdir, "*.pdb")))
    for pdb_path in pdbs:
        tag = os.path.basename(pdb_path)[:-4]
        seq = binder_chain_a_sequence(pdb_path)
        if not seq:
            print(f"SKIP {tag}: no chain A found")
            continue
        inputs = tok([seq], return_tensors="pt", add_special_tokens=False)
        inputs = {k: v.cuda() for k, v in inputs.items()}
        t0 = time.time()
        with torch.no_grad():
            out = model(**inputs)
        plddt = out["plddt"][0, ..., 1].mean().item() if out["plddt"].dim() > 2 else out["plddt"].mean().item()
        dt = time.time() - t0
        rows.append({"tag": tag, "backbone": backbone_tag(tag, args.mpnn_suffix), "plddt": plddt, "len": len(seq)})
        print(f"{tag}\tplddt={plddt:.1f}\tlen={len(seq)}\t{dt:.1f}s")

    with open(args.out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["tag", "backbone", "plddt", "len"])
        writer.writeheader()
        writer.writerows(rows)

    best_per_backbone = {}
    for r in rows:
        b = r["backbone"]
        if b not in best_per_backbone or r["plddt"] > best_per_backbone[b]["plddt"]:
            best_per_backbone[b] = r

    with open(args.out_runlist, "w") as f:
        for b, r in best_per_backbone.items():
            if r["plddt"] >= args.min_plddt:
                f.write(r["tag"] + "\n")

    n_pass = sum(1 for r in best_per_backbone.values() if r["plddt"] >= args.min_plddt)
    print(f"{len(rows)} sequences across {len(best_per_backbone)} backbones; "
          f"{n_pass}/{len(best_per_backbone)} backbones have a best-sequence plddt>={args.min_plddt}; "
          f"wrote {args.out_runlist}")


if __name__ == "__main__":
    main()
