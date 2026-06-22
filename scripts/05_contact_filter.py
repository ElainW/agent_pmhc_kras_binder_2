#!/usr/bin/env python
"""Cheap geometric pre-filter for RFdiffusion backbones (or MPNN-designed PDBs):
contact molecular surface between the binder (chain A) and the peptide (chain B,
residues 181-189), replicating paper/pMHCI_binder_design/software/contact/just_contact_patch.xml
directly in PyRosetta (no compiled rosetta_scripts binary needed).

Reports both the whole-peptide contact and the Asp(p5)=185-only contact so low-cost
backbones that miss the mutant residue specifically can be deprioritized even if they
touch other peptide positions.

Usage:
    python 05_contact_filter.py --pdbdir <dir> --out_csv <path> [--glob '*.pdb']
"""
import argparse
import csv
import glob
import os

import pyrosetta
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector, ResidueIndexSelector
from pyrosetta.rosetta.protocols.simple_filters import ContactMolecularSurfaceFilter

PEPTIDE_RANGE = ",".join(f"{i}B" for i in range(181, 190))  # PDB-numbered, chain B
P5_RESIDUE = "185B"


def make_filter(target_sel, binder_sel):
    f = ContactMolecularSurfaceFilter()
    f.selector1(binder_sel)
    f.selector2(target_sel)
    f.distance_weight(0.5)
    return f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdbdir", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--glob", default="*.pdb")
    args = ap.parse_args()

    pyrosetta.init("-mute all -beta_nov16")

    binder_sel = ChainSelector("A")
    peptide_sel = ResidueIndexSelector(PEPTIDE_RANGE)
    p5_sel = ResidueIndexSelector(P5_RESIDUE)
    peptide_filter = make_filter(peptide_sel, binder_sel)
    p5_filter = make_filter(p5_sel, binder_sel)

    rows = []
    pdbs = sorted(glob.glob(os.path.join(args.pdbdir, args.glob)))
    for pdb_path in pdbs:
        tag = os.path.basename(pdb_path)[:-4]
        try:
            pose = pyrosetta.pose_from_pdb(pdb_path)
            peptide_cms = peptide_filter.compute(pose)
            p5_cms = p5_filter.compute(pose)
        except Exception as e:
            print(f"FAILED {tag}: {e}")
            continue
        rows.append({"tag": tag, "peptide_cms": peptide_cms, "p5_cms": p5_cms})
        print(f"{tag}\tpeptide_cms={peptide_cms:.2f}\tp5_cms={p5_cms:.2f}")

    with open(args.out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["tag", "peptide_cms", "p5_cms"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {args.out_csv}")


if __name__ == "__main__":
    main()
