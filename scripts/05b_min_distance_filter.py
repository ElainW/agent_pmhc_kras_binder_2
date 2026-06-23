#!/usr/bin/env python
"""Cheap geometric pre-filter for raw RFdiffusion backbones: minimum distance from the
binder (chain A, poly-Gly placeholder -- no real side chains yet) to the peptide p5
residue (chain B, residue 185 -- real Asp side chain, since RFdiffusion never touches
the target).

Replaces 05_contact_filter.py's PyRosetta ContactMolecularSurface for the pre-MPNN gate:
CMS needs real side-chain atoms to mean anything (surface area), which the binder doesn't
have yet at this stage -- a poly-Gly structure understates true contact in a way that's
hard to interpret. A raw minimum heavy-atom distance has no such dependency and is a more
honest measure of "does this backbone's trace even reach the peptide" at the backbone-only
stage. (CMS/Rosetta scoring remains appropriate *after* ProteinMPNN + AF2, once there's a
real designed sequence and a real predicted structure to compute it on.)

RFdiffusion's raw output PDB has only backbone atoms (N, CA, C, O) throughout -- including
the "fixed" target -- so a side-chain-based distance isn't available at this stage. Reports:
  - min_ca_dist: minimum CA(chain A)-CA(p5) distance -- backbone-trace proxy
  - min_bb_dist: minimum distance over all backbone atoms (N/CA/C/O) on both sides --
    slightly tighter/more honest "closest approach" than CA-only

Usage:
    python 05b_min_distance_filter.py --pdbdir <dir> --out_csv <path> [--glob '*.pdb']
"""
import argparse
import csv
import glob
import os

P5_RESNUM = 185
BACKBONE_ATOMS = {"N", "CA", "C", "O"}


def parse_atoms(pdb_path):
    chain_a_atoms = []
    p5_atoms = []
    for line in open(pdb_path):
        if not line.startswith("ATOM"):
            continue
        chain = line[21]
        atom_name = line[12:16].strip()
        if atom_name not in BACKBONE_ATOMS:
            continue
        resnum = int(line[22:26])
        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        if chain == "A":
            chain_a_atoms.append((atom_name, (x, y, z)))
        elif chain == "B" and resnum == P5_RESNUM:
            p5_atoms.append((atom_name, (x, y, z)))
    return chain_a_atoms, p5_atoms


def dist(p1, p2):
    return ((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2 + (p1[2] - p2[2]) ** 2) ** 0.5


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdbdir", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--glob", default="*.pdb")
    args = ap.parse_args()

    rows = []
    for pdb_path in sorted(glob.glob(os.path.join(args.pdbdir, args.glob))):
        tag = os.path.basename(pdb_path)[:-4]
        chain_a_atoms, p5_atoms = parse_atoms(pdb_path)
        if not chain_a_atoms or not p5_atoms:
            print(f"SKIP {tag}: missing chain A or p5 atoms")
            continue
        min_ca = min(
            dist(a[1], b[1]) for a in chain_a_atoms if a[0] == "CA" for b in p5_atoms if b[0] == "CA"
        )
        min_bb = min(dist(a[1], b[1]) for a in chain_a_atoms for b in p5_atoms)
        rows.append({"tag": tag, "min_ca_dist": min_ca, "min_bb_dist": min_bb})
        print(f"{tag}\tmin_ca_dist={min_ca:.2f}\tmin_bb_dist={min_bb:.2f}")

    with open(args.out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["tag", "min_ca_dist", "min_bb_dist"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {args.out_csv}")


if __name__ == "__main__":
    main()
