#!/usr/bin/env python
"""Compare binder-to-p5 (mutant Asp) min heavy-atom distance between the original
pose-templated design and an independently-folded structure (e.g. AF3 with no
binder templating). AF2-initial-guess and pmhc_fold.py are both pose-templated --
they measure whether the *intended* pose is stable, not whether the model
independently discovers the contact. A contact that's marginal in the original
design and grows much longer in an untemplated refold is a sign the binding/
specificity signal from the templated tools may be a pose-bias artifact rather
than a real, robust interaction (see docs/03_design_log.md 2026-06-25 entry).

Usage:
    python 12_check_p5_contact.py <pdb_or_cif> [--two_chain]
--two_chain: input is the raw RFdiffusion/MPNN design (chain A=binder, chain B=
             MHC+peptide combined, p5 = chain B residue 185). Default: AF3-style
             4-chain output (A=binder, B=MHC, C=peptide, D=b2m; p5 = chain C res 5).
"""
import sys
import argparse
import pyrosetta
from pyrosetta import pose_from_file

pyrosetta.init("-mute all -detect_disulf false -ignore_unrecognized_res true")


def find_p5(pose, two_chain):
    if two_chain:
        for r in range(1, pose.size() + 1):
            if pose.pdb_info().chain(r) == "B" and pose.pdb_info().number(r) == 185:
                return r
    else:
        pep_res = [r for r in range(1, pose.size() + 1) if pose.pdb_info().chain(r) == "C"]
        return pep_res[4]
    return None


def min_heavy_atom_dist(pose, binder_res, p5):
    p5res = pose.residue(p5)
    min_dist, min_pair = 999.0, None
    for r in binder_res:
        res = pose.residue(r)
        for ai in range(1, res.natoms() + 1):
            if res.atom_type(ai).element() == "H":
                continue
            xa = res.xyz(ai)
            for bi in range(1, p5res.natoms() + 1):
                if p5res.atom_type(bi).element() == "H":
                    continue
                d = (xa - p5res.xyz(bi)).norm()
                if d < min_dist:
                    min_dist = d
                    min_pair = (pose.pdb_info().number(r), res.name(), res.atom_name(ai), p5res.atom_name(bi))
    return min_dist, min_pair


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("structure")
    ap.add_argument("--two_chain", action="store_true")
    args = ap.parse_args()

    pose = pose_from_file(args.structure)
    p5 = find_p5(pose, args.two_chain)
    binder_res = [r for r in range(1, pose.size() + 1) if pose.pdb_info().chain(r) == "A"]
    min_dist, min_pair = min_heavy_atom_dist(pose, binder_res, p5)
    print(f"{args.structure}: p5={pose.residue(p5).name()} min_dist={min_dist:.2f} at {min_pair}")


if __name__ == "__main__":
    main()
