#!/usr/bin/env python
"""Turn 05_contact_filter.py scores into a runlist of backbone tags that pass the cheap
peptide-contact pre-filter, for use with dl_interface_design.py / predict.py's -runlist.

Gate (round 1, set from the observed score distribution over 851 round-1 backbones --
see docs/03_design_log.md): peptide_cms > 10 and p5_cms > 1. This is a generous gate
(drops ~7% of backbones) since RFdiffusion's poly-Gly placeholder sequence understates
true contact area relative to the final ProteinMPNN-designed sequence -- it's meant to
kill clear non-contacts cheaply, not to be the primary selectivity filter (AF2 is).

Usage:
    python 06_select_contact_pass.py --csv contact_filter_backbones.csv \
        --min_peptide_cms 10 --min_p5_cms 1 \
        --out_runlist passing_backbones.txt [--mpnn_suffix _dldesign_0_cycle1]
"""
import argparse
import csv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--min_peptide_cms", type=float, default=10.0)
    ap.add_argument("--min_p5_cms", type=float, default=1.0)
    ap.add_argument("--out_runlist", required=True)
    ap.add_argument("--mpnn_suffix", default="", help="appended to each tag, e.g. _dldesign_0_cycle1")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.csv)))
    passing = [
        r["tag"]
        for r in rows
        if float(r["peptide_cms"]) > args.min_peptide_cms and float(r["p5_cms"]) > args.min_p5_cms
    ]
    with open(args.out_runlist, "w") as f:
        for tag in passing:
            f.write(tag + args.mpnn_suffix + "\n")

    print(f"{len(passing)}/{len(rows)} backbones pass (peptide_cms>{args.min_peptide_cms}, "
          f"p5_cms>{args.min_p5_cms}); wrote {args.out_runlist}")


if __name__ == "__main__":
    main()
