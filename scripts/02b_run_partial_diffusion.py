#!/usr/bin/env python
"""Round-2 backbone generation via RFdiffusion partial diffusion, seeded from round-1
backbones with the best peptide contact (per docs/02_methods_and_inputs.md's "if too few
designs pass, loosen gates ... carry forward and improve via partial diffusion" guidance).

Each seed's binder chain (A) is noised/denoised at diffuser.partial_T (paper range 12-25 of
50 steps); the target chain (B) stays fixed. Unlike full de novo generation, the contig must
match the seed's *existing* binder length exactly (no sampled range).

Usage:
    python 02b_run_partial_diffusion.py --seeds_csv <csv with 'tag' column> \
        --backbone_dir <dir of seed pdbs> --out_prefix <prefix> \
        --variants_per_seed 10 --partial_T 18 --hotspot_res B184,B186,B187

If seeds_csv has a 'partial_T' column, it overrides --partial_T per-seed (e.g. to scale
partial_T by each seed's round-1 AF2 quality -- better seeds get a lighter touch).
"""
import argparse
import csv
import os
import subprocess


def chain_a_length(pdb_path):
    seen = set()
    for line in open(pdb_path):
        if line.startswith("ATOM") and line[21] == "A" and line[12:16].strip() == "CA":
            seen.add(line[22:27])
    return len(seen)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds_csv", required=True)
    ap.add_argument("--backbone_dir", required=True)
    ap.add_argument("--out_prefix", required=True)
    ap.add_argument("--variants_per_seed", type=int, default=10)
    ap.add_argument("--partial_T", type=int, default=18)
    ap.add_argument("--hotspot_res", default="B184,B186,B187")
    ap.add_argument("--n_seeds", type=int, default=None, help="limit to top N rows of seeds_csv")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)

    rows = list(csv.DictReader(open(args.seeds_csv)))
    if args.n_seeds:
        rows = rows[: args.n_seeds]

    for row in rows:
        tag = row["tag"]
        pdb_path = os.path.join(args.backbone_dir, tag + ".pdb")
        if not os.path.exists(pdb_path):
            print(f"SKIP {tag}: pdb not found at {pdb_path}")
            continue
        binder_len = chain_a_length(pdb_path)
        contig = f"[B1-189/0 {binder_len}-{binder_len}]"
        out_prefix = f"{args.out_prefix}_{tag}"
        partial_T = int(row["partial_T"]) if row.get("partial_T") else args.partial_T

        cmd = [
            "/venv/SE3nv/bin/python", "/workspace/RFdiffusion/scripts/run_inference.py",
            f"inference.output_prefix={out_prefix}",
            f"inference.input_pdb={pdb_path}",
            "inference.write_trajectory=False",
            f"contigmap.contigs={contig}",
            f"ppi.hotspot_res=[{args.hotspot_res}]",
            f"inference.num_designs={args.variants_per_seed}",
            f"diffuser.partial_T={partial_T}",
        ]
        print(f"=== seed {tag} (binder_len={binder_len}, partial_T={partial_T}) ===")
        subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
