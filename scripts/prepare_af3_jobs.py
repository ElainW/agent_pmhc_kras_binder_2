#!/usr/bin/env python3
"""Build AlphaFold Server (AF3) batch-upload JSON files for the specificity screen
and final confirmation steps of the KRAS(G12D) pMHC binder campaign.

For each candidate binder sequence, emits two jobs:
  - "<id>_on"  : binder + 9UV8 on-target complex (HLA-A*11:01 + KRAS G12D peptide VVGADGVGK + b2m)
  - "<id>_off" : binder + 8I5E off-target complex (same HLA + WT KRAS peptide VVGAGGVGK + b2m)
HLA heavy chain and b2m sequences are held identical between on/off jobs (docs/02_methods_and_inputs.md
"Sequences (9UV8 Construct)") so the peptide identity is the only variable -- isolating the
G12D-vs-WT specificity test.

The binder chain is folded with no MSA (single-sequence unpairedMsa/pairedMsa, no templates) per
the validated-design convention observed in paper_output/af3_nomsa/*/*_data.json. Target chains are
left MSA-free for input (AF3 Server computes MSA/templates server-side at submission).

Jobs are grouped into batches of <= --batch_size (default 30, matching the AF3 Server daily quota)
for manual upload.

Usage:
    python prepare_af3_jobs.py --candidates candidates.csv --out_dir designs/round1/af3_jobs
candidates.csv must have columns: candidate_id,binder_sequence
"""
import argparse
import csv
import json
import os

HLA_A1101_HEAVY = (
    "MGSHSMRYFYTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQTDRVDLGTLRGYYNQSEDGSHTIQI"
    "MYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAAHAAEQQRAYLEGRCVEWLRRYLENGKETLQRTDPPKTHMTHHPIS"
    "DHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWE"
)
B2M = (
    "MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM"
)
PEPTIDE_ON = "VVGADGVGK"   # 9UV8 KRAS G12D (Asp at p5)
PEPTIDE_OFF = "VVGAGGVGK"  # 8I5E KRAS WT (Gly at p5)


def no_msa_protein(chain_id, sequence):
    return {
        "protein": {
            "id": chain_id,
            "sequence": sequence,
            "modifications": [],
            "unpairedMsa": f">Original query\n{sequence}\n",
            "pairedMsa": f">Original query\n{sequence}\n",
            "templates": [],
        }
    }


def target_protein(chain_id, sequence):
    # No MSA/templates supplied on input -- AF3 Server fetches its own at submission time.
    return {"protein": {"id": chain_id, "sequence": sequence}}


def build_job(name, binder_seq, peptide_seq, seeds):
    return {
        "dialect": "alphafold3",
        "version": 2,
        "name": name,
        "modelSeeds": seeds,
        "sequences": [
            no_msa_protein("A", binder_seq),
            target_protein("B", HLA_A1101_HEAVY),
            target_protein("C", peptide_seq),
            target_protein("D", B2M),
        ],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True, help="CSV with columns candidate_id,binder_sequence")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--batch_size", type=int, default=30, help="AF3 Server daily job quota")
    ap.add_argument("--seeds", default="1,2,3,4,5", help="comma-separated modelSeeds per job")
    args = ap.parse_args()

    seeds = [int(s) for s in args.seeds.split(",")]
    os.makedirs(args.out_dir, exist_ok=True)

    jobs = []
    with open(args.candidates) as f:
        for row in csv.DictReader(f):
            cid = row["candidate_id"].strip()
            seq = row["binder_sequence"].strip()
            jobs.append(build_job(f"{cid}_on", seq, PEPTIDE_ON, seeds))
            jobs.append(build_job(f"{cid}_off", seq, PEPTIDE_OFF, seeds))

    manifest = []
    for i in range(0, len(jobs), args.batch_size):
        batch = jobs[i : i + args.batch_size]
        batch_path = os.path.join(args.out_dir, f"af3_batch_{i // args.batch_size + 1:03d}.json")
        with open(batch_path, "w") as f:
            json.dump(batch, f, indent=1)
        manifest.append((batch_path, [j["name"] for j in batch]))

    print(f"Wrote {len(jobs)} jobs ({len(jobs)//2} candidates) across {len(manifest)} batch file(s):")
    for path, names in manifest:
        print(f"  {path}  ({len(names)} jobs: {', '.join(names)})")


if __name__ == "__main__":
    main()
