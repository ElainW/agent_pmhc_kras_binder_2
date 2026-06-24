import os
import sys

import numpy as np
import predict_utils

AA3to1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def read_chain_seq_and_lines(pdb_path, chain_id):
    seq = []
    seen = set()
    lines = []
    for line in open(pdb_path):
        if line.startswith("ATOM") and line[21] == chain_id:
            lines.append(line)
            if line[12:16].strip() == "CA":
                rn = line[22:27]
                if rn not in seen:
                    seen.add(rn)
                    seq.append(AA3to1.get(line[17:20].strip(), "X"))
    return "".join(seq), lines


def build_combined_template_pdb(binder_pdb, target_pdb, target_chains, out_path):
    """chain A = binder (from binder_pdb), chain B = concatenated target_chains (from
    target_pdb), renumbered continuously so residue numbers don't collide across the
    concatenated chains (e.g. MHC 1-275 then peptide 1-9 would otherwise both have '1'-'9',
    causing load_pdb_coords to silently drop the second chain's overlapping residues)."""
    _, binder_lines = read_chain_seq_and_lines(binder_pdb, "A")
    out_lines = list(binder_lines)
    counter = 0
    prev_key = None
    for ch in target_chains:
        _, lines = read_chain_seq_and_lines(target_pdb, ch)
        for line in lines:
            key = (ch, line[22:27])
            if key != prev_key:
                counter += 1
                prev_key = key
            newline = line[:21] + "B" + f"{counter:>4}" + line[26:]
            out_lines.append(newline)
    with open(out_path, "w") as f:
        f.writelines(out_lines)
        f.write("END\n")


def main():
    binder_pdb = sys.argv[1]
    target_pdb = sys.argv[2]  # e.g. 9UV8.pdb or 8I5E.pdb
    mhc_chain = sys.argv[3]   # e.g. 'A' for 9UV8, 'H' for 8I5E
    peptide_chain = sys.argv[4]  # e.g. 'C' for 9UV8, 'P' for 8I5E
    out_prefix = sys.argv[5]
    weights = "/workspace/designs/prep/pmhc_fold_weights/mixed_mhc_pae_run6_af_mhc_params_20640.pkl"

    binder_seq, _ = read_chain_seq_and_lines(binder_pdb, "A")
    mhc_seq, _ = read_chain_seq_and_lines(target_pdb, mhc_chain)
    peptide_seq, _ = read_chain_seq_and_lines(target_pdb, peptide_chain)

    query_sequence = binder_seq + mhc_seq + peptide_seq
    chainbreak_sequence = f"{binder_seq}/{mhc_seq}/{peptide_seq}"
    print("query:", chainbreak_sequence)

    combo_pdb = out_prefix + "_template.pdb"
    build_combined_template_pdb(binder_pdb, target_pdb, [mhc_chain, peptide_chain], combo_pdb)

    identity_alignment = {i: i for i in range(len(query_sequence))}
    tf = predict_utils.create_single_template_features(
        query_sequence, combo_pdb, identity_alignment, "template1",
        allow_chainbreaks=True, allow_skipped_lines=True,
    )
    template_features = predict_utils.compile_template_features([tf])

    model_runners = predict_utils.load_model_runners(
        model_names=["model_2_ptm_ft"],
        crop_size=len(query_sequence),
        data_dir="",
        model_params_files=[weights],
    )

    msa = [query_sequence]
    deletion_matrix = [[0] * len(query_sequence)]

    metrics = predict_utils.run_alphafold_prediction(
        query_sequence=query_sequence,
        msa=msa,
        deletion_matrix=deletion_matrix,
        chainbreak_sequence=chainbreak_sequence,
        template_features=template_features,
        model_runners=model_runners,
        out_prefix=out_prefix,
    )
    print(metrics)
    np.save(out_prefix + "_lens.npy", np.array([len(binder_seq), len(mhc_seq), len(peptide_seq)]))


if __name__ == "__main__":
    main()
