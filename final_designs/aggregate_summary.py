#!/usr/bin/env python3
"""Aggregate all per-candidate stats for the 6 ipsae>0 final designs into one summary TSV."""

import csv
import sys

STATS_DIR = "/workspace/final_designs/ipsae_gt0_stats"
CANDIDATES_CSV = "/workspace/final_designs/candidates_with_sequences.csv"
OUT_TSV = "/workspace/final_designs/top_candidate_summary_stats.tsv"


def read_tsv(path):
    rows = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            key = row["design"]
            rows[key] = row
    return rows


def read_candidates(path):
    rows = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows[row["candidate_id"]] = row
    return rows


def main():
    cand = read_candidates(CANDIDATES_CSV)
    af3 = read_tsv(f"{STATS_DIR}/af3_design_stats.tsv")
    relax = read_tsv(f"{STATS_DIR}/fastrelax_af3_scores.tsv")
    audit = read_tsv(f"{STATS_DIR}/redesign_audit.tsv")

    # The 6 designs in ipsae_gt0_stats
    designs = list(af3.keys())

    out_cols = [
        # Identity
        "design",
        "binder_sequence",
        "round",
        "seed_lineage",
        "binder_len",
        # AF3 ipTM + ipSAE (on-target)
        "af3_iptm",
        "af3_ptm",
        "af3_ranking_score",
        "af3_fraction_disordered",
        "ipsae_binder_peptide",          # current run (binder-pep chain pair)
        "ipsae_bp",                       # same as above (bp direction)
        "ipsae_pb",                       # reverse direction
        "ipsae_n_contacts",
        "ipsae_d0",
        "mean_pae_bp",
        "mean_pae_contact",
        # ipSAE outward-asp conformer re-run
        "ipsae_bp_outward_run",
        "cms_peptide_total_outward_run",
        # Off-target from candidates.csv
        "af3_off_iptm",
        "af3_off_binder_pep_iptm",
        "ipsae_af3_binder_pep_candidates",  # from candidates CSV
        "their_ipsae_binder_pep_on",
        "their_ipsae_binder_pep_off",
        "passes_af3_gates",
        # Secondary structure
        "ss_helix_pct",
        "ss_sheet_pct",
        "ss_loop_pct",
        "n_helix",
        "n_sheet",
        "n_loop",
        "dssp_string",
        # CMS per-peptide-position (Asp inward conformation)
        "asp_p5_chi1",
        "asp_p5_conformation",
        "cms_p1",
        "cms_p2",
        "cms_p3",
        "cms_p4",
        "cms_p5",                         # CMS for p5 (whichever conformation is observed)
        "cms_p5_inward",                  # CMS when Asp chi1 is inward
        "cms_p5_outward",                 # CMS when Asp chi1 is outward
        "cms_p6",
        "cms_p7",
        "cms_p8",
        "cms_p9",
        "cms_peptide_total",
        # Contacts — overall
        "n_all_contacts",
        "n_polar_contacts",
        "n_hydrophobic_contacts",
        "n_hbonds",
        "hbond_available",
        "binder_res_in_contact",
        # Hotspot (all peptide positions) contacts
        "n_contacts_hotspot_total_all",
        "n_contacts_hotspot_total_polar",
        "n_contacts_hotspot_total_hydro",
        "n_contacts_hotspot_total_hbond",
        # p5 contacts — all conformations
        "n_contacts_hotspot_p5_all",
        "n_contacts_hotspot_p5_polar",
        "n_contacts_hotspot_p5_hydro",
        "n_contacts_hotspot_p5_hbond",
        "n_saltbridge_hotspot_p5",
        "saltbridge_hotspot_p5_binder_resnums",
        "saltbridge_hotspot_p5_binder_aas",
        # p5 contacts — inward/outward separated
        "n_contacts_p5_all_inward",
        "n_contacts_p5_all_outward",
        "n_contacts_p5_polar_inward",
        "n_contacts_p5_polar_outward",
        "n_contacts_p5_hbond_inward",
        "n_contacts_p5_hbond_outward",
        "contacts_polar_p5_inward",       # which binder residues make polar contacts (inward)
        "contacts_polar_p5_outward",      # which binder residues make polar contacts (outward)
        # per-position contacts (all + hbond) for each peptide position
        "contacts_all_p1",  "contacts_all_p2",  "contacts_all_p3",  "contacts_all_p4",
        "contacts_all_p5",  "contacts_all_p6",  "contacts_all_p7",  "contacts_all_p8",  "contacts_all_p9",
        "contacts_hbond_p1","contacts_hbond_p2","contacts_hbond_p3","contacts_hbond_p4",
        "contacts_hbond_p5","contacts_hbond_p6","contacts_hbond_p7","contacts_hbond_p8","contacts_hbond_p9",
        "contacts_polar_p1","contacts_polar_p2","contacts_polar_p3","contacts_polar_p4",
        "contacts_polar_p5","contacts_polar_p6","contacts_polar_p7","contacts_polar_p8","contacts_polar_p9",
        "contacts_hydro_p1","contacts_hydro_p2","contacts_hydro_p3","contacts_hydro_p4",
        "contacts_hydro_p5","contacts_hydro_p6","contacts_hydro_p7","contacts_hydro_p8","contacts_hydro_p9",
        # Rosetta interface scores (from fastrelax)
        "rosetta_score_before",
        "rosetta_score_after",
        "rosetta_score_delta",
        "dG_separated",
        "dG_separated/dSASAx100",
        "dSASA_int",
        "dSASA_polar",
        "dSASA_hphobic",
        "sc_value",
        "packstat",
        "per_residue_energy_int",
        "hbonds_int",
        "hbond_E_fraction",
        "nres_int",
        "nres_all",
        "complex_normalized",
        # Rosetta detail energy terms (from fastrelax)
        "fa_atr", "fa_rep", "fa_sol",
        "fa_elec", "hbond_sr_bb", "hbond_lr_bb", "hbond_bb_sc", "hbond_sc",
        "fa_dun", "p_aa_pp", "ref", "rama_prepro",
        # BindCraft filters
        "buns_delta_unsat",
        "delta_unsatHbonds",
        "buns_pass",
        "surface_hydrophobicity",
        "interface_n_K",
        "interface_n_M",
        "pass_all_bindcraft",
        # Charge / Cys audit
        "net_charge",
        "n_cys",
        "n_met",
        "flag_charge_redesign",
        "flag_cys_redesign",
        "pass_charge",
        "pass_cys",
        "n_interface_fixed",
        "interface_fixed_resnums",
        "pass_all",
        # From candidates.csv
        "af2_pae_interaction",
        "pmhcfold_delta_on_minus_off",
        "mpnn_spec_score",
        "real_cms_p5",
        "notes",
    ]

    rows_out = []
    for d in designs:
        a = af3[d]
        r = relax.get(d, {})
        au = audit.get(d, {})
        c = cand.get(d, {})

        row = {}
        for col in out_cols:
            if col == "binder_sequence":
                row[col] = c.get("binder_sequence", "")
            elif col == "round":
                row[col] = c.get("round", "")
            elif col == "seed_lineage":
                row[col] = c.get("seed_lineage", "")
            elif col == "binder_len":
                row[col] = c.get("binder_len", "")
            elif col == "af3_off_iptm":
                row[col] = c.get("af3_off_iptm", "")
            elif col == "af3_off_binder_pep_iptm":
                row[col] = c.get("af3_off_binder_pep_iptm", "")
            elif col == "ipsae_af3_binder_pep_candidates":
                row[col] = c.get("ipsae_af3_binder_pep", "")
            elif col == "their_ipsae_binder_pep_on":
                row[col] = c.get("their_ipsae_binder_pep_on", "")
            elif col == "their_ipsae_binder_pep_off":
                row[col] = c.get("their_ipsae_binder_pep_off", "")
            elif col == "passes_af3_gates":
                row[col] = c.get("passes_af3_gates", "")
            elif col == "af2_pae_interaction":
                row[col] = c.get("af2_pae_interaction", "")
            elif col == "pmhcfold_delta_on_minus_off":
                row[col] = c.get("pmhcfold_delta_on_minus_off", "")
            elif col == "mpnn_spec_score":
                row[col] = c.get("mpnn_spec_score", "")
            elif col == "real_cms_p5":
                row[col] = c.get("real_cms_p5", "")
            elif col == "notes":
                row[col] = c.get("notes", "")
            elif col in r:
                row[col] = r[col]
            elif col in a:
                row[col] = a[col]
            elif col in au:
                row[col] = au[col]
            else:
                row[col] = ""
        rows_out.append(row)

    with open(OUT_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_cols, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"Wrote {len(rows_out)} rows x {len(out_cols)} columns to {OUT_TSV}")


if __name__ == "__main__":
    main()
