"""
16_compile_ipsae_gt0_stats.py

Compile complete stats for the 6 ipSAE > 0 designs from candidates_with_sequences.csv.
Merges af3_stats_13 (negdelta AF3 run) with af3_outward_asp_stats (outward-Asp AF3 run
for dldesign_6 and _7) and labels CMS p5 and polar contacts at p5 as "inward" or
"outward" based on the chi1 dihedral of the mutant Asp (p5) in each relaxed structure.

Chi1 classification: |chi1| < 120° → inward (gauche), else → outward (trans).

Outputs to /workspace/final_designs/ipsae_gt0_stats/:
  af3_design_stats.tsv       — full per-design stats with inward/outward p5 columns
  fastrelax_af3_scores.tsv   — Rosetta interface metrics
  redesign_audit.tsv         — charge / Cys redesign audit

Usage:
    /venv/dl_binder_design/bin/python3 16_compile_ipsae_gt0_stats.py
"""

import os
import sys
import math
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.vectors import calc_dihedral

# ── Paths ─────────────────────────────────────────────────────────────────────

WORKSPACE   = "/workspace"
STATS_13    = os.path.join(WORKSPACE, "designs", "af3_stats_13")
STATS_OUT   = os.path.join(WORKSPACE, "designs", "af3_outward_asp_stats")
BMI_DIR     = os.path.join(WORKSPACE, "designs", "binder_metrics_input")   # relaxed PDBs
OUT_DIR     = os.path.join(WORKSPACE, "final_designs", "ipsae_gt0_stats")

# ── The 6 ipSAE > 0 designs (ipsae_af3_binder_pep > 0.1) ────────────────────
# Maps canonical design name → name in af3_stats_13 (has _on suffix and
# _control_ infix for r3_r1b_273_28)

DESIGNS = {
    "r1b_403_dldesign_2":            "r1b_403_dldesign_2_on",
    "r3_r1b_273_28_dldesign_0":      "r3_r1b_273_28_dldesign_0_control_on",
    "r3_r1b_403_65_dldesign_4":      "r3_r1b_403_65_dldesign_4_on",
    "r3_r1b_870_87_dldesign_2":      "r3_r1b_870_87_dldesign_2_on",
    "r3_r1b_870_87_dldesign_6":      "r3_r1b_870_87_dldesign_6_on",
    "r3_r1b_870_87_dldesign_7":      "r3_r1b_870_87_dldesign_7_on",
}

# Designs that also have an "outward Asp" AF3 run in af3_outward_asp_stats
OUTWARD_DESIGNS = {
    "r3_r1b_870_87_dldesign_6",
    "r3_r1b_870_87_dldesign_7",
}

# ── Relaxed PDB path ──────────────────────────────────────────────────────────

def relaxed_pdb_path(canonical_name):
    """
    Return path to the relaxed PDB in binder_metrics_input.
    For r3_r1b_273_28_dldesign_0, the file is named with _control_.
    """
    if canonical_name == "r3_r1b_273_28_dldesign_0":
        fname = "r3_r1b_273_28_dldesign_0_control_relaxed.pdb"
    else:
        fname = f"{canonical_name}_relaxed.pdb"
    return os.path.join(BMI_DIR, fname)


# ── Asp p5 chi1 calculation ───────────────────────────────────────────────────

def get_p5_asp_chi1(pdb_path, peptide_chain="C"):
    """
    Return (resnum, chi1_degrees) for the Asp at position 5 in the peptide chain.
    Chi1 = N-CA-CB-CG dihedral (standard definition).
    Returns (None, float('nan')) on failure.
    """
    parser = PDBParser(QUIET=True)
    try:
        struct = parser.get_structure("s", pdb_path)
    except Exception as e:
        print(f"  ERROR loading {pdb_path}: {e}")
        return None, float("nan")

    for model in struct:
        for chain in model:
            if chain.id != peptide_chain:
                continue
            std_res = [r for r in chain if r.id[0] == " "]
            if len(std_res) < 5:
                continue
            r = std_res[4]           # 0-indexed position 4 = p5
            if r.resname != "ASP":
                print(f"  WARNING: p5 residue is {r.resname} (expected ASP) in {pdb_path}")
            try:
                chi1 = float(np.degrees(calc_dihedral(
                    r["N"].get_vector(),
                    r["CA"].get_vector(),
                    r["CB"].get_vector(),
                    r["CG"].get_vector(),
                )))
                return r.id[1], chi1
            except KeyError as e:
                print(f"  ERROR computing chi1 in {pdb_path}: missing atom {e}")
                return r.id[1], float("nan")
        break
    return None, float("nan")


def classify_chi1(chi1_deg):
    """'inward' if |chi1| < 120 (gauche); 'outward' if |chi1| >= 120 (trans)."""
    if math.isnan(chi1_deg):
        return "unknown"
    return "inward" if abs(chi1_deg) < 120.0 else "outward"


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # ── 1. Load af3_stats_13 data ───────────────────────────────────────────
    stats13_path   = os.path.join(STATS_13, "af3_design_stats.tsv")
    fr13_path      = os.path.join(STATS_13, "fastrelax_af3_scores.tsv")
    audit13_path   = os.path.join(STATS_13, "redesign_audit.tsv")

    df13    = pd.read_csv(stats13_path,  sep="\t")
    df_fr13 = pd.read_csv(fr13_path,     sep="\t")
    df_a13  = pd.read_csv(audit13_path,  sep="\t")

    # ── 2. Load af3_outward_asp_stats data (for dldesign_6 and _7 outward) ──
    stats_out_path  = os.path.join(STATS_OUT, "af3_design_stats.tsv")
    fr_out_path     = os.path.join(STATS_OUT, "fastrelax_af3_scores.tsv")
    df_out    = pd.read_csv(stats_out_path, sep="\t")
    df_fr_out = pd.read_csv(fr_out_path,    sep="\t")

    # ── 3. Build the 6-design stats table ───────────────────────────────────

    stats_rows   = []
    fr_rows      = []
    audit_rows   = []

    for canonical, stats13_name in DESIGNS.items():
        print(f"\n[{canonical}]")

        # --- a. Main stats from af3_stats_13 ---
        row13 = df13[df13["design"] == stats13_name]
        if row13.empty:
            print(f"  ERROR: {stats13_name} not found in af3_stats_13")
            continue
        rec = row13.iloc[0].to_dict()

        # Normalise design name to canonical (without _on / _control_)
        rec["design"] = canonical

        # Update relaxed_pdb to point at binder_metrics_input
        bmi_pdb = relaxed_pdb_path(canonical)
        if os.path.exists(bmi_pdb):
            rec["relaxed_pdb"] = bmi_pdb

        # --- b. Asp p5 chi1 and conformation ---
        _, chi1 = get_p5_asp_chi1(bmi_pdb)
        conformation = classify_chi1(chi1)
        rec["asp_p5_chi1"]       = round(chi1, 2) if not math.isnan(chi1) else float("nan")
        rec["asp_p5_conformation"] = conformation
        print(f"  chi1 = {chi1:.1f}°  → {conformation}")

        # --- c. Inward / outward CMS p5 and polar contacts (from af3_stats_13) ---
        cms_p5   = rec.get("cms_hotspot_p5",          float("nan"))
        polar_p5 = rec.get("n_contacts_hotspot_p5_polar", float("nan"))
        all_p5   = rec.get("n_contacts_hotspot_p5_all",   float("nan"))
        hbond_p5 = rec.get("n_contacts_hotspot_p5_hbond", float("nan"))
        cpolar_p5 = rec.get("contacts_polar_p5",           float("nan"))

        if conformation == "inward":
            rec["cms_p5_inward"]              = cms_p5
            rec["cms_p5_outward"]             = float("nan")
            rec["n_contacts_p5_all_inward"]   = all_p5
            rec["n_contacts_p5_all_outward"]  = float("nan")
            rec["n_contacts_p5_polar_inward"] = polar_p5
            rec["n_contacts_p5_polar_outward"]= float("nan")
            rec["n_contacts_p5_hbond_inward"] = hbond_p5
            rec["n_contacts_p5_hbond_outward"]= float("nan")
            rec["contacts_polar_p5_inward"]   = cpolar_p5
            rec["contacts_polar_p5_outward"]  = float("nan")
        else:
            rec["cms_p5_inward"]              = float("nan")
            rec["cms_p5_outward"]             = cms_p5
            rec["n_contacts_p5_all_inward"]   = float("nan")
            rec["n_contacts_p5_all_outward"]  = all_p5
            rec["n_contacts_p5_polar_inward"] = float("nan")
            rec["n_contacts_p5_polar_outward"]= polar_p5
            rec["n_contacts_p5_hbond_inward"] = float("nan")
            rec["n_contacts_p5_hbond_outward"]= hbond_p5
            rec["contacts_polar_p5_inward"]   = float("nan")
            rec["contacts_polar_p5_outward"]  = cpolar_p5

        # --- d. Outward conformation stats (dldesign_6 and _7 only) ---
        if canonical in OUTWARD_DESIGNS:
            out_name = f"{canonical}_on"
            row_out = df_out[df_out["design"] == out_name]
            if not row_out.empty:
                out = row_out.iloc[0]
                out_pdb = os.path.join(STATS_OUT, "relaxed_pdbs", f"{out_name}_relaxed.pdb")
                _, chi1_out = get_p5_asp_chi1(out_pdb)
                conf_out = classify_chi1(chi1_out)
                rec["asp_p5_chi1_outward_run"]        = round(chi1_out, 2)
                rec["asp_p5_conformation_outward_run"] = conf_out
                rec["cms_p5_outward"]             = out.get("cms_hotspot_p5",               float("nan"))
                rec["n_contacts_p5_all_outward"]  = out.get("n_contacts_hotspot_p5_all",    float("nan"))
                rec["n_contacts_p5_polar_outward"]= out.get("n_contacts_hotspot_p5_polar",  float("nan"))
                rec["n_contacts_p5_hbond_outward"]= out.get("n_contacts_hotspot_p5_hbond",  float("nan"))
                rec["contacts_polar_p5_outward"]  = out.get("contacts_polar_p5",             float("nan"))
                rec["ipsae_bp_outward_run"]        = out.get("ipsae_bp",                     float("nan"))
                rec["cms_peptide_total_outward_run"]= out.get("cms_peptide_total",            float("nan"))
                print(f"  outward run: chi1={chi1_out:.1f}°  cms_p5={rec['cms_p5_outward']:.2f}"
                      f"  polar_p5={rec['n_contacts_p5_polar_outward']}")
            else:
                print(f"  WARNING: {out_name} not found in af3_outward_asp_stats")

        print(f"  cms_p5_inward={rec['cms_p5_inward']}  cms_p5_outward={rec['cms_p5_outward']}")
        print(f"  polar_p5_inward={rec['n_contacts_p5_polar_inward']}  "
              f"polar_p5_outward={rec['n_contacts_p5_polar_outward']}")
        stats_rows.append(rec)

        # --- e. Fastrelax row ---
        fr13 = df_fr13[df_fr13["design"] == stats13_name]
        if not fr13.empty:
            fr_rec = fr13.iloc[0].to_dict()
            fr_rec["design"] = canonical
            # Update relaxed_pdb
            if os.path.exists(bmi_pdb):
                fr_rec["relaxed_pdb"] = bmi_pdb
            fr_rows.append(fr_rec)
        else:
            print(f"  WARNING: {stats13_name} not found in fastrelax_af3_scores.tsv")

        # --- f. Redesign audit row ---
        audit13 = df_a13[df_a13["design"] == stats13_name]
        if not audit13.empty:
            a_rec = audit13.iloc[0].to_dict()
            a_rec["design"] = canonical
            audit_rows.append(a_rec)
        else:
            print(f"  WARNING: {stats13_name} not found in redesign_audit.tsv")

    # ── 4. Build DataFrames and reorder columns ──────────────────────────────
    # Recompute pass_all_bindcraft correctly (source data skipped hydrophobicity check)
    BINDCRAFT = {
        "buns_delta_unsat":       4,    # < 4
        "surface_hydrophobicity": 0.35, # <= 0.35
        "interface_n_K":          3,    # <= 3
        "interface_n_M":          3,    # <= 3
    }

    def recompute_bindcraft(rec):
        buns  = rec.get("buns_delta_unsat",       float("nan"))
        hydro = rec.get("surface_hydrophobicity",  float("nan"))
        n_K   = rec.get("interface_n_K",           float("nan"))
        n_M   = rec.get("interface_n_M",           float("nan"))
        checks = []
        for val, thresh in [(buns, 4), (n_K, 3), (n_M, 3)]:
            if not math.isnan(float(val)):
                checks.append(float(val) < thresh)
        if not math.isnan(float(hydro)):
            checks.append(float(hydro) <= 0.35)
        return all(checks) if checks else None

    for rec in stats_rows:
        rec["pass_all_bindcraft"] = recompute_bindcraft(rec)
        rec["buns_pass"]          = (float(rec.get("buns_delta_unsat", float("nan"))) < 4
                                     if not math.isnan(float(rec.get("buns_delta_unsat", float("nan"))))
                                     else None)



    df_stats = pd.DataFrame(stats_rows)
    df_fr    = pd.DataFrame(fr_rows)
    df_audit = pd.DataFrame(audit_rows)

    # Move inward/outward columns after the standard p5 columns
    def reorder(df, anchor_cols, new_cols):
        """Insert new_cols after the last anchor_col that exists in df."""
        cols = list(df.columns)
        insert_after = -1
        for a in anchor_cols:
            if a in cols:
                insert_after = max(insert_after, cols.index(a))
        if insert_after == -1:
            return df
        present_new = [c for c in new_cols if c in cols]
        remaining   = [c for c in cols if c not in present_new]
        pos         = min(insert_after + 1, len(remaining))
        ordered     = remaining[:pos] + present_new + remaining[pos:]
        return df[[c for c in ordered if c in df.columns]]

    inward_outward_cols = [
        "asp_p5_chi1", "asp_p5_conformation",
        "cms_p5_inward", "cms_p5_outward",
        "n_contacts_p5_all_inward",   "n_contacts_p5_all_outward",
        "n_contacts_p5_polar_inward", "n_contacts_p5_polar_outward",
        "n_contacts_p5_hbond_inward", "n_contacts_p5_hbond_outward",
        "contacts_polar_p5_inward",   "contacts_polar_p5_outward",
        "asp_p5_chi1_outward_run", "asp_p5_conformation_outward_run",
        "ipsae_bp_outward_run", "cms_peptide_total_outward_run",
    ]
    df_stats = reorder(df_stats, ["cms_hotspot_p5", "contacts_polar_p5"], inward_outward_cols)

    # ── 5. Save outputs ──────────────────────────────────────────────────────

    stats_out   = os.path.join(OUT_DIR, "af3_design_stats.tsv")
    fr_out      = os.path.join(OUT_DIR, "fastrelax_af3_scores.tsv")
    audit_out   = os.path.join(OUT_DIR, "redesign_audit.tsv")

    df_stats.to_csv(stats_out, sep="\t", index=False)
    df_fr   .to_csv(fr_out,    sep="\t", index=False)
    df_audit.to_csv(audit_out, sep="\t", index=False)

    print(f"\nSaved {len(df_stats)} rows → {stats_out}")
    print(f"Saved {len(df_fr)}    rows → {fr_out}")
    print(f"Saved {len(df_audit)} rows → {audit_out}")

    # ── 6. Summary ───────────────────────────────────────────────────────────

    print("\n── Summary ─────────────────────────────────────────────────────────")
    sum_cols = [
        "design", "asp_p5_conformation",
        "ipsae_binder_peptide",
        "cms_hotspot_p5", "cms_p5_inward", "cms_p5_outward",
        "n_contacts_p5_polar_inward", "n_contacts_p5_polar_outward",
        "n_contacts_p5_hbond_inward", "n_contacts_p5_hbond_outward",
        "dG_separated", "packstat", "sc_value",
        "buns_delta_unsat", "surface_hydrophobicity",
        "pass_all_bindcraft", "pass_all",
    ]
    present = [c for c in sum_cols if c in df_stats.columns]
    print(df_stats[present].to_string(index=False))


if __name__ == "__main__":
    main()
