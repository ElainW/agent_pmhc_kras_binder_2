#!/usr/bin/env python
"""
Calibrate Rosetta interface-metric cutoffs from the experimentally VALIDATED designs.

Run this on the GPU host inside the PyRosetta conda env (see envs_export/). It
loads each validated binder-pMHC model in paper_output/af3_nomsa/ (and optionally
the 9O5S crystal), runs InterfaceAnalyzer for the binder-vs-(MHC+peptide)
interface, and reports ddG / shape-complementarity / packstat / buried-unsat
H-bonds + a contact-molecular-surface (CMS) proxy on peptide position 5.

The resulting distribution across the 19 hits IS the reference: set the scoring-plan
cutoffs (docs/02) to roughly the 10th-percentile (permissive) / median (target) of
what validated binders achieve, instead of guessing.

Usage (PyRosetta lives in the `dl_binder_design` conda env):
    conda activate dl_binder_design
    python scripts/rosetta_reference_metrics.py \
        --inputs 'post_filter/outputs/author_design_stats/fastrelax_af3_nomsa/relaxed_pdbs/*_relaxed.pdb' \
        --out designs/validated_reference_metrics.csv
    # raw AF3 cifs also work:  --inputs 'paper_output/af3_nomsa/*/*_model.cif'
    # optionally also:         --crystal paper_output/mmdb_9O5S.pdb --crystal_binder_chain Q

Chain detection (AF3 models): peptide = shortest chain (8-11 res); binder =
50-120 res and NOT ~99 (that's beta-2 microglobulin); everything else = "target".
Verify on a couple of structures before trusting the batch.
"""
import argparse, glob, os, re, csv, statistics as st
import pyrosetta
from pyrosetta import pose_from_file
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.relax import FastRelax

pyrosetta.init("-mute all -detect_disulf false -ignore_unrecognized_res true")

# Raw AF3 models (no side-chain relaxation/repacking) score with huge spurious ddG from
# clashes -- the validated-design reference numbers were computed on FastRelax'd structures
# (see --inputs default above), so raw AF3 cif input needs --relax for comparable numbers,
# rather than silently giving bogus ddG.
_FASTRELAX = None


def get_fastrelax():
    global _FASTRELAX
    if _FASTRELAX is None:
        scorefxn = pyrosetta.get_fa_scorefxn()  # same default IA/analyze() implicitly use
        fr = FastRelax()
        fr.set_scorefxn(scorefxn)
        fr.max_iter(200)
        _FASTRELAX = fr
    return _FASTRELAX


def relax_pose(pose):
    get_fastrelax().apply(pose)


def chain_lengths(pose):
    """Return {pdb_chain_id: n_residues} using pose chain order."""
    info = {}
    for i in range(1, pose.num_chains() + 1):
        begin = pose.chain_begin(i)
        cid = pose.pdb_info().chain(begin)
        info[cid] = pose.chain_end(i) - begin + 1
    return info


def classify_chains(pose):
    lens = chain_lengths(pose)
    items = sorted(lens.items(), key=lambda kv: kv[1])
    pep = items[0][0]                      # shortest = peptide
    binder = None
    for cid, n in items[1:]:
        if 50 <= n <= 120 and not (95 <= n <= 101):   # skip beta-2-microglobulin (~99)
            binder = cid
            break
    if binder is None:
        binder = items[1][0]
    target = [c for c in lens if c != binder]   # MHC (+ b2m + peptide) = everything else
    return binder, pep, target, lens


def cms_proxy_p5(pose, binder, pep):
    """Lightweight CMS proxy: count binder heavy-atom contacts (<4.5 A) to the
    5th residue of the peptide chain. Replace with the repo contact_filter CMS for
    the production number; this is just for relative calibration here."""
    pep_res = [r for r in range(1, pose.size() + 1)
               if pose.pdb_info().chain(r) == pep]
    if len(pep_res) < 5:
        return None
    p5 = pep_res[4]
    p5res = pose.residue(p5)
    contacts = 0
    for r in range(1, pose.size() + 1):
        if pose.pdb_info().chain(r) != binder:
            continue
        res = pose.residue(r)
        # quick CA-CA prefilter
        if (res.xyz("CA") - p5res.xyz("CA")).norm() > 14:
            continue
        for ai in range(1, res.natoms() + 1):
            if res.atom_type(ai).element() == "H":
                continue
            xa = res.xyz(ai)
            hit = False
            for bi in range(1, p5res.natoms() + 1):
                if p5res.atom_type(bi).element() == "H":
                    continue
                if (xa - p5res.xyz(bi)).norm() <= 4.5:
                    hit = True
                    break
            if hit:
                contacts += 1
                break
    return contacts


def analyze(pose, binder, target):
    ia = InterfaceAnalyzerMover(f"{binder}_{''.join(target)}")
    ia.set_pack_separated(True)
    ia.set_compute_packstat(True)
    ia.set_compute_interface_sc(True)                 # required, else sc_value is unset
    ia.apply(pose)
    d = ia.get_all_data()
    return dict(
        ddG=ia.get_separated_interface_energy(),     # dG_separated (REU)
        sc=d.sc_value,                                # shape complementarity
        packstat=d.packstat,
        buns=d.delta_unsat_hbonds,                    # delta unsatisfied H-bonds
        dSASA=ia.get_interface_delta_sasa(),
        n_int_res=ia.get_num_interface_residues(),
    )


def main():
    ap = argparse.ArgumentParser()
    # Default = the FastRelax'd AF3 no-MSA PDBs (what we actually scored). Point
    # --inputs at any glob of relaxed PDBs (recommended) or raw AF3 *_model.cif.
    ap.add_argument("--inputs",
                    default="post_filter/outputs/author_design_stats/"
                            "fastrelax_af3_nomsa/relaxed_pdbs/*_relaxed.pdb",
                    help="glob of relaxed *.pdb (preferred) or AF3 *_model.cif")
    ap.add_argument("--out", default="designs/validated_reference_metrics.csv")
    ap.add_argument("--crystal")
    ap.add_argument("--crystal_binder_chain")
    ap.add_argument("--relax", action="store_true",
                    help="FastRelax each pose before scoring -- required for raw AF3 cifs "
                         "(unrelaxed clashes otherwise give nonsense ddG); not needed if "
                         "--inputs already points at pre-relaxed PDBs.")
    args = ap.parse_args()
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    rows = []
    models = sorted(glob.glob(args.inputs))
    for path in models:
        name = re.sub(r"(_relaxed|_model)$", "",
                      os.path.splitext(os.path.basename(path))[0])
        try:
            pose = pose_from_file(path)
            if args.relax:
                relax_pose(pose)
            binder, pep, target, lens = classify_chains(pose)
            m = analyze(pose, binder, target)
            m["cms_p5_proxy"] = cms_proxy_p5(pose, binder, pep)
            m.update(design=name, binder_chain=binder, pep_chain=pep,
                     binder_len=lens[binder])
            rows.append(m)
            print(f"{name:12} binder={binder}({lens[binder]}) ddG={m['ddG']:.1f} "
                  f"sc={m['sc']:.2f} packstat={m['packstat']:.2f} "
                  f"buns={m['buns']} cms_p5={m['cms_p5_proxy']}")
        except Exception as e:
            print(f"{name:12} FAILED: {e}")

    if args.crystal and args.crystal_binder_chain:
        pose = pose_from_file(args.crystal)
        target = [c for c in chain_lengths(pose) if c != args.crystal_binder_chain]
        m = analyze(pose, args.crystal_binder_chain, target)
        m.update(design="9O5S_crystal", binder_chain=args.crystal_binder_chain)
        rows.append(m)

    cols = ["design", "binder_chain", "binder_len", "pep_chain", "ddG", "sc",
            "packstat", "buns", "dSASA", "n_int_res", "cms_p5_proxy"]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    # print the calibration summary
    def stat(key):
        vals = [r[key] for r in rows if isinstance(r.get(key), (int, float))]
        if not vals:
            return "n/a"
        vals.sort()
        p10 = vals[max(0, int(0.1 * len(vals)) - 1)]
        return f"min={vals[0]:.2f} p10={p10:.2f} median={st.median(vals):.2f} max={vals[-1]:.2f}"
    print("\n=== Validated-design reference distributions (use to set cutoffs) ===")
    for k in ["ddG", "sc", "packstat", "buns", "dSASA", "cms_p5_proxy"]:
        print(f"  {k:12} {stat(k)}")
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
