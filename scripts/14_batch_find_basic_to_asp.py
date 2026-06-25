import glob, os
import pyrosetta
from pyrosetta import pose_from_file

pyrosetta.init("-mute all -detect_disulf false -ignore_unrecognized_res true")

BASIC = {"ARG": ["NH1", "NH2", "NE"], "LYS": ["NZ"], "HIS": ["ND1", "NE2"]}

pdbs = sorted(glob.glob("/workspace/designs/round3/spec_scan_winners/*.pdb"))
print(f"scanning {len(pdbs)} backbones...")

hits = []
failed = 0
for i, path in enumerate(pdbs):
    tag = os.path.basename(path)[:-4]
    try:
        pose = pose_from_file(path)
    except Exception:
        failed += 1
        continue
    p5 = None
    for r in range(1, pose.size()+1):
        if pose.pdb_info().chain(r) == 'B' and pose.pdb_info().number(r) == 185:
            p5 = r
            break
    if p5 is None:
        continue
    p5res = pose.residue(p5)
    if p5res.name3() != 'ASP':
        continue
    sidechain_atoms = [bi for bi in range(1, p5res.natoms()+1) if p5res.atom_name(bi).strip() in ("OD1", "OD2")]
    binder_res = [r for r in range(1, pose.size()+1) if pose.pdb_info().chain(r) == 'A']

    best_d = 999
    best_info = None
    for r in binder_res:
        res = pose.residue(r)
        if res.name3() not in BASIC:
            continue
        if (res.xyz("CA") - p5res.xyz("CA")).norm() > 15:
            continue
        for atom_name in BASIC[res.name3()]:
            if not res.has(atom_name):
                continue
            xa = res.xyz(atom_name)
            for bi in sidechain_atoms:
                d = (xa - p5res.xyz(bi)).norm()
                if d < best_d:
                    best_d = d
                    best_info = (pose.pdb_info().number(r), res.name3(), atom_name, p5res.atom_name(bi).strip())
    if best_d <= 6.0:
        hits.append((tag, best_d, best_info))
    if (i+1) % 200 == 0:
        print(f"  ...{i+1}/{len(pdbs)} scanned, {len(hits)} hits so far")

print(f"\nTotal: {len(hits)} hits out of {len(pdbs)} backbones scanned ({failed} failed to load)")
for tag, d, info in sorted(hits, key=lambda x: x[1]):
    print(f"  {tag}: {info[1]}{info[0]}:{info[2]} <-> p5:{info[3]} dist={d:.2f}A")
