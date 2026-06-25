import sys
import pyrosetta
from pyrosetta import pose_from_file

pyrosetta.init("-mute all -detect_disulf false -ignore_unrecognized_res true")

path = sys.argv[1]
tag = sys.argv[2] if len(sys.argv) > 2 else path
two_chain = "--two_chain" in sys.argv

pose = pose_from_file(path)

if two_chain:
    p5 = None
    for r in range(1, pose.size()+1):
        if pose.pdb_info().chain(r) == 'B' and pose.pdb_info().number(r) == 185:
            p5 = r
            break
    binder_res = [r for r in range(1, pose.size()+1) if pose.pdb_info().chain(r) == 'A']
else:
    chain_c = [r for r in range(1, pose.size()+1) if pose.pdb_info().chain(r) == 'C']
    p5 = chain_c[4]
    binder_res = [r for r in range(1, pose.size()+1) if pose.pdb_info().chain(r) == 'A']

p5res = pose.residue(p5)
if p5res.name3() != 'ASP':
    print(f"{tag}: p5 is not ASP ({p5res.name3()}) -- skip")
    sys.exit()

sidechain_atoms = [bi for bi in range(1, p5res.natoms()+1) if p5res.atom_name(bi).strip() in ("OD1", "OD2")]
BASIC = {"ARG": ["NH1", "NH2", "NE"], "LYS": ["NZ"], "HIS": ["ND1", "NE2"]}

hits = []
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
            if d <= 6.0:
                hits.append((pose.pdb_info().number(r), res.name3(), atom_name, p5res.atom_name(bi).strip(), round(d,2)))

if hits:
    best = min(hits, key=lambda x: x[4])
    print(f"{tag}: HIT -- {best[1]}{best[0]}:{best[2]} <-> p5:{best[3]} dist={best[4]}A  (all: {hits})")
else:
    print(f"{tag}: no Arg/Lys/His side-chain atom within 6.0A of p5 OD1/OD2")
