import json
import sys

import numpy as np

pdb_in, json_in, out_prefix = sys.argv[1], sys.argv[2], sys.argv[3]

# Determine binder length (chain A) from the PDB.
binder_resnums = set()
lines_by_chain = {"A": [], "B": []}
for line in open(pdb_in):
    if line.startswith("ATOM"):
        chain = line[21]
        resnum = int(line[22:26])
        if chain == "A":
            binder_resnums.add(resnum)
        lines_by_chain.setdefault(chain, []).append(line)
binder_len = len(binder_resnums)
pep_start = binder_len + 181  # peptide = target residues 181-189 (our convention)
pep_end = binder_len + 189

# Write reduced PDB: chain A as-is, peptide residues relabeled to chain C, renumbered 1-9.
out_lines = []
for line in lines_by_chain["A"]:
    out_lines.append(line)
new_pep_num = {}
for line in lines_by_chain["B"]:
    resnum = int(line[22:26])
    if pep_start <= resnum <= pep_end:
        if resnum not in new_pep_num:
            new_pep_num[resnum] = len(new_pep_num) + 1
        nr = new_pep_num[resnum]
        newline = line[:21] + "C" + f"{nr:>4}" + line[26:]
        out_lines.append(newline)
with open(out_prefix + ".pdb", "w") as f:
    f.writelines(out_lines)
    f.write("END\n")

# Slice PAE/plddt: keep binder indices [0,binder_len) + peptide indices [pep_start-1, pep_end)
d = json.load(open(json_in))
pae = np.array(d["pae"])
plddt = np.array(d["plddt"])
keep = list(range(0, binder_len)) + list(range(pep_start - 1, pep_end))
pae_reduced = pae[np.ix_(keep, keep)]
plddt_reduced = plddt[keep]
with open(out_prefix + ".json", "w") as f:
    json.dump({"plddt": plddt_reduced.tolist(), "pae": pae_reduced.tolist()}, f)

print(f"binder_len={binder_len} peptide_orig_range=({pep_start},{pep_end}) reduced_n={len(keep)}")
