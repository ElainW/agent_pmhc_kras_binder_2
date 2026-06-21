from Bio.PDB import PDBParser, PDBIO, Select

parser = PDBParser(QUIET=True)
struct = parser.get_structure("9UV8", "input/9UV8.pdb")
model = struct[0]

io = PDBIO()

class CombineSelect(Select):
    def __init__(self):
        self.counter = 0
    def accept_chain(self, chain):
        return chain.id in ("A", "C")
    def accept_residue(self, residue):
        chain_id = residue.get_parent().id
        resnum = residue.id[1]
        if chain_id == "A" and 1 <= resnum <= 180:
            return True
        if chain_id == "C" and 1 <= resnum <= 9:
            return True
        return False

# Manually renumber+relabel since PDBIO.Select can't rewrite ids; build new structure
from Bio.PDB import Structure, Model, Chain, Residue

new_struct = Structure.Structure("9UV8_B")
new_model = Model.Model(0)
new_chain = Chain.Chain("B")

new_resnum = 1
for chain_id, lo, hi in [("A", 1, 180), ("C", 1, 9)]:
    chain = model[chain_id]
    for res in chain:
        resnum = res.id[1]
        if res.id[0] != " ":
            continue
        if not (lo <= resnum <= hi):
            continue
        new_res = Residue.Residue((" ", new_resnum, " "), res.resname, res.segid)
        for atom in res:
            new_res.add(atom.copy())
        new_chain.add(new_res)
        new_resnum += 1

new_model.add(new_chain)
new_struct.add(new_model)

io.set_structure(new_struct)
io.save("designs/prep/9UV8_target_chainB.pdb")
print("wrote designs/prep/9UV8_target_chainB.pdb, residues 1-%d (peptide = %d-%d)" % (new_resnum-1, 181, new_resnum-1))
