# python

from chempy import io
from chempy import protein
from chempy import protein_mmff
from chempy import bond_mmff

#
#print 'normal'
#protein_mmff.check_sum(protein_mmff.normal)
#print 'n_terminal'
#protein_mmff.check_sum(protein_mmff.n_terminal)
#print 'c_terminal'
#protein_mmff.check_sum(protein_mmff.c_terminal)
                       
model= io.pdb.fromFile("../../test/dat/pept.pdb")

model= protein.generate(model,forcefield=protein_mmff,bondfield=bond_mmff)

for a in model.atom:
   a.numeric_type = protein_mmff.alpha_map[a.text_type]
   
sm = 0
for a in model.atom:
   sm = sm + a.partial_charge

print " prot: net partial charge on protein is %8.3f" % sm
print " prot: (this should be integral)!"

io.pkl.toFile(model,"generate_mmff.pkl")

