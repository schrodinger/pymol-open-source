# python

from chempy import io
from chempy import protein
from chempy import protein_amber99

model= io.pdb.fromFile("../../test/dat/pept.pdb")

model= protein.generate(model,forcefield=protein_amber99)

sm = 0
for a in model.atom:
   sm = sm + a.partial_charge

print " prot: net partial charge on protein is %8.4f" % sm
print " prot: (this should be integral)!"

io.pkl.toFile(model,"generate_amber.pkl")

