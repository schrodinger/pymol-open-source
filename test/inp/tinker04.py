# python

# test of chempy forcefield for comparison with the others.

from chempy import io
from chempy import tinker
from chempy import protein
from chempy.tinker.state import State
from chempy import feedback

import os

#model= io.pdb.fromFile("dat/il2.pdb")
model= io.pdb.fromFile("dat/pept.pdb")

model = protein.generate(model)

state = State()

state.echo = 0

state.load_model(model)
print " test: atom 0 position:",model.atom[0].coord
state.energy(kw=["debug\n"])

print " test: energy is ->",state.energy

for a in state.summary:
   print " test: summary ",a

os.system("touch .no_fail tinker_*")
os.system("/bin/rm .no_fail tinker_*")

