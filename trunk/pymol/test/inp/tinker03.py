# python

# dynamics test

from chempy import io
from chempy import tinker
from chempy import protein
from chempy.tinker.state import State
from chempy import feedback

import os

m = io.pdb.fromFile("dat/pept.pdb")

m = protein.generate(m)

s = State()

s.echo = 0

s.load_model(m)
print " test: atom 0 position:",m.atom[0].coord
s.dynamics(steps=50,interval=10)

print " test:",len(s.frames),"frames returned."
print " test: atom 0 position:",m.atom[0].coord

s.dynamics(steps=50,interval=10)
print " test: atom 0 position:",m.atom[0].coord

print " test: energy is ->",s.energy

for a in s.summary:
   print " test: summary ",a

os.system("touch .no_fail tinker_*")
os.system("/bin/rm .no_fail tinker_*")

