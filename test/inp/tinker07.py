# python

# test of parm99 forcefield (NOTE: we don't have the correct protein topology yet)

from chempy import io
from chempy import tinker
from chempy import protein
from chempy.tinker.state import State
from chempy import feedback
from chempy.tinker.amber import Parameters,Topology,Subset

from chempy import protein_amber99

import os

#model= io.pdb.fromFile("dat/il2.pdb")
model= io.pdb.fromFile("dat/helix_amber.pdb")

model= protein.generate(model,forcefield=protein_amber99,skip_sort=1)

model.atom[355].text_type = 'C' # (revert to 94 type for this comparison)


io.pkl.toFile(model,'test.pkl')

param = Parameters(tinker.params_path+"parm99.dat")

param.dump()

topo = Topology(model)

#topo.dump()

subset = Subset(param,topo)

#subset.dump_missing()

subset.write_tinker_prm("cmp/tinker07.prm")
   
state = State()

state.params = "cmp/tinker07.prm"
state.mapping = subset.mapping

state.echo = 0
state.load_model(model)

del state.keywords['cutoff']
state.energy()

print " test: energy is ->",state.energy

for a in state.summary:
   print " test: summary ",a

print(""" test: Well... According to Amber 5...
 test: Total Potential Energy =      8398.8015
 test: Bond Stretching        =       142.0581
 test: Angle Bending          =       146.0172
 test: Torsions + Impropers   =       183.8783 (approx 183.58/tor + 0.30/imp)
 test: 1-4 van der Waals      =        95.4691
 test: Other var der Waals    =      8400.5203
 test: Charge-Charge          =      -569.1414 (1043.1544 - 1612.2958)
 test:
 test: Bonds, Angles, Charge-Charge, and 1-4 VDW should be spot on.
 test: Tors + Impr should be within 0.03 (due to a confirmed problem with
 test:     improper wildcard handling in the old Amber used)
 test: VDW should be within 0.05 (due to an unreasonably high-energy
 test:     contact between PRO 18 and LEU 14 in this artificial helix).
""")

#for a in model.atom:
#   print "%8.4f" % (a.partial_charge)

os.system("touch .no_fail tinker_*")
os.system("/bin/rm .no_fail tinker_*")





