from pymol import cmd
from chempy.models import Indexed
from chempy import Atom,Bond

# here's how you build a chemical python model from scratch...

# you can run this script as follows:
#    "pymol chempy_model01.py" from unix
#    "run chempy_model01.py" from within PyMOL

# first we need some raw material -- this stuff could come from
# anywhere, but I'll just start with some Python lists...

atoms = [
   ['C', 'CH3' ,'1','ACE'],
   ['C', 'C'   ,'1','ACE'],
   ['O', 'O'   ,'1','ACE'],
   ['H', '1HH3','1','ACE'],
   ['H', '2HH3','1','ACE'],
   ['H', '3HH3','1','ACE'],
   ['N', 'N'   ,'2','ALA'],
   ['C', 'CA'  ,'2','ALA'],
   ['C', 'CB'  ,'2','ALA'],
   ['C', 'C'   ,'2','ALA'],
   ['O', 'O'   ,'2','ALA'],
   ['H', '1HB' ,'2','ALA'],
   ['H', '2HB' ,'2','ALA'],
   ['H', '3HB' ,'2','ALA'],
   ['H', 'H'   ,'2','ALA'],
   ['H', 'HA'  ,'2','ALA'],
   ['N', 'N'   ,'3','NME'],
   ['C', 'CH3' ,'3','NME'],
   ['H', '1HH3','3','NME'],
   ['H', '2HH3','3','NME'],
   ['H', '3HH3','3','NME'],
   ['H', 'H'   ,'3','NME'],
   ]

coords = [
   [  -1.862, -11.290,  18.059],
   [  -0.948, -10.106,  18.286],
   [  -1.142,  -9.030,  17.704],
   [  -1.851, -11.572,  17.015],
   [  -1.532, -12.131,  18.651],
   [  -2.873, -11.038,  18.343],
   [   0.186, -10.234,  19.211],
   [   0.862,  -8.940,  19.211],
   [   0.354,  -8.148,  20.428],
   [   2.362,  -9.114,  19.211],
   [   2.929,  -9.926,  19.953],
   [  -0.742,  -7.998,  20.392],
   [   0.578,  -8.662,  21.383],
   [   0.810,  -7.143,  20.485],
   [   0.696, -11.147,  18.941],
   [   0.594,  -8.401,  18.283],
   [   3.174,  -8.292,  18.302],
   [   4.560,  -8.671,  18.519],
   [   4.679,  -9.084,  19.509],
   [   4.859,  -9.412,  17.791],
   [   5.199,  -7.806,  18.421],
   [   2.880,  -7.315,  18.440],
   ]

bonds = [
   [0,1,1],
   [0,3,1],
   [0,4,1],
   [0,5,1],
   [1,2,2],
   [6,1,1],
   [6,7,1],
   [6,14,1],
   [7,8,1],
   [7,9,1],
   [7,15,1],
   [9,10,2],
   [9,16,1],
   [11,8,1],
   [12,8,1],
   [13,8,1],
   [16,17,1],
   [16,21,1],
   [17,18,1],
   [17,19,1],
   [17,20,1],
   ]

# okay, now we'll build the object from scratch...

# create a model instance

model = Indexed()

# append the atoms onto it 

for a in atoms:
   new_atom = Atom()
   new_atom.symbol = a[0]     # elemental symbol
   new_atom.name = a[1]       # atom name
   new_atom.resi = a[2]       # residue identifier
   new_atom.resn = a[3]       # residue name
   model.atom.append(new_atom)
# (note that there are a bunch of other fields we're not using -- and none are required)

# add coordinates onto the atoms

for a in model.atom: # now assign coordinates
   a.coord = coords.pop(0)

# now specify the bonds

for a in bonds:
   new_bond = Bond()
   new_bond.index = [a[0],a[1]]  # atom indices (zero-based)
   new_bond.order = a[2]         # bond order
   model.bond.append(new_bond)

# finally, load the model into PyMOL

cmd.load_model(model,"example")


   


