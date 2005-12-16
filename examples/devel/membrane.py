from pymol import cmd
from chempy.models import Indexed
from chempy import Atom,Bond
import math


# create a model instance

model = Indexed()

# append the atoms onto it 

for x in range(-63,63,1):
   for y in range(-63,63,1):
      new_atom = Atom()
      new_atom.symbol = 'O'
      new_atom.coord = [ x*2, y*2, 30 * math.cos(math.sqrt(x*x+y*y)/60.0) ]
      model.atom.append(new_atom)


cmd.load_model(model,"membrane")
cmd.hide("everything","membrane")
cmd.show("spheres","membrane")
cmd.color("lightblue","membrane")

cmd.set_view( (\
     0.736728907,   -0.144400939,    0.660589039,\
     0.675238073,    0.208899528,   -0.707400322,\
    -0.035847016,    0.967217624,    0.251406968,\
     0.000008686,    0.000009686, -332.961212158,\
    -2.366872311,   -1.122793436,   23.127344131,\
    11.627288818,  654.294433594,    0.000000000 ))


# uncomment this if you programmable shaders
# cmd.set("sphere_mode",5)

# uncomment if you have lotsa RAM
# cmd.set("hash_max",240)
