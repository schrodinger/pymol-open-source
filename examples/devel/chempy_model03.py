from pymol import cmd
from random import random
import time

# this shows how you can efficiently update the coordinates
# of an existing model for viewing as a PyMOL trajectory (from RAM)

# first we need a model

cmd.load("$PYMOL_PATH/test/dat/pept.pdb","demo")

# let's dress it up a little bit 

cmd.show("sticks","demo")

cmd.show("spheres","resi 10")

cmd.color("yellow","resi 5 and element c")

# now loop, updating the coordinates and appending the model
# onto 99 subsequent frames...

m = cmd.get_model()
for a in range(1,100):
   for a in m.atom:
      a.coord[0]+=(random()-0.5)*0.1
      a.coord[1]+=(random()-0.5)*0.1
      a.coord[2]+=(random()-0.5)*0.1
   cmd.load_model(m,"demo") # NOTE: no state number provided -> appends

# now define the movie with short pauses at beginning and and

cmd.mset("1 x15 1 -100 100 x15")

# now play the movie...

cmd.mplay()

# by default, PyMOL plays ~30 fps.
# "set movie_delay=0" to see maximum speed...

