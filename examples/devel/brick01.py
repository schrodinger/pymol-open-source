


# ===============================
# brick creation

from chempy.brick import Brick

brik = Brick()

brik.setup_from_min_max(
   [0.0,0.0,0.0],
   [2.0,2.0,2.0],
   [0.2,0.2,0.2])

# ===============================
# brick data

import math

for x in range(11):
   for y in range(11):   
      for z in range(11):
          mx=1.5*(x-5)
          my=2*(y-5)
          mz=3*(z-5)
          brik.lvl[x,y,z]=1/math.exp(math.sqrt(mx*mx+my*my+mz*mz)/10.0)

# ===============================
# easy pickling code

#from chempy import io

#io.pkl.toFile(brik,"brick.pkb")

#brik = io.pkl.fromFile("brick.pkb")

# ===============================
# for display

#from pymol import cmd

cmd.load_brick(brik,"brick")

for a in range(6):
	lvl = a/10.0+0.3
	cmd.do("isomesh lvl%1.0f=brick,%0.1f" %(lvl*10,lvl))

cmd.disable()
cmd.enable("brick")
cmd.enable("lvl6")
cmd.color("red","brick")
