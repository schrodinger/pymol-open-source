from pymol.cgo import *
from pymol import cmd
import math

# this is an example of creating a cgo object consisting of multiple
# states which form a movie

for a in range(1,62):

   # create a state with some unique geometry

   obj = [
      BEGIN, LINES,
      COLOR, 1.0, 1.0, 1.0,

      VERTEX, 0.0, 0.0, 0.0,
      VERTEX, 1.0, 0.0, 0.0,

      VERTEX, 0.0, 0.0, 0.0,
      VERTEX, 0.0, 2.0, 0.0,

      VERTEX, 0.0, 0.0, 0.0,
      VERTEX, 00, 0.0, 3.0,

      END,

      BEGIN, TRIANGLES,

      COLOR, 1.0, 0.2, 0.2,
      
      VERTEX, 0.0, 0.0, 0.0,
      VERTEX, math.cos((a-1)/10.0), math.sin((a-1)/10.0), 2.0,
      VERTEX, 0.0, 0.0, 2.0,
      
      END
      ]

   # load state into PyMOL

   cmd.load_cgo(obj,'cgo02',a)

# this moves the rear clipping plane back a bit

cmd.clip('far',-5)

# now, start the movie

cmd.mplay()
   
   
