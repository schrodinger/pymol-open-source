from pymol.cgo import *
from pymol import cmd

# this is a trivial example of creating a cgo object consisting of a
# a single state 

# first we create a list of floats containing a GL rendering stream

obj = [
   BEGIN, LINES,
   COLOR, 1.0, 1.0, 1.0,
   
   VERTEX, 0.0, 0.0, 0.0,
   VERTEX, 1.0, 0.0, 0.0,
   
   VERTEX, 0.0, 0.0, 0.0,
   VERTEX, 0.0, 2.0, 0.0,
   
   VERTEX, 0.0, 0.0, 0.0,
   VERTEX, 00, 0.0, 3.0,

   END
   ]

# then we load it into PyMOL

cmd.load_cgo(obj,'cgo01')

# move the read clipping plane back a bit to that that is it brighter

cmd.clip('far',-5)

