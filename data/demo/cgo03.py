from pymol.cgo import *
from pymol import cmd
import math

# a more comple example of using cgo's to create a fancy animation.

# In a real script you'd want to use Numeric Python in order to handle
# the math efficiently.

# create a 30x30 array

z = []
for x in range(0,30):
   z.append([0.0]*30)

# create a 63 frame movie inside of PyMOL

for a in range(0,63):

   aa = a/9.973
   ab = a/4.987
   
   # blank list
   
   obj = []

   # generate a flowing mesh
   
   for x in range(0,30):
      for y in range(0,30):
         z[x][y]=(math.cos((aa)+x/2.0)+math.sin((aa)+(y/2.0)))

   obj.extend( [ COLOR, 0.3, 0.3, 1.0 ] )
   
   for x in range(0,30):
      obj.extend( [ BEGIN, LINE_STRIP ] )
      for y in range(0,30):
         obj.extend( [ VERTEX, float(x), float(y), z[x][y] ] )
      obj.append( END )

   for y in range(0,30):
      obj.extend( [ BEGIN, LINE_STRIP ] )
      for x in range(0,30):
         obj.extend( [ VERTEX, float(x), float(y), z[x][y] ] )
      obj.append( END )

   # add into this a couple circulating spheres

   obj.extend( [ COLOR, 1.0, 0.2, 0.2 ] )
   obj.extend( [ SPHERE, 5*math.cos(ab)+15.0, 5*math.sin(ab)+15.0, 6.0, 1.0 ] )

   obj.extend( [ COLOR, 1.0, 1.0, 0.2 ] )
   obj.extend( [ SPHERE, 5*math.cos(aa)+15.0, 5*math.sin(aa)+15.0, 3.0, 2.0 ] )

   # now add a colorful cylinder

   obj.extend( [ CYLINDER,
                 5.0,  15+math.sin(aa)*10, -5.0,      # XYZ 1
                 25.0, 15+math.sin(aa)*10, -5.0,      # XYZ 2
                 2.0,                                     # Radius
                 1.0, (1.0+math.sin(aa))/2.0,(1.0+math.cos(aa))/2.0, # RGB Color 1
                 0.3, (1.0+math.cos(aa))/2.0, 0.5,              # RGB Color 2
                 ] )
                
   # load this state into the PyMOL object

   cmd.load_cgo(obj,'cgo03',a+1)

# this zooms out a bit more than usual

cmd.reset()
cmd.zoom('cgo03',3.0)

# this moves the rear clipping plane back a bit

cmd.clip('far',-12.0)

# give us a nice view

cmd.turn('z',30)
cmd.turn('x',-60)

# start the movie

cmd.mplay()
   
   
