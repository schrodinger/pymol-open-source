from pymol.cgo import *
from pymol import cmd
from random import random, seed
from chempy import cpv

# CGO ellipsoids

# first draw some walls

obj = [
   COLOR, 1.0, 1.0, 1.0,

   BEGIN, TRIANGLE_STRIP,
   NORMAL,  0.0,  0.0,  1.0,
   VERTEX,  0.0,  0.0,  0.0,   
   VERTEX, 10.0,  0.0,  0.0,
   VERTEX,  0.0, 10.0,  0.0,
   VERTEX, 10.0, 10.0,  0.0,
   END,

   BEGIN, TRIANGLE_STRIP,
   NORMAL,  1.0,  0.0,  0.0,
   VERTEX,  0.0,  0.0,  0.0,   
   VERTEX,  0.0, 10.0,  0.0,
   VERTEX,  0.0, 0.0,  10.0,  
   VERTEX,  0.0, 10.0, 10.0, 
   END,

   BEGIN, TRIANGLE_STRIP,
   NORMAL,  0.0,  1.0,  0.0,
   VERTEX,  0.0,  0.0,  0.0,   
   VERTEX,  0.0,  0.0, 10.0, 
   VERTEX, 10.0,  0.0,  0.0, 
   VERTEX, 10.0,  0.0, 10.0, 
   END
   
   ]

seed(0x1)

def random_conic(box, size, min_axis):

    # return a random ellipsoid record of the form:
    # [ ELLIPSOID, x_pos, y_pos, z_pos, size, x0, y0, z0, x1, y1, z2, x2, y2, z2 ]
    # where the xyz vectors are orthogonal and of length 1.0 or less.
    
    box = box - size
    tmp0 = [ size + random() * box, size + random() * box, size + random() * box ]
    tmp1 = cpv.random_vector()
    tmp2 = cpv.scale(tmp1,box/10)
    tmp1 = cpv.add(tmp2,tmp0)
    
    return [ CONE,
             tmp0[0], tmp0[1], tmp0[2], # coordinates 
             tmp1[0], tmp1[1], tmp1[2],
             (abs(random())*0.4+0.2) * size, # radii
             (abs(random())*0.1+0.01) * size,
             random(), random(), random(), # colors
             random(), random(), random(),
             1.0, 1.0 ]

for count in range(50):
    obj.extend( random_conic(10.0, 1.5, 0.2) )

# use more triangles when drawing ellipsoids

cmd.set('cgo_ellipsoid_quality', 2)

# then we load it into PyMOL

cmd.load_cgo(obj,'cgo08')

# rotate the view

cmd.turn('y',-45)
cmd.turn('x',30)
        
# zoom out a bit

cmd.zoom('all', 2)

# move the read clipping plane back a bit to brighten things up

cmd.clip('far',-5)


