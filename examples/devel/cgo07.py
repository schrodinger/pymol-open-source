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

def random_ellipsoid(box, size, min_axis):

    # return a random ellipsoid record of the form:
    # [ ELLIPSOID, x_pos, y_pos, z_pos, size, x0, y0, z0, x1, y1, z2, x2, y2, z2 ]
    # where the xyz vectors are orthogonal and of length 1.0 or less.
    
    box = box - size
    tmp0 = cpv.random_vector()
    tmp1 = cpv.random_vector()
    tmp2 = cpv.cross_product(tmp1, tmp0)
    tmp3 = cpv.cross_product(tmp1, tmp2)
    tmp4 = cpv.cross_product(tmp2, tmp3)
    tmp2 = cpv.normalize(tmp2)
    tmp3 = cpv.normalize(tmp3)
    tmp4 = cpv.normalize(tmp4)
    primary = cpv.scale(tmp2, random())
    secondary = cpv.scale(tmp3,random())
    tertiary = cpv.scale(tmp4,random())
    factor = 1.0 / max( cpv.length(primary), cpv.length(secondary), cpv.length(tertiary))
    primary = cpv.scale(primary, factor)
    secodary = cpv.scale(secondary, factor)
    tertiary = cpv.scale(tertiary, factor)
    return [ ELLIPSOID,
             size + random() * box, size + random() * box, size + random() * box,
             max(random() * size, min_axis),
             ] + primary + secondary + tertiary
             

for count in range(100):
#    obj.extend( [ALPHA, random() ] )
    obj.extend( [COLOR, random(), random(), random()] )
    obj.extend( random_ellipsoid(10.0, 1.5, 0.2) )

# use more triangles when drawing ellipsoids

cmd.set('cgo_ellipsoid_quality', 2)

# then we load it into PyMOL

cmd.load_cgo(obj,'cgo01')

# rotate the view

cmd.turn('y',-45)
cmd.turn('x',30)
        
# zoom out a bit

cmd.zoom('all', 2)

# move the read clipping plane back a bit to brighten things up

cmd.clip('far',-5)


