from pymol.cgo import *
from pymol import cmd
from chempy.cpv import normalize, cross_product, add, scale
from math import sqrt

# hole in a plane

radius = 2.0
edge = 10.0
sampling = 4


obj = []

for basis in [ [ [ 1.0,  0.0,  0.0], [  0.0,  1.0,  0.0], 0],
               [ [ 0.0,  1.0,  0.0], [ -1.0,  0.0,  0.0], 0],
               [ [-1.0,  0.0,  0.0], [  0.0, -1.0,  0.0], 0],
               [ [ 0.0, -1.0,  0.0], [  1.0,  0.0,  0.0], 0],
               [ [ 1.0,  0.0,  0.0], [  0.0,  1.0,  0.0], 1],
               [ [ 0.0,  1.0,  0.0], [ -1.0,  0.0,  0.0], 1],
               [ [-1.0,  0.0,  0.0], [  0.0, -1.0,  0.0], 1],
               [ [ 0.0, -1.0,  0.0], [  1.0,  0.0,  0.0], 1] ]:
    
    hand = basis[2]
    if hand:
        (x,y) = basis[0:2]
        normal = normalize(cross_product(x,y))
    else:
        (y,x) = basis[0:2]
        normal = normalize(cross_product(y,x))

    obj.extend( [BEGIN, TRIANGLE_STRIP] +
                [COLOR, 1.0, 1.0, 1.0] +
                [NORMAL] + normal )

    for i in range(sampling+1):
        x1 = edge
        y1 = edge*i/sampling
        vlen = sqrt(x1*x1+y1*y1)
        x0 = radius*x1/vlen
        y0 = radius*y1/vlen
        v0 = add( scale(x, x0 ), scale(y,y0) )
        v1 = add( scale(x, x1 ), scale(y,y1) )

        if hand:
            obj.extend( [ VERTEX ] + v0 + [ VERTEX ] + v1 )
        else:
            obj.extend( [ VERTEX ] + v1 + [ VERTEX ] + v0 )

    obj.extend( [END] )

# then we load it into PyMOL

cmd.load_cgo(obj,'cgo05')

# move the read clipping plane back a bit to that that is it brighter

cmd.clip('far',-5)

#cmd.set("two_sided_lighting")

