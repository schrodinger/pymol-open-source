from pymol.cgo import *
from pymol import cmd
from chempy.cpv import normalize, cross_product, add, scale
from math import sqrt

# pore through a membrane

radius = 19.0
edge = 70.0
sampling = 6
depth = 20.0

obj = []

for z1 in [0.0, depth]:
    for basis in [ [ [ 1.0,  0.0,  0.0], [  0.0,  1.0,  0.0], 0],
                   [ [ 0.0,  1.0,  0.0], [ -1.0,  0.0,  0.0], 0],
                   [ [-1.0,  0.0,  0.0], [  0.0, -1.0,  0.0], 0],
                   [ [ 0.0, -1.0,  0.0], [  1.0,  0.0,  0.0], 0],
                   [ [ 1.0,  0.0,  0.0], [  0.0,  1.0,  0.0], 1],
                   [ [ 0.0,  1.0,  0.0], [ -1.0,  0.0,  0.0], 1],
                   [ [-1.0,  0.0,  0.0], [  0.0, -1.0,  0.0], 1],
                   [ [ 0.0, -1.0,  0.0], [  1.0,  0.0,  0.0], 1],
                   ]:

        hand = basis[2]
        if hand:
            if z1 == 0.0:
                (x,y) = basis[0:2]
            else:
                (y,x) = basis[0:2]
            normal = normalize(cross_product(x,y))
        else:
            if z1 == 0.0:
                (y,x) = basis[0:2]
            else:
                (x,y) = basis[0:2]
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
            v0 = add( add( scale(x, x0 ), scale(y,y0) ), scale(normal,z1) )
            v1 = add( add( scale(x, x1 ), scale(y,y1) ), scale(normal,z1) )

            if hand:
                obj.extend( [ VERTEX ] + v0 + [ VERTEX ] + v1 )
            else:
                obj.extend( [ VERTEX ] + v1 + [ VERTEX ] + v0 )

        obj.extend( [END] )


        obj.extend( [BEGIN, TRIANGLE_STRIP] +
                    [COLOR, 1.0, 1.0, 1.0])

        for i in range(sampling+1):
            x1 = edge
            y1 = edge*i/sampling
            vlen = sqrt(x1*x1+y1*y1)
            n0 = add( scale(x,-x1/vlen), scale(y, -y1/vlen) )
            x0 = radius*x1/vlen
            y0 = radius*y1/vlen
            v0 = add( scale(x, x0), scale(y,y0) )
            v1 = add( add( scale(x, x0), scale(y,y0) ), scale(normal,z1) )
            obj.extend( [ NORMAL] + n0 )
            if hand:
                obj.extend( [ VERTEX ] + v0 + [ VERTEX ] + v1 )
            else:
                obj.extend( [ VERTEX ] + v1 + [ VERTEX ] + v0 )

        obj.extend( [END] )

        
        obj.extend( [BEGIN, TRIANGLE_STRIP] +
                    [COLOR, 1.0, 1.0, 1.0])

        for i in range(sampling+1):
            x1 = edge
            y1 = edge*i/sampling
            vlen = sqrt(x1*x1+y1*y1)
            n0 = scale(x, 1.0/edge) 
            v0 = add( scale(x, x1), scale(y,y1) )
            v1 = add( add( scale(x, x1), scale(y,y1) ), scale(normal,z1) )
            obj.extend( [ NORMAL] + n0 )
            if not hand:
                obj.extend( [ VERTEX ] + v0 + [ VERTEX ] + v1 )
            else:
                obj.extend( [ VERTEX ] + v1 + [ VERTEX ] + v0 )

        obj.extend( [END] )

# then we load it into PyMOL

cmd.load_cgo(obj,'cgo06')
                            
# position haemolysin through pore
if 1:
    
    cmd.fetch("7ahl",async=0)

    cmd.transform_selection("7ahl",
                            (0.70349078502033213, 0.19033556773811347, -0.68474258132841503, -11.993190038310882,
                             -0.15869362855324043, 0.98121371667870472, 0.10970571819737528, -29.939406659327624,
                             0.69276001846262825, 0.03148755047390385, 0.72048024614206196, -26.250180846754432,
                             0.0, 0.0, 0.0, 1.0))
                            
    cmd.show_as("cartoon","7ahl")
    cmd.do("util.cbc")
    cmd.set_view((\
    -0.589197695,   -0.440498680,    0.677344978,\
     0.574851871,   -0.817646086,   -0.031699535,\
     0.567794442,    0.370693058,    0.734975874,\
     0.000004128,   -0.000099868, -651.343872070,\
    -4.874452114,    8.360315323,   11.387301445,\
   528.689147949,  773.963684082,    0.000000000 ))
                 

