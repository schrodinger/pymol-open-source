from pymol import cmd
from pymol.cgo import *

# This example shows one way to create a plane behind the current
# molecule.
#
# Note that because PyMOL's ray-tracer doesn't currently support
# perspective, the plane will not look right if the edges are
# showing. Thus, it is best to zoom the image so that the edges can't
# be seen and that the plane appears infinite
#
# To use this script, setup your molecule and then "run cgo_plane.py".
# This will create a plane in space about 80% of the way to the far
# clipping plane.  You can then rotate the camera around to get the
# desired shadowing.
#
# If the plane is too close to the molecule, move the rear clipping
# plane back and then re-run cgo_plane.py.  
#
# NOTE that once the plane is created, there is no easy way to move it
# (other than recreating it).  However, you can move you molecule
# relative to the plane using PyMOL's molecular editing features.
# 
# (Remember, you can't click on cartoons, but you can click on ribbons).
#
# Good luck,
# Warren

view = cmd.get_view()

# locate plane most of the way to the rear clipping plane

plane_z = - (view[11] + (view[15]+5*view[16])/6.0)

# choose the size of the plane

plane_size = abs(view[11])/3.0
   
obj = []

# now create a plane in camera space

plane = [
   [ -plane_size,  plane_size, plane_z ],
   [  plane_size,  plane_size, plane_z ],
   [ -plane_size, -plane_size, plane_z ],
   [  plane_size, -plane_size, plane_z ]]


normal = [ 0.0, 0.0, 1.0 ]


# then transform plane coordinates into model space 

plane = map( lambda p,v=view: [
   v[0] * p[0] + v[1] * p[1] + v[2]* p[2],
   v[3] * p[0] + v[4] * p[1] + v[5]* p[2],
   v[6] * p[0] + v[7] * p[1] + v[8]* p[2]
   ], plane )

normal = apply( lambda p,v=view:[
   v[0] * p[0] + v[1] * p[1] + v[2]* p[2],
   v[3] * p[0] + v[4] * p[1] + v[5]* p[2],
   v[6] * p[0] + v[7] * p[1] + v[8]* p[2]
   ], (normal,) )

# and position relative to the camera 

plane = map( lambda p,v=view: [
   p[0] + v[9 ] + v[12],
   p[1] + v[10] + v[13],
   p[2] +       + v[14],
   ], plane )

# set color
obj.extend( [ COLOR, 0.8, 0.8, 0.8 ] ) # greyish

# begin triangle strip 
obj.extend( [ BEGIN, TRIANGLE_STRIP ] )

# set normal

obj.append( NORMAL )
obj.extend( normal )

# draw the plane
for a in plane:
   obj.append( VERTEX)
   obj.extend(a)
obj.append( END )

# delete existing object (if any)
cmd.delete("cgo_plane")

# now load the new object without zooming
auto_zoom = cmd.get('auto_zoom')
cmd.set('auto_zoom', 0, quiet=1)
cmd.load_cgo(obj,'cgo_plane')
cmd.set('auto_zoom', auto_zoom, quiet=1)

