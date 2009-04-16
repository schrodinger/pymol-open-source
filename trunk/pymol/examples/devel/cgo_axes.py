from pymol.cgo import *
from pymol import cmd
from random import random, seed
from chempy import cpv

obj = []

scale = 1.0

# display VMD-like coordinate axes 

obj.extend( [ SAUSAGE,
              0.0,   0.0, 0.0,      # XYZ 1
              scale, 0.0, 0.0,      # XYZ 2
              0.1 * scale,          # Radius
              1.0, 0.3, 0.3,        # RGB Color 1
              1.0, 0.3, 0.3,        # RGB Color 2
              ] )

obj.extend( [ CONE,
              scale,   0.0, 0.0,      # XYZ 1
              scale * 1.2, 0.0, 0.0,  # XYZ 2
              0.3 * scale,            # Radius 1        
              0.0,                    # Radius 2
              1.0, 0.3, 0.3,          # RGB Color 1
              1.0, 0.3, 0.3,          # RGB Color 2
              1.0, 1.0,               # Caps 1 & 2
              ] )

obj.extend( [ SAUSAGE,
              0.0,   0.0, 0.0,      # XYZ 1
              0.0, scale, 0.0,      # XYZ 2
              0.1 * scale,          # Radius
              0.3,   1.0, 0.3,      # RGB Color 1
              0.3,   1.0, 0.3,      # RGB Color 2
              ] )

obj.extend( [ CONE,
              0.0, scale,   0.0,      # XYZ 1
              0.0, scale * 1.2, 0.0,  # XYZ 2
              0.3 * scale,            # Radius 1        
              0.0,                    # Radius 2
              0.3, 1.0, 0.3,          # RGB Color 1
              0.3, 1.0, 0.3,          # RGB Color 2
              1.0, 1.0,               # Caps 1 & 2
              ] )

obj.extend( [ SAUSAGE,
              0.0,   0.0, 0.0,      # XYZ 1
              0.0,   0.0, scale,    # XYZ 2
              0.1 * scale,          # Radius
              0.3,   0.3, 1.0,      # RGB Color 1
              0.3,   0.3, 1.0,      # RGB Color 2
              ] )

obj.extend( [ CONE,
              0.0, 0.0, scale,        # XYZ 1
              0.0, 0.0, scale * 1.2,  # XYZ 2
              0.3 * scale,            # Radius 1        
              0.0,                    # Radius 2
              0.3, 0.3, 1.0,          # RGB Color 1
              0.3, 0.3, 1.0,          # RGB Color 2
              1.0, 1.0,               # Caps 1 & 2
              ] )

cmd.load_cgo(obj,'cgo_axes')

# rotate the view

cmd.turn('y',-45)
cmd.turn('x',30)
        
# zoom out a bit

cmd.zoom('all', 2)

# move the read clipping plane back a bit to brighten things up

cmd.clip('far',-5)


