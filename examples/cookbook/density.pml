# alter settings for publication quality

set mesh_radius = 0.02
set antialias = 1
set stick_radius = 0.1

# load pdb and map file

load 1DN2.pdb,1dn2
load 2fofc.xplor

# display region of interest

hide
show sticks,(byres ((c;f & i;5,6) x;4))

# show e-density nearby

isomesh den, 2fofc, 1.0, (c;F & i;6), 8.0
color marine,den

# zoom in, turn, and clip

zoom (c;f & i;5,6),2
turn x,10
turn y,-10
clip slab, 6.5

ray
png density.png

