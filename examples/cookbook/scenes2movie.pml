# first, create some scenes

load $TUT/1hpv.pdb
util.cbc

orient
as ribbon
scene new, store

orient organic
show sticks, organic
scene new, store

show lines, byres organic around 5
orient rep lines
scene new, store

# now define a movie

mset 1 x180

# here's what our movie will look like:

# 1-30, scene 001 static
# 31-60, scene 002 interpolation
# 61-90, scene 002 static
# 91-120, scene 003 interpolation
# 121-150, scene 003 static
# 151-180, scene 001 interpolation

# first we define the scene content

mdo 1: scene 001, view=0, animate=0, quiet=1
mdo 31: scene 002, view=0, animate=0, quiet=1
mdo 91: scene 003, view=0, animate=0, quiet=1
mdo 151: scene 001, view=0, animate=0, quiet=1

# then we define the camera waypoints

scene 001, animate=0
mview store, 1
mview store, 31
mview store, 180

scene 002, animate=0
mview store, 61
mview store, 91

scene 003, animate=0
mview store, 121
mview store, 151

# now we interpolate

mview interpolate

# and play the movie

rewind
mplay

