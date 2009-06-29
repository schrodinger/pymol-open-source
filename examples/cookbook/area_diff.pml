

load $TUT/1hpv.pdb

# get rid of solvent, ligands, etc.

remove not polymer

# measure as a complex

set dot_solvent,1

get_area all, load_b=1

# save areas in occupancy field

alter all, q=b

# then split apart into separate objects

extract chA, chain A

extract chB, chain B

# now re-measure

get_area chA, load_b=1

get_area chB, load_b=1

# subtract areas

alter all, b=b-q

# now color based on change in exposure

spectrum b, blue_red

# and splay apart

reset

orient

mset 1 x150

mview store, object=chA
mview store, object=chB
mview store, 15, object=chA
mview store, 15, object=chB

rotate y, -90, object=chA

translate [-15,0,0], object=chA

rotate y, 90, object=chB

translate [15,0,0], object=chB

mview store, 75, object=chA
mview store, 75, object=chB
mview store, 90, object=chA
mview store, 90, object=chB

zoom

show surface
set surface_solvent

mplay


