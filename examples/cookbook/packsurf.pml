# load pdb file which has a CRYST1 record

load 1DN2.pdb,1dn2

# remove waters 

remove (resn hoh)

# create neighbors with contacts

symexp s,1dn2,(all),4.0

# display surfaces 

hide 
show surface

# rotate and zoom

reset
turn y,90
move z,30
turn x,5
turn y,5

# color symmetry-related copies distinctly

color white
color yellow, (1dn2)
color cyan,   (s01-10100)
color marine, (s0000-100 or s00000100)
color red,    (s01000000)
color violet, (s00000001 or s000000-1)
color orange, (s01000001)
color pink,   (s01000101)
color green,  (s01000100)

# now render 

set antialias=1
ray
png packsurf.png
