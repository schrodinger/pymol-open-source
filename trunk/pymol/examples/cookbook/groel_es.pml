# run molscript on the PDB file

system molauto -nocentre 1AON.pdb | molscript -r > 1AON.r3d

# load raster3D input file

load 1AON.r3d

# zoom in a little

reset
move z,100

# render

set antialias=1
ray
png groel_es1.png

# render

turn y,90
ray
png groel_es2.png

