# -c

from pymol import cmd

/print "BEGIN-LOG"

load dat/small01.mol

save tmp/small01.mol
save tmp/small01.pkl
save tmp/small02.pkl,format=pkla
save tmp/small01.mmod
save tmp/small01.mmd
save tmp/small01.pdb

multisave tmp/small01.pmo,small01

dele all
load tmp/small01.pmo

dele all
load tmp/small01.pdb
count_atoms

dele all
load tmp/small01.mmd
count_atoms

dele all
load tmp/small01.mmod
count_atoms

dele all
load tmp/small01.pkl
count_atoms

dele all
load tmp/small02.pkl
count_atoms

dele all
load tmp/small01.mol
count_atoms

feedback disable,ray,details
ray
png tmp/tmp.png
load_png tmp/tmp.png

/print "END-LOG"





