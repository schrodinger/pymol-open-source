# -c

/print "BEGIN-LOG"

load dat/small03.mol2

print cmd.get_names()
print cmd.count_states()
print cmd.select("(all)")

set ray_default_renderer=2

set valence

frame 1
ray
frame 10
ray

hide
show sticks
ray
unset valence
ray

dele all

load dat/small03.mol2, discrete=0
count_atoms
dele all

load dat/small03.mol2, multiplex=1
print cmd.get_names()
print cmd.count_states()
print cmd.select("(all)")

dele all
set multiplex, 1
load dat/small03.mol2
print cmd.get_names()
print cmd.count_states()
print cmd.select("(all)")


/print "END-LOG"
