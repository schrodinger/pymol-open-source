# -c

/print "BEGIN-LOG"

load dat/ligs3d.sdf

print cmd.get_names()
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

/print "END-LOG"
