# -c

/print "BEGIN-LOG"

load dat/tiny.pdb
print cmd.get_names()

show sticks
refresh

show spheres
refresh

show dots
refresh

show surface
refresh

show mesh
refresh

/print "END-LOG"
