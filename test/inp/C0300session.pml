# -c

/print "BEGIN-LOG"

load dat/pept.pdb
select cas,name ca
set line_width,5
show spheres,name n

save tmp/session.pse
reinitialize

print cmd.get_names()
get line_width
count_atoms rep spheres

load tmp/session.pse

print cmd.get_names()
get line_width
count_atoms rep spheres

count_atoms cas

/print "END-LOG"
