# -c

/print "BEGIN-LOG"

load dat/pept.pkl

print cmd.get_names()
print cmd.select("(all)")

save cmp/C0030pkl.pdb

/print "END-LOG"
