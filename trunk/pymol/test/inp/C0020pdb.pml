# -c

/print "BEGIN-LOG"

load dat/pept.pdb

print cmd.get_names()
print cmd.select("(all)")

save cmp/C0020pdb.pdb

/print "END-LOG"
