# -c

/print "BEGIN-LOG"

load dat/small01.mol

print cmd.get_names()
print cmd.select("(all)")

/print "END-LOG"
