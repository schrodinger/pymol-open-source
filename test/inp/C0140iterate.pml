# -c

/print "BEGIN-LOG"

load dat/tiny.pdb


names = []

iterate (all),names.append(name)

print names

bfct = []

iterate (resn pro),bfct.append(b)
/for a in bfct: print "%8.4f"%a

# the following should have no effect
iterate (all),resn = 'NON'
iterate (all),b = b + 10

bfct = []
iterate (resn pro),bfct.append(b)
/for a in bfct: print "%8.4f"%a

/print "END-LOG"
