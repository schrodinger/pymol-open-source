# -c

/print "BEGIN-LOG"

load dat/tiny.pdb


names = []

iterate (all),names.append(name)

print names

bfct = []

iterate (resn pro),bfct.append(b)
print bfct

# the following should have no effecit
iterate (all),resn = 'NON'
iterate (all),b = b + 10

bfct = []
iterate (resn pro),bfct.append(b)
print bfct

/print "END-LOG"
