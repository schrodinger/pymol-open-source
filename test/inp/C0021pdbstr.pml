# -c

/print "BEGIN-LOG"

f=open("dat/pept.pdb","rb")
l=f.read()
f.close()

cmd.read_pdbstr(l,"test")

print cmd.count_atoms()

/print "END-LOG"
