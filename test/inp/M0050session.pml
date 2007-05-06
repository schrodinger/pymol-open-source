# -c

/print "BEGIN-LOG"

from pymol2 import PyMOL

a = PyMOL()
print a != None
a.start()
a.cmd.load("dat/pept.pdb",quiet=0)
a.cmd.save("tmp/M0050session.pse",quiet=0)
a.cmd.reset()
print a.cmd.count_atoms()
a.cmd.load("tmp/M0050session.pse",quiet=0)
print a.cmd.count_atoms()

print cmd.count_atoms()
a.stop()
del a

/print "END-LOG"
