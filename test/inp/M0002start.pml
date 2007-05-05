# -c

/print "BEGIN-LOG"

from pymol2 import PyMOL

a = PyMOL()
print a != None
a.start()
print a.cmd.count_atoms()
a.stop()
del a

/print "END-LOG"
