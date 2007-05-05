# -c

/print "BEGIN-LOG"

from pymol2 import PyMOL

lst = []

for i in range(0,3): lst.append( PyMOL() )

for pi in lst: print pi.__class__

for pi in lst: pi.start()

for pi in lst: pi.stop()

del lst

/print "END-LOG"
