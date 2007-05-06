# -c

/print "BEGIN-LOG"

from pymol2 import PyMOL

lst = []

for i in range(0,3): lst.append( PyMOL() )

for pi in lst: print pi.__class__

for pi in lst: pi.start()

load dat/helix_amber.pdb
lst[0].cmd.do("load dat/pept.pdb")
lst[1].cmd.do("load dat/3al1.pdb")
lst[2].cmd.do("load dat/1tii.pdb")

count_atoms

for pi in lst: print pi.cmd.count_atoms()

for pi in lst: pi.stop()

del lst

/print "END-LOG"
