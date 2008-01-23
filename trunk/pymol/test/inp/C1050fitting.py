# -c

import pymol
from pymol import cmd

from random import random,seed

print "BEGIN-LOG"

seed(123)

print "the 1st random number should be %8.3f\notherwise the rest of the test is meaningless...\n"%pymol.random()

pymol.random = random

cmd.load("dat/pept.pdb","ref")

for a in xrange(1,11):
   cmd.create("trg","ref",1,a,quiet=0)
   cmd.alter_state(a,"trg","x=x+random()/2")
   cmd.alter_state(a,"trg","y=y+random()/2")
   cmd.alter_state(a,"trg","z=z+random()/2",quiet=0)


cmd.frame(1)
print "%8.3f"%cmd.fit("ref","trg")

# asdf


for a in xrange(1,14):
   print a,
   print "%8.3f"%cmd.fit("ref and resi %d"%a,"trg"),
   print "%8.3f"%cmd.rms("ref","trg"),
   print "%8.3f"%cmd.rms_cur("ref","trg")

cmd.frame(10)

print "%8.3f"%cmd.fit("ref","trg")
for a in xrange(1,14):
   print a,
   print "%8.3f"%cmd.fit("ref and resi %d"%a,"trg"),
   print "%8.3f"%cmd.rms("ref","trg"),
   print "%8.3f"%cmd.rms_cur("ref","trg")


a = 1
print "%8.3f"%cmd.fit("ref","trg")
for b in xrange(1,11):
   cmd._dump_floats(cmd.intra_fit("trg and resi %d"%a,b))
   cmd._dump_floats(cmd.intra_rms("trg",b))
   cmd._dump_floats(cmd.intra_rms_cur("trg",b))

cmd.do("intra_fit (trg),10")
cmd.do("intra_fit (trg and resi 1),10")
cmd.do("intra_rms (trg),10")
cmd.do("intra_rms_cur (trg),10")

print "END-LOG"





