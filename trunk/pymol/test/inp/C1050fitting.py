# -c

from whrandom import random,seed

print "BEGIN-LOG"

print random()

seed(123)

cmd.load("dat/pept.pdb","ref")

for a in xrange(1,11):
   cmd.create("trg","ref",1,a)
   cmd.alter_state(a,"trg","x=x+random()/3")
   cmd.alter_state(a,"trg","y=y+random()/3")
   cmd.alter_state(a,"trg","z=z+random()/3")

util.cbay("ref")

cmd.button("l","shft","pkat")
cmd.button("r","ctrl","move")

cmd.wizard("pair_fit")

print "END-LOG"





