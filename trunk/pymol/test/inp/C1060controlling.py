# -c

from pymol import cmd

print "BEGIN-LOG"

print cmd.set_key('F1',lambda :cmd.turn('x',10))
print cmd.set_key('F12',lambda :cmd.turn('x',10))

print cmd.set_key('left',lambda :cmd.turn('x',10))
print cmd.set_key('right',lambda :cmd.turn('x',10))
print cmd.set_key('pgup',lambda :cmd.turn('x',10))
print cmd.set_key('pgdn',lambda :cmd.turn('x',10))
print cmd.set_key('home',lambda :cmd.turn('x',10))
print cmd.set_key('insert',lambda :cmd.turn('x',10))

print cmd.set_key('ALT-A',lambda :cmd.turn('y',10))

print cmd.set_key('CTRL-C',lambda :cmd.turn('z',10))

print cmd.set_key('SHFT-F1', lambda :cmd.turn('z',10))

print cmd.set_key('ALT-F1', lambda :cmd.turn('y',10))

print cmd.set_key('CTRL-F8', lambda :cmd.move('x',1))


print "END-LOG"





