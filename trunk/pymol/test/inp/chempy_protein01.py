# python

from chempy import io
from chempy import protein


def dump(a):
   print 'atoms:',a.nAtom
   c = 0
   for at in a.atom:
      print at.name,
      c = c + 1
      if c == 15:
         print
         c = 0
   print
   print 'bonds:',a.nBond
   c = 0
   for b in a.bond:
      print b.index[0],b.index[1],b.order,
      c = c + 1
      if c == 6:
         print
         c = 0
   print

a = io.pdb.fromFile("dat/il2.pdb")
print 'before'
dump(a)

protein.add_bonds(a)
print 'add_bonds'
dump (a)

a = io.pdb.fromFile("dat/il2.pdb")

protein.assign_types(a)
print 'assign_types'
dump (a)


protein.add_bonds(a)
print 'add_bonds'
dump (a)

