# pymol

from chempy import io
from chempy import protein
from pymol import cmd

from chempy import feedback

feedback['verbose']=1

a = io.pdb.fromFile("dat/helix_amber.pdb")
protein.add_bonds(a)
protein.assign_types(a)
c = 0
a.sort()
io.pdb.toFile(a,"cmp/pc_io01.01.pdb")
io.xyz.toFile(a,"cmp/pc_io01.01.xyz")
for at in a.atom:
   print "%5s%5s%5s %6.3f%3s\n" % (
      at.name,at.resi,at.resn,at.partial_charge,at.text_type)
   c = c + 1
   if c == 100:
      break;
cmd.load_model(a,"test")
b = cmd.get_model()
c = 0
b.sort()
for at in b.atom:
   print "%5s%5s%5s %6.3f%3s\n" % (
      at.name,at.resi,at.resn,at.partial_charge,at.text_type)
   c = c + 1
   if c == 100:
      break;
io.pdb.toFile(b,"cmp/pc_io01.02.pdb")
io.xyz.toFile(b,"cmp/pc_io01.02.xyz")
cmd.quit()
