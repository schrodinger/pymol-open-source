# pymol -c

from pymol import cmd

try:
   from chempy import io
   from chempy.fast import FastModel

   cmd.load("dat/pept.pdb")
   cmd.flag("2","(elem c)")
   m1 = cmd.get_model()
   
   f1 = FastModel()
   f1.from_indexed(m1)
   m2 = f1.convert_to_indexed()

   cmd.load_model(m2,'conv')
   c = 0
   for a in m2.atom:
      c = c + 1
      print a.name,a.resi,a.segi,a.coord[0],a.flags
      if c == 20:
         break

   io.pdb.toFile(m2,"cmp/fast02.01.pdb")
   
except ImportError:
   print ' error: this test requires the Numerical Python package'


