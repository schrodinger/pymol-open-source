# python

try:
   from chempy import io
   from chempy.fast import FastModel

   m1 = io.mol.fromFile("dat/small01.mol")
   f1 = FastModel()
   f1.from_indexed(m1)
   m2 = f1.convert_to_indexed()

   io.mol.toFile(m2,"cmp/fast01.01.mol")
   
except ImportError:
   print ' error: this test requires the Numerical Python package'


