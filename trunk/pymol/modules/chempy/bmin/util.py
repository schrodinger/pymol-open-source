# pymol

from chempy import Bond
import copy
import whrandom

# cap gaps in chain as aldehydes and neutral amines...
# (for instance, when

def cap(object):
   from pymol import cmd
   
   model = cmd.get_model(object)
   # guarantee identical ordering
   cmd.delete(object)
   cmd.load_model(model,object)
   n_list = cmd.identify("(n;n &!(n;c a;2.0))")
   c_list = cmd.identify("(n;c &!(n;n a;2.0))")
   print n_list
   print c_list
   for a in n_list:
      newat = copy.deepcopy(model.atom[a])
      newat.coord = [
         newat.coord[0] + whrandom.random(),
         newat.coord[1] + whrandom.random(),
         newat.coord[2] + whrandom.random(),
         ]
      newat.symbol = 'H'
      newat.name = 'HN'
      newat.numeric_type = 43
      bond = Bond()
      bond.order = 1
      bond.stereo = 0
      bond.index = [ a, model.nAtom ]
      print "adding",newat.name,bond.index
      model.add_atom(newat)
      model.add_bond(bond)
   for a in c_list:
      newat = copy.deepcopy(model.atom[a])
      newat.coord = [
         newat.coord[0] + whrandom.random(),
         newat.coord[1] + whrandom.random(),
         newat.coord[2] + whrandom.random(),
         ]
      newat.symbol = 'H'
      newat.name = 'HC'
      newat.numeric_type = 41
      bond = Bond()
      bond.order = 1
      bond.stereo = 0
      bond.index = [ a, model.nAtom ]
      print "adding",newat.name,bond.index
      model.add_atom(newat)
      model.add_bond(bond)
   # reload
   cmd.delete(object)
   cmd.load_model(model,object)
   cmd.sort(object)
   

