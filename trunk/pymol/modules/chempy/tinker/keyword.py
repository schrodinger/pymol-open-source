
import chempy

def get_partial_charge(model):
   if chempy.feedback['verbose']:
      print ' '+str(__name__)+': generating partial charge keywords...'
   list = []
   c = -1
   for a in model.atom:
      list.append("CHARGE %d %6.4f\n" %(c,a.partial_charge))
      c = c - 1
   return list



