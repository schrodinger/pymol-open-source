
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


def get_restrain_positions(model,flag,f_cnst):
   list = []
   n = 0
   c = 1
   for a in model.atom:
      if (a.flags&flag):
         list.append("RESTRAIN-POSITION %5d %8.3f %8.3f %8.3f %8.3f\n" %
                     (c,a.coord[0],a.coord[1],a.coord[2],f_cnst))
         n = n + 1
      c = c + 1
   if chempy.feedback['actions']:
      print ' '+str(__name__)+': %d atoms restrained...' % n
      
   return list

def get_inactive(model,flag):
   list = []
   n = 0
   c = 1
   for a in model.atom:
      if (a.flags&flag):
         list.append("INACTIVE %d\n" % (c))
         n = n + 1
      c = c + 1
   if chempy.feedback['actions']:
      print ' '+str(__name__)+': %d atoms fixed...' % n
   return list

