from chempy import Storage,Atom
from chempy.models import Indexed,Connected
import string
import copy

class ARC(Storage):
   
#---------------------------------------------------------------------------------
   def fromFile(self,fname):
      list = []
      f = open(fname)
      while 1:
         hdr = f.readline()
         if not hdr: break
         hdr = hdr[0:6]
         if hdr:
            lst = []
            for b in xrange(int(hdr)):
               a = f.readline()
               lst.append([float(a[11:23]),float(a[23:35]),float(a[35:47])])
            list.append(lst)
      f.close()
      return list
   
   



