
from chempy import Storage

import pickle

class PKL(Storage):

   def fromFile(self,fname,**params):
      fp = open(fname)
      result = load(fp)
      fp.close()
      return result

#---------------------------------------------------------------------------
   def toFile(self,indexed,fname,**params):
      fp = open(fname,'w')
      result = dump(indexed,fp)
      fp.close()
      
#---------------------------------------------------------------------------
   def fromString(self,st):   
      return loads(st)

#---------------------------------------------------------------------------
   def toString(self,model):
      return dumps(st)

