
from chempy import Storage

import cPickle

class PKL(Storage):

   def fromFile(self,fname,**params):
      fp = open(fname)
      result = cPickle.load(fp)
      fp.close()
      return result

#---------------------------------------------------------------------------
   def toFile(self,indexed,fname,**params):
      fp = open(fname,'w')
      result = cPickle.dump(indexed,fp,1)
      fp.close()
      
#---------------------------------------------------------------------------
   def fromString(self,st):   
      return cPickle.loads(st)

#---------------------------------------------------------------------------
   def toString(self,model):
      return cPickle.dumps(st)

