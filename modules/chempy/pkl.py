
from chempy import Storage

import pickle
import cPickle

class PKL(Storage):

   def fromFile(self,fname,**params):
      fp = open(fname,'rb')
      result = cPickle.load(fp)
      fp.close()
      return result

#---------------------------------------------------------------------------
   def toFile(self,indexed,fname,**params):
      fp = open(fname,'wb')
      if(not params.has_key('bin')):
         result = cPickle.dump(indexed,fp,1)
      else:
         result = cPickle.dump(indexed,fp,params['bin'])         
      fp.close()

#---------------------------------------------------------------------------
   def fromStream(self,fp,**params):
      try:
         return cPickle.load(fp)
      except EOFError:
         return None

#---------------------------------------------------------------------------
   def toStream(self,indexed,fp,**params):
      if(not params.has_key('bin')):
         result = cPickle.dump(indexed,fp,1)
      else:
         result = cPickle.dump(indexed,fp,params['bin'])         
 
#---------------------------------------------------------------------------
   def fromString(self,st):
      return cPickle.loads(st)

#---------------------------------------------------------------------------
   def toString(self,model):
      return cPickle.dumps(st)

