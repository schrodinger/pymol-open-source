
from chempy import Storage

class LST(Storage):

   def fromFile(self,fname,**params):
      fp = open(fname)
      result = fp.readlines()
      fp.close()
      return result

#---------------------------------------------------------------------------
   def toFile(self,list,fname,**params):
      fp = open(fname,'w')
      fp.writelines(list)
      fp.close()
      

