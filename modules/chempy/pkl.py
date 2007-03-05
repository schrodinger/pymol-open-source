#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information. 
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-* 
#-* 
#-*
#Z* -------------------------------------------------------------------

from chempy import Storage

import pickle
import cPickle

class PKL(Storage):

    def fromFile(self,fname,**params):
        fp = self.my_open(fname,'rb')
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
        return cPickle.dumps(model)

