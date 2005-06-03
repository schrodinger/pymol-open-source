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

class LST(Storage):

    def fromFile(self,fname,**params):
        fp = open(fname)
        result = fp.readlines()
        fp.close()
        return result

#---------------------------------------------------------------------------
    def toFile(self,list,fname,**params):
        fp = open(fname,'w')
        try:
            fp.writelines(list)
        except TypeError:
            for a in list:
                fp.write(str(a)+"\n")
        fp.close()
        

