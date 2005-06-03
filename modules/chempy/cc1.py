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

from chempy.models import Indexed
from chempy import Storage,Atom,Bond

import string

class CC1(Storage): # ChemDraw3D 5.0 std., cartesian coordinates

    def fromList(self,molList):

        model = Indexed()

        # read header information
        nAtom = int(molList[0][0:3])

        # read atoms and bonds
        id_dict = {}
        irec = 1
        cnt = 0
        for a in range(nAtom):
            at = Atom()
            at.index = cnt
            id_dict[string.strip(molList[irec][3:8])] = at.index
            at.coord = [float(molList[irec][8:20]), 
                float(molList[irec][20:32]),float(molList[irec][32:44])]
            at.symbol = string.strip(molList[irec][0:3])
            at.numeric_type = int(molList[irec][44:49])
            lst = string.split(string.strip(molList[irec][49:]))
            at.bonds = lst
            irec = irec + 1
            cnt = cnt + 1
            model.atom.append(at)
            
        # interpret bonds
        cnt = 0
        for a in model.atom:
            lst = a.bonds
            del a.bonds
            for b in lst:
                if a.index<id_dict[b]:
                    bnd = Bond()
                    bnd.index = [ a.index,id_dict[b]]
                    model.bond.append(bnd)

        # obtain formal charges from M  CHG record
        return model

#------------------------------------------------------------------------------
    def toList(self,model):

        return []

        # not implemented yet

