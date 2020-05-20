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

from chempy import Storage,Atom,feedback
from chempy.models import Indexed,Connected
import string
import copy

generics =  {
    'CT' : 1,
    'C'  : 2,
    'CA' : 3,
    'CM' : 4,
    'CC' : 5,
    'CV' : 6,
    'CW' : 7,
    'CR' : 8,
    'CB' : 9,
    'C*' : 10,
    'CN' : 11,
    'CK' : 12,
    'CQ' : 13,
    'N'  : 14,
    'NA' : 15,
    'NB' : 16,
    'NC' : 17,
    'N*' : 18,
    'N2' : 19,
    'N3' : 20,
    'OW' : 21,
    'OH' : 22,
    'OS' : 23,
    'O'  : 24,
    'O2' : 25,
    'S'  : 26,
    'SH' : 27,
    'P'  : 28,
    'H'  : 29,
    'HW' : 30,
    'HO' : 31,
    'HS' : 32,
    'HA' : 33,
    'HC' : 34,
    'H1' : 35,
    'H2' : 36,
    'H3' : 37,
    'HP' : 38,
    'H4' : 39,
    'H5' : 40,
}

class XYZ(Storage):

#------------------------------------------------------------------------------
    def toList(self,model,mapping=None):

        conn = copy.deepcopy(model)
        conn = conn.convert_to_connected()

        list = []

        if len(model.atom):
            if not model.atom[0].has('numeric_type'):
                if not mapping:
                    mapping = generics
                    if feedback['warnings']:
                        print(' '+str(self.__class__)+': no numeric atom types found, using defaults.')
            list.append("%6d\n" % conn.nAtom)
            c = 0
            for a in conn.atom:
                if mapping:
                    n_type = mapping[a.text_type]
                else:
                    n_type = a.numeric_type
                if n_type<0:
                    print(str(self.__class__)+\
                            '-WARNING: negative numeric type (%d) for atom %d'% (n_type,c))
                st = "%6d  %-3s%12.6f%12.6f%12.6f%6d" % (
                    c+1,a.text_type,a.coord[0],a.coord[1],a.coord[2],int(n_type))
                for b in conn.bond[c]:
                    idx = b.index[0]
                    if idx == c:
                        idx = b.index[1]
                    st = st + "%6d" % (idx+1)
                st = st + "\n"
                list.append(st)
                c = c + 1
        return(list)

#----------------------------------------------------------------------------
    def updateFromList(self,model,list):

        c = 0
        for a in list[1:]:
            model.atom[c].coord = [ float(a[11:23]),float(a[23:35]),float(a[35:47])]
            c = c + 1
