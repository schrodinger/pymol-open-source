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

from chempy import Storage,Atom
from chempy.models import Indexed
import string

class PDB(Storage):
    
#---------------------------------------------------------------------------------
    def fromList(self,list):   # currently no handling of conect records 

        model = Indexed()
        
# read atoms
        cnt = 0
        at = None
        for rec in list:
            if (rec[0:4]=='ATOM') or (rec[0:6]=='HETATM'):
                at = Atom()
                if rec[0]=='A': at.hetatm=0  # default is 1
                at.index = cnt
                at.name = string.strip(rec[12:16])
                at.alt = string.strip(rec[16:17])
                at.resn = string.strip(rec[17:20])
                at.chain = string.strip(rec[21:22])
                at.resi = string.strip(rec[22:27]) # note: insertion is part of resi
                at.resi_number = int(rec[22:26])
                at.coord = [float(rec[30:38]), 
                                float(rec[38:46]),
                                float(rec[46:54])]
                try:
                    at.q = float(rec[54:60])
                except ValueError:
                    at.q = 1.0
                try:               
                    at.b = float(rec[60:66])
                except ValueError:
                    at.b = 0.0
                at.segi = string.strip(rec[72:76])
                at.symbol = string.strip(rec[76:78])
                if not len(at.symbol):
                    at.symbol = at.name[0:1]
                    if at.symbol in '012345678':
                        at.symbol = at.name[1:2]
                cnt = cnt + 1
                model.add_atom(at)
            elif (rec[0:3]=='TER'):
                if at:
                    at.ter=1
        return(model)

#---------------------------------------------------------------------------------
    def toList(self,model):

        list = []

        cnt = 1
        for at in model.atom:
            if at.hetatm:
                het = 'HETATM'
            else:
                het = 'ATOM  '
            latn = len(at.name)
            if latn and (latn<4):
                if at.name[0] in '0123456789':
                    name=at.name
                else:
                    name=' '+at.name
            else:
                name=at.name
            lrsi = len(at.resi)
            if lrsi and (lrsi<5):
                if at.resi[lrsi-1] in '0123456789':
                    resi=at.resi+' '
                else:
                    resi=at.resi
            list.append(
                "%6s%5i %-4s%1s%3s %1s%5s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s\n" % 
                (het,cnt,name,at.alt,at.resn,at.chain,resi,
                at.coord[0],at.coord[1],at.coord[2],at.q,at.b,at.segi,at.symbol))
            if hasattr(at,'ter'):
                if at.ter:
                    cnt = cnt + 1
                    list.append(
                        "%6s%5i      %3s %1s%4s\n" %
                        ('TER   ',cnt,at.resn,at.chain,at.resi))
            cnt = cnt + 1
        list.append("END\n")
        return(list)

