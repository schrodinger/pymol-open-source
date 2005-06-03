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

from chempy.models import Indexed,Connected
from chempy import Storage,Atom,Bond

import string
import copy

class MMD(Storage):
    
#---------------------------------------------------------------------------------
    def fromList(self,MMODList):

        model = Connected()
        
# get header information
        nAtom = int(MMODList[0][1:6])
        model.molecule.title = string.strip(MMODList[0][8:])
        irec = 1

# loop through atoms
        cnt = 0
        for a in range(nAtom):
            model.bond.append([])
            
        for a in range(nAtom):
            at = Atom()
            at.numeric_type = int(MMODList[irec][1:4])

# extract connectivity information
            tokens = string.splitfields(MMODList[irec][5:52])
            at.neighbor = []
            at.bondorder = []

            for i in range(6):
                if tokens[2*i] != "0":
                    a2 = int(tokens[2*i])-1
                    if (a2>cnt):
                        b = Bond()
                        b.index = [cnt,a2]
                        b.order = int(tokens[2*i+1])
                        model.bond[b.index[0]].append(b) # note two refs to same object
                        model.bond[b.index[1]].append(b) # note two refs to same object 
                else:
                    break
            
# extract other information
            at.coord = [float(MMODList[irec][53:64]), 
                float(MMODList[irec][65:76]), float(MMODList[irec][77:88])]
            at.resi = string.strip(MMODList[irec][89:94])
            at.resi_number = int(at.resi)
            resn_code = string.strip(MMODList[irec][94:95])
            if len(resn_code): at.resn_code = resn_code
            color_code = string.strip(MMODList[irec][96:100])
            if color_code!='':
                at.color_code = int(color_code)
            else:
                at.color_code = 0
            chain = string.strip(MMODList[irec][95:96])
            if len(chain): at.chain = chain
            at.partial_charge = float(MMODList[irec][100:109])
            at.resn = MMODList[irec][119:123]
            name = string.strip(MMODList[irec][124:128])
            if len(name): at.name = name
            model.atom.append(at)
            irec = irec + 1
            cnt = cnt + 1
            
# fill in remaining datatypes
        cnt = 1
        for a in model.atom:
            a.text_type = MMOD_atom_data[a.numeric_type][0]
            a.symbol = MMOD_atom_data[a.numeric_type][1]
            a.formal_charge = MMOD_atom_data[a.numeric_type][4]
            cnt = cnt + 1
            
        return(model.convert_to_indexed())

#---------------------------------------------------------------------------------
    def updateFromList(self,model,list): # updates charges and coordinates
        nAtom = int(list[0][1:6])
        try:
            model.molecule.energy = float(list[0][58:68])/4.184 # convert to kcal
        except:
            if hasattr(model.molecule,'energy'):
                del model.molecule.energy
        if nAtom!=model.nAtom:
            raise RuntimeError(" mmd: atom count mismatch")
        c = 0
        for a in list[1:]:
            mac = model.atom[c]
            mac.coord = [float(a[53:64]), 
                float(a[65:76]), float(a[77:88])]
            mac.partial_charge = float(a[100:109])         
            c = c + 1
            
#---------------------------------------------------------------------------------
    def toList(self,model,no_blank_names=1):

        conn = copy.deepcopy(model)
        conn = conn.convert_to_connected()
        
        MMODList = []
        MMODList.append(" %5d  %-70s\n" %(conn.nAtom,conn.molecule.title))
        c = 0
        for i in conn.atom:
            
# construct neighbor list
            neighbors = 6*[0]
            bondorders = 6*[0]
            j = 0
            for b in conn.bond[c]:
                n = b.index[0]
                if n == c:
                    n = b.index[1]
                neighbors[j] = n + 1
                bondorders[j] = b.order
                j = j + 1
                
# assemble output line
            if i.numeric_type>0:
                tline = " %3d" % (i.numeric_type)
            else:
                tline = " %3d" % 64
            for j in range(6): 
                tline = tline + " %5d %1d" % (neighbors[j], bondorders[j])
            tline = tline + " %11.6f %11.6f %11.6f " % (i.coord[0], 
                i.coord[1], i.coord[2])
            name = i.name
            if not len(name):
                if no_blank_names:
                    name = "%s%d" % (i.symbol,c+1)
                    name = name[0:4]
            tline = tline + "%5d%1s%1s%4d%9.5f%9.5f %-4s %4s\n" % \
                (i.resi_number,i.resn_code, i.chain, i.color_code,
                 i.partial_charge, i.partial_charge, i.resn, name)
            MMODList.append(tline)
            c = c + 1
            
        return(MMODList)


#---------------------------------------------------------------------------------
'#Ntype Atype Elem Hybr Att Chg\n',
MMOD_atom_data = {
    1: ['C1','C' ,'sp' , 2, 0],
    2: ['C2','C' ,'sp2', 3, 0],
    3: ['C3','C' ,'sp3', 4, 0],
    4: ['CA','C' ,'sp3', 3, 0],
    5: ['CB','C' ,'sp3', 2, 0],
    6: ['CC','C' ,'sp3', 1, 0],
    7: ['CD','C' ,'sp2', 2, 0],
    8: ['CE','C' ,'sp2', 1, 0],
    9: ['CF','C' ,'sp' , 1, 0],
  10: ['CM','C' ,'unk',-1,-1],
  11: ['CP','C' ,'unk',-1, 1],
  12: ['CR','C' ,'unk',-1, 0],
  14: ['C0','C' ,'unk',-1, 0],
  15: ['O2','O' ,'sp2', 1, 0],
  16: ['O3','O' ,'sp3', 2, 0],
  17: ['OA','O' ,'sp3', 1, 0],
  18: ['OM','O' ,'sp3', 1,-1],
  19: ['OW','O' ,'sp3', 0, 0],
  20: ['OP','O' ,'sp2', 2, 1],
  21: ['OQ','O' ,'sp3', 3, 1],
  23: ['O0','O' ,'unk',-1, 0],
  24: ['N1','N' ,'sp' , 1, 0],
  25: ['N2','N' ,'sp2', 2, 0],
  26: ['N3','N' ,'sp3', 3, 0],
  27: ['NA','N' ,'sp3', 2, 0],
  28: ['NB','N' ,'sp3', 1, 0],
  29: ['NC','N' ,'sp3', 0, 0],
  30: ['ND','N' ,'sp2', 1, 0],
  31: ['N4','N' ,'sp2', 3, 1],
  32: ['N5','N' ,'sp3', 4, 1],
  33: ['NE','N' ,'sp3', 3, 1],
  34: ['NF','N' ,'sp3', 2, 1],
  35: ['NG','N' ,'sp3', 1, 1],
  36: ['NH','N' ,'sp2', 2, 1],
  37: ['NI','N' ,'sp2', 1, 1],
  40: ['N0','N' ,'unk',-1, 0],
  41: ['H1','H' ,'s'  , 1, 0],
  42: ['H2','H' ,'s'  , 1, 0],
  43: ['H3','H' ,'s'  , 1, 0],
  44: ['H4','H' ,'s'  , 0, 0],
  45: ['H5','H' ,'s'  , 0, 0],
  48: ['H0','H' ,'s'  ,-1, 0],
  49: ['S1','S' ,'sp3', 2, 0],
  50: ['SA','S' ,'sp3', 1, 0],
  51: ['SM','S' ,'sp3', 0,-1],
  52: ['S0','S' ,'unk',-1, 0],
  53: ['P0','P' ,'unk',-1, 0],
  54: ['B2','B' ,'sp2', 2, 0],
  55: ['B3','B' ,'sp3', 3, 0],
  56: ['F0','F' ,'sp3', 1, 0],
  57: ['Cl','Cl','sp3', 1, 0],
  58: ['Br','Br','sp3', 1, 0],
  59: ['I0','I' ,'sp3', 1, 0],
  60: ['Si','Si','unk',-1, 0],
  61: ['Du','Du','unk',-1, 0],
  62: ['Du','Du','unk',-1, 0],
  63: ['Lp','Lp','unk', 1, 0],
  64: ['Du','Du','unk',-1, 0]};

