#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information. 
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-* Scott Dixon, Metaphorics, LLC
#-* 
#-*
#Z* -------------------------------------------------------------------

from chempy import feedback

import string
import copy

default_extra = { # atomic number, and normal valency (for tinker)
    # NOTE: THIS SET IS ONLY USED IF THERE ARE NO TINKER RECORDS
    # IN THE PARAMETER FILE
    # bromine 
    'BR' : [                 35 ,               1 ],
    # carbon                           
    'C'  : [                  6 ,               3 ],
    'CA' : [                  6 ,               3 ],
    'CB' : [                  6 ,               3 ],
    'CC' : [                  6 ,               3 ],
    'CF' : [                  6 ,               3 ],
    'CK' : [                  6 ,               3 ],
    'CM' : [                  6 ,               3 ],
    'CN' : [                  6 ,               3 ],
    'CQ' : [                  6 ,               3 ],
    'CR' : [                  6 ,               3 ],
    'CT' : [                  6 ,               4 ],
    'CV' : [                  6 ,               3 ],
    'CW' : [                  6 ,               3 ], 
    'CY' : [                  6 ,               2 ],
    'CX' : [                  6 ,               2 ],
    'C5' : [                  6 ,               3 ],
    'C*' : [                  6 ,               3 ],
    # calcium
    'C0' : [                 20 ,               3 ],
    # chloride
    'Cl' : [                 17 ,               1 ], 
    # fluorine
    'F'  : [                  9 ,               1 ],
    # hydrogen
    'H'  : [                  1 ,               1 ],
    'H1' : [                  1 ,               1 ],
    'H2' : [                  1 ,               1 ],
    'H3' : [                  1 ,               1 ],
    'H4' : [                  1 ,               1 ],
    'H5' : [                  1 ,               1 ],
    'HA' : [                  1 ,               1 ],
    'HC' : [                  1 ,               1 ],
    'HO' : [                  1 ,               1 ],
    'HP' : [                  1 ,               1 ],
    'HS' : [                  1 ,               1 ],
    'HW' : [                  1 ,               2 ],
    # iodine
    'I'  : [                 53 ,               1 ],
    # chloride anion
    'IM' : [                 17 ,               0 ],
    # sodium cation
    'IP' : [                 11 ,               0 ],
    # magnesium
    'MG' : [                 12 ,               0 ],
    # nitrogen
    'N'  : [                  7 ,               3 ],
    'N*' : [                  7 ,               2 ],
    'N2' : [                  7 ,               3 ],
    'N3' : [                  7 ,               4 ],
    'NA' : [                  7 ,               3 ],
    'NB' : [                  7 ,               2 ],
    'NC' : [                  7 ,               2 ],
    'NX' : [                  7 ,               2 ],
    'NZ' : [                  7 ,               3 ],
    #
    'OW' : [                  8 ,               2 ],
    'OH' : [                  8 ,               2 ],
    'OS' : [                  8 ,               2 ],
    'O'  : [                  8 ,               1 ],
    'O2' : [                  8 ,               1 ],
    'OM' : [                  8 ,               1 ],
    'OZ' : [                  8 ,               1 ],   
    'P'  : [                 15 ,               4 ],
    'S'  : [                 16 ,               2 ],
    'SO'  : [                16 ,               4 ],
    'SX'  : [                16 ,               2 ],
    'SH' : [                 16 ,               2 ],
    # copper
    'CU' : [                 29 ,               0 ],
    # iron
    'FE' : [                 26 ,               0 ],
    # lithium
    'Li' : [                  3 ,               0 ],
    # potassium
    'K'  : [                 19 ,               0 ],
    # rubidium
    'Rb' : [                 37 ,               0 ],
    # cesium
    'Cs' : [                 55 ,               0 ],
    # calcium (alternate)
    'Ca' : [                 20 ,               0 ],
    # dummy
    'X ' : [                  6 ,               1 ],
    # zinc
    'Z0' : [                 30 ,               0 ],
    'Z4' : [                 30 ,               4 ],
    'Z5' : [                 30 ,               5 ],
    }

class Parameters:

    def __init__(self,fname):
        if feedback['actions']:
            print ' '+str(self.__class__)+': loading from "%s"...' % fname
        f = open(fname)
        # skip
        l = f.readline()
        # read names & molecular weights
        self.type = []
        self.mw = {}
        while 1:
            l = string.strip(f.readline())
            if not len(l): break
            a2 = string.strip(l[0:2])
            self.type.append(a2)
            self.mw[a2] = [float(l[3:13]),string.strip(l[34:])]
        # skip 1
        l = f.readline()
        # read bonds
        self.bond = {}
        while 1:
            l = string.strip(f.readline())
            if not len(l): break
            a5 = l[0:5]
            self.bond[a5] = [float(l[5:12]),float(l[12:22]),
                                  string.strip(l[22:])]
        # read angles
        self.angle = {}
        while 1:
            l = string.strip(f.readline())
            if not len(l): break
            a5 = l[0:8]
            self.angle[a5] = [float(l[8:16]),float(l[16:28]),
                                  string.strip(l[28:])]
        # read torsion 
        self.torsion = {}
        while 1:
            l = string.strip(f.readline())
            if not len(l): break
            a5 = l[0:11]
            if self.torsion.has_key(a5):
                self.torsion[a5].extend(
                    [int(l[11:15]),
                     float(l[15:27]),
                     float(l[27:36]),
                     abs(int(float(l[40:52]))),
                     string.strip(l[52:])])
            else:
                self.torsion[a5] = [int(l[11:15]),
                                          float(l[15:27]),
                                          float(l[27:36]),
                                          abs(int(float(l[40:52]))),
                                          string.strip(l[52:])]
        # read impropers
        self.improper = {}
        while 1:
            l = string.strip(f.readline())
            if not len(l): break
            a5 = l[0:11]
            self.improper[a5] = [float(l[15:27]),
                                        float(l[27:40]),
                                        abs(int(float(l[40:51]))),
                                        string.strip(l[51:])]
        # skip
        while 1:
            l = string.strip(f.readline())
            if not len(l): break
        # read vdw equivalents 
        self.vdw_eq = {}
        while 1:
            l = string.strip(f.readline())
            if not len(l): break
            a4 = string.strip(l[0:4])
            l = l[4:]
            while len(l):
                self.vdw_eq[string.strip(l[0:4])] = a4
                l = l[4:]
        # skip
        l = string.strip(f.readline())
        # read vdw parameters
        self.vdw = {}
        while 1:
            l = string.strip(f.readline())
            if not len(l): break
            l = '  ' + l
            a4 = string.strip(l[0:4])
            self.vdw[a4] =  [float(l[4:20]),
                                  float(l[20:37]),
                                  string.strip(l[37:])]
            
        # read extra tinker information if present
        self.extra = {}
        while 1:
            l = f.readline()
            if not l: break
            if l[0:6] == 'TINKER':
                self.extra[string.strip(l[6:12])]  = [
                    int(l[12:18]),
                    int(l[18:24])]
        if not len(self.extra.keys()):
                self.extra = default_extra
        # now generate redundant equivalents
        for a in self.vdw_eq.keys():
            self.vdw[a] = self.vdw[self.vdw_eq[a]]
        for a in self.bond.keys():
            k = a[3:5]+'-'+a[0:2]
            self.bond[k] = self.bond[a]
        for a in self.angle.keys():
            k = a[6:8]+'-'+a[3:5]+'-'+a[0:2]
            self.angle[k] = self.angle[a]
            
    def dump(self):
        for b in self.type:
            print b
        kees = self.mw.keys()
        kees.sort()
        for b in kees:
            print b,self.mw[b]
        kees = self.vdw_eq.keys()
        kees.sort()
        for b in kees:
            print b,self.vdw_eq[b]
        kees = self.vdw.keys()
        kees.sort()
        for b in kees:
            print b,self.vdw[b]
        kees = self.bond.keys()
        kees.sort()
        for b in kees:
            print b,self.bond[b]
        kees = self.angle.keys()
        kees.sort()
        for b in kees:
            print b,self.angle[b]
        kees = self.torsion.keys()
        kees.sort()
        for b in kees:
            print b,self.torsion[b]
        kees = self.improper.keys()
        kees.sort()
        for b in kees:
            print b,self.improper[b]
        kees = self.vdw.keys()
        kees.sort()
        for b in kees:
            print b,self.vdw[b]
            

class Topology:

    def __init__(self,model):
        self.model = model
        if feedback['actions']:
            print ' '+str(self.__class__)+': searching...'
        # get a connected version too
        cmodel = copy.deepcopy(model).convert_to_connected()
        # find atom types in molecule
        self.present = {}
        for a in model.atom:
            self.present[a.text_type] = 1
        # copy bonds
        self.bond = {}
        for b in model.bond:
            a0 = b.index[0]
            a1 = b.index[1]
            if a0 < a1:
                self.bond[(a0,a1)] = 1
            else:
                self.bond[(a1,a0)] = 1
        # find angles
        self.angle = {}
        ang = self.angle
        for b in model.bond:
            a0 = b.index[0]
            a1 = b.index[1]
            for c in cmodel.bond[a0]: # a0 in center
                a2 = c.index[0] 
                if a2 not in (a0,a1): # outside atom
                    if a1 < a2:
                        an = (a1,a0,a2)
                    else:
                        an = (a2,a0,a1)                  
                    ang[an] = 1
                a2 = c.index[1]
                if a2 not in (a0,a1): # outside atom
                    if a1 < a2:
                        an = (a1,a0,a2)
                    else:
                        an = (a2,a0,a1)                  
                    ang[an] = 1
            for c in cmodel.bond[a1]: # a1 in center
                a2 = c.index[0] 
                if a2 not in (a0,a1): # outside atom
                    if a0 < a2:
                        an = (a0,a1,a2)
                    else:
                        an = (a2,a1,a0)                  
                    ang[an] = 1
                a2 = c.index[1]
                if a2 not in (a0,a1): # outside atom
                    if a0 < a2:
                        an = (a0,a1,a2)
                    else:
                        an = (a2,a1,a0)                  
                    ang[an] = 1
        # find torsions
        self.torsion = {}
        tors = self.torsion
        for b in model.bond: # use bond as center of torsion
            a1 = b.index[0]
            a2 = b.index[1]
            for c in cmodel.bond[a1]: 
                a0 = c.index[0] 
                if a0 not in (a1,a2): # outside atom
                    for d in cmodel.bond[a2]:
                        a3 = d.index[0] 
                        if a3 not in (a0,a1,a2): # outside atom
                            if a0 < a3:
                                to = (a0,a1,a2,a3)
                            else:
                                to = (a3,a2,a1,a0)                        
                            tors[to] = 1
                        a3 = d.index[1] 
                        if a3 not in (a0,a1,a2): # outside atom
                            if a0 < a3:
                                to = (a0,a1,a2,a3)
                            else:
                                to = (a3,a2,a1,a0)
                            tors[to] = 1
                a0 = c.index[1] 
                if a0 not in (a1,a2): # outside atom
                    for d in cmodel.bond[a2]:
                        a3 = d.index[0] 
                        if a3 not in (a0,a1,a2): # outside atom
                            if a0 < a3:
                                to = (a0,a1,a2,a3)
                            else:
                                to = (a3,a2,a1,a0)                        
                            tors[to] = 1
                        a3 = d.index[1] 
                        if a3 not in (a0,a1,a2): # outside atom
                            if a0 < a3:
                                to = (a0,a1,a2,a3)
                            else:
                                to = (a3,a2,a1,a0)                        
                            tors[to] = 1
        # find impropers (only autogenerates for atoms with 3 bonds)
        self.improper = {}
        impr = self.improper
        a2 = 0
        for a in model.atom: # a2 is the center atom 
            bnd = cmodel.bond[a2]
            if len(bnd) == 3:
                lst = []
                for b in bnd:
                    at = b.index[0]
                    if at == a2:
                        at = b.index[1]
                    lst.append(at)
                lst.sort()
                impr[(lst[0],lst[1],a2,lst[2])] = 1
            a2 = a2 + 1
            
        if feedback['actions']:
            print ' '+str(self.__class__)+': found:'
            print ' '+str(self.__class__)+':    types       %6d' % (
                len(self.present.keys()))
            print ' '+str(self.__class__)+':    bonds       %6d' % (
                len(self.bond.keys()))
            print ' '+str(self.__class__)+':    angles      %6d' % (
                len(self.angle.keys()))
            print ' '+str(self.__class__)+':    torsions    %6d' % (
                len(self.torsion.keys()))
            print ' '+str(self.__class__)+':    impropers   %6d' % (
                len(self.improper.keys()))

    def dump(self):
        kees = self.present.keys()
        kees.sort()
        for b in kees:
            print b
        kees = self.bond.keys()
        kees.sort()
        for b in kees:
            print b
        kees = self.angle.keys()
        kees.sort()
        for b in kees:
            print b
        kees = self.torsion.keys()
        kees.sort()
        for b in kees:
            print b
        kees = self.improper.keys()
        kees.sort()
        for b in kees:
            print b

    def get_list(self):
        return [self.bond.keys(),
                  self.angle.keys(),
                  self.torsion.keys(),
                  self.improper.keys()]

class Subset:
    
    def __init__(self,par,top):
        if feedback['actions']:
            print ' '+str(self.__class__)+': applying parameter set to topology...'
        self.model = top.model
        self.present = copy.deepcopy(top.present)
        atype = []
        for a in self.model.atom:
            atype.append(a.text_type)
        # extra tinker info (atomic number, valency)
        self.miss_extra = []
        self.extra = {}
        s_extra = self.extra
        present = top.present
        p_extra = par.extra
        for kee in present.keys():
            if p_extra.has_key(kee):
                s_extra[kee] = p_extra[kee]
            else:
                self.miss_extra.append(kee)
        # molecular weight
        self.miss_mw = []
        self.mw = {}
        s_mw = self.mw
        present = top.present
        p_mw = par.mw
        for kee in present.keys():
            if p_mw.has_key(kee):
                s_mw[kee] = p_mw[kee]
            else:
                self.miss_mw.append(kee)
                
        # van der waals
        self.miss_vdw = []
        self.vdw = {}
        s_vdw = self.vdw
        present = top.present
        p_vdw = par.vdw
        for kee in present.keys():
            if p_vdw.has_key(kee):
                s_vdw[kee] = p_vdw[kee]
            else:
                self.miss_vdw.append(kee)
        # bonds
        self.miss_bond = []
        self.bond = {}
        s_bond = self.bond
        bond = top.bond
        p_bond = par.bond
        for a in bond.keys():
            kee = "%-2s-%-2s" % (atype[a[0]],atype[a[1]])
            if p_bond.has_key(kee):
                s_bond[kee] = p_bond[kee]
            else:
                self.miss_bond.append((a,kee))
        # angles
        self.miss_angle = []
        self.angle = {}
        s_angle = self.angle
        angle = top.angle
        p_angle = par.angle
        for a in angle.keys():
            kee = "%-2s-%-2s-%-2s" % (atype[a[0]],atype[a[1]],atype[a[2]])
            if p_angle.has_key(kee):
                s_angle[kee] = p_angle[kee]
            else:
                self.miss_angle.append((a,kee))
        # torsions      
        self.miss_torsion = []
        self.torsion = {}
        s_torsion = self.torsion
        torsion = top.torsion
        p_torsion = par.torsion
        for a in torsion.keys():
            at0=atype[a[0]]
            at1=atype[a[1]]
            at2=atype[a[2]] # center atom
            at3=atype[a[3]]
            while 1: # not a real loop, just "else" avoidance
                kee1 = "%-2s-%-2s-%-2s-%-2s" % (at0,at1,at2,at3)
                if p_torsion.has_key(kee1):
                    s_torsion[kee1] = p_torsion[kee1]
                    break
                kee2 = "%-2s-%-2s-%-2s-%-2s" % (at3,at2,at1,at0)
                if p_torsion.has_key(kee2):
                    s_torsion[kee2] = p_torsion[kee2]
                    break
                kee = "X -%-2s-%-2s-X " % (at1,at2)
                if p_torsion.has_key(kee):
                    s_torsion[kee1] = p_torsion[kee]
                    break
                kee = "X -%-2s-%-2s-X " % (at2,at1)
                if p_torsion.has_key(kee): 
                    s_torsion[kee2] = p_torsion[kee]
                    break
                self.miss_torsion.append((a,kee1))
                break
        # impropers      
        self.miss_improper = []
        self.improper = {}
        s_improper = self.improper
        improper = top.improper
        p_improper = par.improper
        for a in improper.keys():
            at0=atype[a[0]]
            at1=atype[a[1]]
            at2=atype[a[2]] # center atom
            at3=atype[a[3]]
            while 1: # not a real loop, just "else" avoidance

                kee1 = "%-2s-%-2s-%-2s-%-2s" % (at1,at3,at2,at0)
                if p_improper.has_key(kee1):
                    s_improper[kee1] = p_improper[kee1]
                    break
                kee2 = "%-2s-%-2s-%-2s-%-2s" % (at3,at1,at2,at0)
                if p_improper.has_key(kee2):
                    s_improper[kee2] = p_improper[kee2]
                    break

                kee3 = "%-2s-%-2s-%-2s-%-2s" % (at0,at3,at2,at1)
                if p_improper.has_key(kee3):
                    s_improper[kee3] = p_improper[kee3]
                    break
                kee4 = "%-2s-%-2s-%-2s-%-2s" % (at3,at0,at2,at1)
                if p_improper.has_key(kee4):
                    s_improper[kee4] = p_improper[kee4]
                    break

                kee5 = "%-2s-%-2s-%-2s-%-2s" % (at0,at1,at2,at3)
                if p_improper.has_key(kee5):
                    s_improper[kee5] = p_improper[kee5]
                    break
                kee6 = "%-2s-%-2s-%-2s-%-2s" % (at1,at0,at2,at3)
                if p_improper.has_key(kee6):
                    s_improper[kee6] = p_improper[kee6]
                    break

                kee = "X -%-2s-%-2s-%-2s" % (at3,at2,at0)
                if p_improper.has_key(kee):
                    s_improper[kee1] = p_improper[kee]
                    break
                kee = "X -%-2s-%-2s-%-2s" % (at1,at2,at0)
                if p_improper.has_key(kee):
                    s_improper[kee2] = p_improper[kee]
                    break
                kee = "X -%-2s-%-2s-%-2s" % (at3,at2,at1)
                if p_improper.has_key(kee):
                    s_improper[kee3] = p_improper[kee]
                    break
                kee = "X -%-2s-%-2s-%-2s" % (at0,at2,at1)
                if p_improper.has_key(kee):
                    s_improper[kee4] = p_improper[kee]
                    break
                kee = "X -%-2s-%-2s-%-2s" % (at1,at2,at3)
                if p_improper.has_key(kee):
                    s_improper[kee5] = p_improper[kee]
                    break
                kee = "X -%-2s-%-2s-%-2s" % (at0,at2,at3)
                if p_improper.has_key(kee):
                    s_improper[kee6] = p_improper[kee]
                    break

                kee = "X -X -%-2s-%-2s" % (at2,at0)
                if p_improper.has_key(kee):
                    s_improper[kee1] = p_improper[kee]
                    break
                kee = "X -X -%-2s-%-2s" % (at2,at1)
                if p_improper.has_key(kee):
                    s_improper[kee3] = p_improper[kee]
                    break
                kee = "X -X -%-2s-%-2s" % (at2,at3)
                if p_improper.has_key(kee):
                    s_improper[kee5] = p_improper[kee]
                    break
                self.miss_improper.append((a,kee1))
                break

        if feedback['actions']:
            print ' '+str(self.__class__)+': missing:'
            print ' '+str(self.__class__)+':    mol. wts.         %6d' % (
                len(self.miss_mw))
            print ' '+str(self.__class__)+':    vdw               %6d' % (
                len(self.miss_vdw))
            print ' '+str(self.__class__)+':    bonds             %6d' % (
                len(self.miss_bond))
            print ' '+str(self.__class__)+':    angles            %6d' % (
                len(self.miss_angle))
            print ' '+str(self.__class__)+':    torsions          %6d' % (
                len(self.miss_torsion))
            print ' '+str(self.__class__)+':    impropers         %6d (usually okay)' % (
                len(self.miss_improper))
            print ' '+str(self.__class__)+':    extra tinker info %6d' % (
                len(self.miss_extra))
        
    def dump(self):
        kees = self.mw.keys()
        kees.sort()
        for b in kees:
            print b,self.mw[b]
        kees = self.bond.keys()
        kees.sort()
        for b in kees:
            print b,self.bond[b]
        kees = self.angle.keys()
        kees.sort()
        for b in kees:
            print b,self.angle[b]
        kees = self.torsion.keys()
        kees.sort()
        for b in kees:
            print b,self.torsion[b]
        kees = self.improper.keys()
        kees.sort()
        for b in kees:
            print b,self.improper[b]
        kees = self.vdw.keys()
        kees.sort()
        for b in kees:
            print b,self.vdw[b]

    def dump_missing(self,impropers=0):
        if len(self.miss_mw):
            print ' '+str(self.__class__)+': missing molecular weights...'
        for b in self.miss_mw:
            print " ",b
        if len(self.miss_vdw):
            print ' '+str(self.__class__)+': missing van der Waalss...'
        for b in self.miss_vdw:
            print " ",b
        if len(self.miss_bond):
            print ' '+str(self.__class__)+': missing bonds...'
        for b in self.miss_bond:
            print " ",b
        if len(self.miss_angle):
            print ' '+str(self.__class__)+': missing angless...'
        for b in self.miss_angle:
            print " ",b
        if len(self.miss_torsion):
            print ' '+str(self.__class__)+': missing torsions...'
        for b in self.miss_torsion:
            print " ",b
        if impropers:
            if len(self.miss_improper):
                print ' '+str(self.__class__)+': missing impropers...'
            for b in self.miss_improper:
                print " ",b

    def complete(self):
        if not ( len(self.miss_mw) or
                    len(self.miss_vdw) or
                    len(self.miss_bond) or
                    len(self.miss_angle) or
                    len(self.miss_torsion) ):
            return 1
        else:
            return 0

    def assign_types(self):
        c = 0
        self.mapping = {}
        map = self.mapping
        kees = self.present.keys()
        label = []
        type = []
        kees.sort()
        for a in kees:
            c = c + 1
            st = str(c)
            label.append(st)
            type.append(a)
            map[a] = st
        # assign numeric types 
        for a in self.model.atom:
            a.numeric_type = map[a.text_type]
        
    def write_tinker_prm(self,fname,proofread=None,smooth=None):
        c = 0
        self.mapping = {}
        map = self.mapping
        kees = self.present.keys()
        label = []
        type = []
        kees.sort()
        if proofread:
            for a in kees:
                label.append(a)
                type.append(a)
                map[a] = a
        else:
            for a in kees:
                c = c + 1
                st = str(c)
                label.append(st)
                type.append(a)
                map[a] = st
        # assign numeric types as per the parameter we're writing
        for a in self.model.atom:
            a.numeric_type = map[a.text_type]
        f = open(fname,'w')
        if smooth:
            f.write('''
forcefield              SMOOTH AMBER
vdwtype                 GAUSSIAN
gausstype               LJ-2
''')
        else:
            f.write('''
forcefield              AMBER95
vdwtype                 LENNARD-JONES
''')
        f.write('''
radiusrule              ARITHMETIC
radiustype              R-MIN
radiussize              RADIUS
epsilonrule             GEOMETRIC
vdw-14-use
vdw-scale               2.0
chg-14-use
chg-scale               1.2
dielectric              1.0
''')

        # types
        for c in range(len(self.present)):
            at = type[c]
            f.write("atom %6s %6s    %-2s      %-25s %3d %10.3f%6d\n" % (
                label[c],label[c],at,'"'+self.mw[at][1][0:23]+'"',
                self.extra[at][0],
                self.mw[at][0],self.extra[at][1]))
        # van der waals
        for c in range(len(self.present)):
            at = type[c]
            f.write("vdw   %8s          %10.4f %10.4f\n" % (
                map[at],self.vdw[at][0],self.vdw[at][1]))
        # bonds
        bond = {}
        for a in self.bond.keys():
            kee = (map[string.strip(a[0:2])],map[string.strip(a[3:5])])
            bond[kee] = a
        kees = bond.keys()
        kees.sort()
        for a in kees:
            kee = bond[a]
            f.write("bond    %6s%5s     %10.1f %10.4f\n" % (
                a[0],a[1],self.bond[kee][0],self.bond[kee][1]))
        # angles
        angle = {}
        for a in self.angle.keys():
            kee = (map[string.strip(a[0:2])],
                     map[string.strip(a[3:5])],
                     map[string.strip(a[6:8])],                
                     )
            angle[kee] = a
        kees = angle.keys()
        kees.sort()
        for a in kees:
            kee = angle[a]
            f.write("angle    %5s%5s%5s%9.2f    %8.2f\n" % (
                a[0],a[1],a[2],self.angle[kee][0],self.angle[kee][1]))
        # impropers
        if not smooth:
            improper = {}
            for a in self.improper.keys():
                kee = (map[string.strip(a[0:2])],
                         map[string.strip(a[3:5])],
                         map[string.strip(a[6:8])],
                         map[string.strip(a[9:11])],                                
                         )
                improper[kee] = a
            kees = improper.keys()
            kees.sort()
            for a in kees:
                kee = improper[a]
                f.write("imptors  %5s%5s%5s%5s     %10.3f%7.1f%3d\n" % (
                    a[0],a[1],a[2],a[3],self.improper[kee][0],
                    self.improper[kee][1],self.improper[kee][2]))
        else:
            improper = {}
            for a in self.improper.keys():
                kee = (map[string.strip(a[0:2])],
                         map[string.strip(a[3:5])],
                         map[string.strip(a[6:8])],
                         map[string.strip(a[9:11])],                                
                         )
                improper[kee] = a
            kees = improper.keys()
            kees.sort()
            for a in kees:
                kee = improper[a]
                f.write("improper %5s%5s%5s%5s      %10.2f       0.00\n" % (
                    a[2],a[0],a[1],a[3],self.improper[kee][0]*3.5))
        # torsions
        torsion = {}
        for a in self.torsion.keys():
            kee = (map[string.strip(a[0:2])],
                     map[string.strip(a[3:5])],
                     map[string.strip(a[6:8])],
                     map[string.strip(a[9:11])],                                
                     )
            torsion[kee] = a
        kees = torsion.keys()
        kees.sort()
        for a in kees:
            kee = torsion[a]
            st = "torsion  %5s%5s%5s%5s     " % (
                a[0],a[1],a[2],a[3])
            lst = self.torsion[kee]
            div = lst[0]
            lst = lst[1:]
            while len(lst):
                st = st + "%10.3f%7.1f%3d  " % (
                    lst[0]/div,lst[1],lst[2])
                lst = lst[5:]
            while len(st)>79:
                st = string.replace(st,'  ',' ')
            f.write(st+"\n")
        # null charge records
        for c in range(len(self.present)):
            a = label[c]
            f.write("charge  %6s    %10.4f\n" % (a,0.0))
        f.close()


    def get_list(self):
        list = []
        c = 0
        self.mapping = {}
        map = self.mapping
        kees = self.present.keys()
        label = []
        type = []
        kees.sort()
        for a in kees:
            c = c + 1
            st = str(c)
            type.append(a)
            map[a] = c
        # assign numeric types as per the parameter we're writing
        for a in self.model.atom:
            a.numeric_type = map[a.text_type]
        # types
        typ_list = []
        for c in range(len(self.present)):
            at = type[c]
            typ_list.append((map[at],self.extra[at][0],self.mw[at][0],self.extra[at][1]))
        list.append(typ_list)
        # van der waals
        vdw_list = []
        for c in range(len(self.present)):
            at = type[c]
            vdw_list.append((map[at],self.vdw[at][0],self.vdw[at][1]))
        list.append(vdw_list)
        # bonds
        bnd_list = []
        bond = {}
        for a in self.bond.keys():
            kee = (map[string.strip(a[0:2])],map[string.strip(a[3:5])])
            bond[kee] = a
        kees = bond.keys()
        kees.sort()
        for a in kees:
            kee = bond[a]
            bnd_list.append((a[0],a[1],self.bond[kee][0],self.bond[kee][1]))
        list.append(bnd_list)
        # angles
        ang_list = []
        angle = {}
        for a in self.angle.keys():
            kee = (map[string.strip(a[0:2])],
                     map[string.strip(a[3:5])],
                     map[string.strip(a[6:8])],                
                     )
            angle[kee] = a
        kees = angle.keys()
        kees.sort()
        for a in kees:
            kee = angle[a]
            ang_list.append((a[0],a[1],a[2],self.angle[kee][0],self.angle[kee][1]))
        list.append(ang_list)
        # impropers
        imp_list = []
        improper = {}
        for a in self.improper.keys():
            kee = (map[string.strip(a[0:2])],
                     map[string.strip(a[3:5])],
                     map[string.strip(a[6:8])],
                     map[string.strip(a[9:11])],                                
                     )
            improper[kee] = a
        kees = improper.keys()
        kees.sort()
        for a in kees:
            kee = improper[a]
            imp_list.append((a[0],a[1],a[2],a[3],self.improper[kee][0],
                self.improper[kee][1],self.improper[kee][2]))
        list.append(imp_list)
        # torsions
        tor_list = []
        torsion = {}
        for a in self.torsion.keys():
            kee = (map[string.strip(a[0:2])],
                     map[string.strip(a[3:5])],
                     map[string.strip(a[6:8])],
                     map[string.strip(a[9:11])],                                
                     )
            torsion[kee] = a
        kees = torsion.keys()
        kees.sort()
        for a in kees:
            kee = torsion[a]
            tor = [ a[0],a[1],a[2],a[3] ]
            lst = self.torsion[kee]
            div = lst[0]
            lst = lst[1:]
            while len(lst):
                tor.extend([lst[0]/div,lst[1],lst[2]])
                lst = lst[5:]
            tor_list.append(tor)
        list.append(tor_list)
        return list
