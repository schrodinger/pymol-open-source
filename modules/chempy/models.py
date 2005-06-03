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

import chempy
import copy
from cpv import *
import operator

class Base:

#------------------------------------------------------------------------------
    def update_index(self):
        if chempy.feedback['verbose']:
            print " "+str(self.__class__)+": updating indexes..."
        c = 0
        self.index = {}
        idx = self.index
        for a in self.atom:
            idx[id(a)] = c
            c = c + 1
            
#------------------------------------------------------------------------------
    def get_residues(self):
        list = []
        if self.nAtom:
            last = self.atom[0]
            c = 0
            start = 0
            for a in self.atom:
                if not a.in_same_residue(last):
                    list.append((start,c))
                    start = c
                    last = a
                c = c + 1
            if (c-start>1):
                list.append((start,c))
        return list

#------------------------------------------------------------------------------

    def get_coord_list(self):
        lst = []
        for a in self.atom:
            lst.append(a.coord)
        return lst

#------------------------------------------------------------------------------
    def get_mass(self):
        sm = 0.0
        for a in self.atom:
            sm = sm + a.get_mass()
        return sm

#------------------------------------------------------------------------------
    def get_nuclear_charges(self):
        '''Return the sum of nuclear charges of all atoms in a molecule.'''
        sm = 0
        for a in self.atom:
            sm = sm + a.get_number()
        return sm

#------------------------------------------------------------------------------
    def list(self):
        for a in self.atom:
            print a.symbol, a.name,  a.coord
        for a in self.bond:
            print a.index
            
#------------------------------------------------------------------------------
    def get_implicit_mass(self):
        # mass calculation for implicit models

        valence = [0]*len(self.atom)
        implicit = [0]*len(self.atom)
        
        for a in self.bond:
            ai0 = a.index[0]
            ai1 = a.index[1]
            valence[ai0] = valence[ai0] + a.order
            valence[ai1] = valence[ai1] + a.order
        c = 0
        for a in self.atom:
            valence[c] = valence[c] - a.formal_charge
            implicit[c] = a.get_implicit_valence()[valence[c]]
            c = c + 1
        h_count = reduce(operator.__add__,implicit)
        hydrogen = chempy.Atom()
        hydrogen.symbol='H'
        return self.get_mass()+hydrogen.get_mass()*h_count

#------------------------------------------------------------------------------
    def assign_names(self,preserve=1):
        dct = {}
        cnt = {}
        if preserve:
            for a in self.atom:
                if a.has('name'):
                    dct[a['name']] = 1
        else:
            for a in self.atom:
                if hasattr(a,'name'):
                    del a.name
        for a in self.atom:
            if not a.has('name'):
                if not cnt.has_key(a.symbol):
                    c = 1
                else:
                    c = cnt[a.symbol]
                while 1:
                    nam = a.symbol+str(c)
                    c = c + 1
                    if not dct.has_key(nam):
                        break
                cnt[a.symbol]=c
                a.name=nam
                dct[nam]=1
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

class Indexed(Base):

    attr_value = {
        'nAtom' : compile('len(self.atom)','Indexed','eval'),
        'nBond' : compile('len(self.bond)','Indexed','eval'),
        }

    def __getattr__(self,attr):
        if Indexed.attr_value.has_key(attr):
            return eval(Indexed.attr_value[attr])
        else:
            raise AttributeError(attr)

#------------------------------------------------------------------------------
    def __init__(self):
        self.reset()
        
#------------------------------------------------------------------------------
    def reset(self):
        self.index = None
        self.molecule = chempy.Molecule()
        self.atom = []
        self.bond = []

#------------------------------------------------------------------------------
    def get_min_max(self):
        if len(self.atom):
            mn = copy.deepcopy(self.atom[0].coord)
            mx = copy.deepcopy(self.atom[0].coord)
            for a in self.atom:
                ac = a.coord
                if mn[0]>ac[0]: mn[0]=ac[0]
                if mn[1]>ac[1]: mn[1]=ac[1]
                if mn[2]>ac[2]: mn[2]=ac[2]
                if mx[0]<ac[0]: mx[0]=ac[0]
                if mx[1]<ac[1]: mx[1]=ac[1]
                if mx[2]<ac[2]: mx[2]=ac[2]
            return [mn,mx]
        else:
            return [[0.0,0.0,0.0],[0.0,0.0,0.0]]
        
#------------------------------------------------------------------------------
    def merge(self,other): # steals atom objects from 'other' and resets 'other'
        if chempy.feedback['actions']:
            print " "+str(self.__class__)+": merging models..."
        nAtom=self.nAtom
        self.atom.extend(other.atom)
        for b in other.bond:
            b.index[0]=b.index[0]+nAtom
            b.index[1]=b.index[1]+nAtom
            self.bond.append(b)
        other.reset()
        if self.index:
            self.update_index()

#--------------------------------------------------------------------------------
    def delete_atom(self,index):
        if chempy.feedback['atoms']:
            print " "+str(self.__class__)+": deleting atom %d." % index
    
        nAtom=self.nAtom

# update index if it exists
        if self.index:
            idx = self.index
            for k in idx.keys():
                if idx[k] > index:
                    idx[k] = idx[k] - 1
            del idx[id(self.atom[index])]

# delete atom
        del self.atom[index]

# delete bonds associated with this atom

        nBond = len(self.bond)
        templist = []
        for i in range(nBond):
            if index in self.bond[i].index:
                templist.append(i)
        for i in range(len(templist)):
            j = templist[i] - i
            del self.bond[j]

# re-index bond table
        for b in self.bond:
            if b.index[0] > index: 
                b.index[0] = b.index[0] - 1
            if b.index[1] > index: 
                b.index[1] = b.index[1] - 1

#--------------------------------------------------------------------------------
    def delete_list(self,list): # delete a list of indexed atoms

        if chempy.feedback['atoms']:
            print " "+str(self.__class__)+": deleting atoms %s." % str(list)

        nAtom=self.nAtom

        lrev = copy.deepcopy(list)
        lrev.sort()
        lrev.reverse()
        
        # generate cross-reference tables
    
        o2n = {} # old to new
        if len(lrev):
            nxt = lrev.pop()
        else:
            nxt = None
        shft = 0
        for i in range(nAtom):
            if i == nxt:
                o2n[i]=-1
                if len(lrev):
                    nxt = lrev.pop()
                else:
                    nxt = None
                shft = shft - 1
            else:
                ixs = i + shft
                o2n[i] = ixs

        if shft:

            # delete atoms

            new_atom = []
            for i in range(nAtom):
                if o2n[i]>=0:
                    new_atom.append(self.atom[i])
            self.atom = new_atom
            
            # delete bonds 

            new_bond = []   
            for b in self.bond:
                b0 = b.index[0]
                b1 = b.index[1]
                if (o2n[b0]>=0) and (o2n[b1]>=0):
                    b.index[0] = o2n[b0]
                    b.index[1] = o2n[b1]
                    new_bond.append(b)
            self.bond = new_bond

            # update index if it exists
            if self.index:
                self.index = {}
                idx = self.index
                i = 0
                for a in self.atom:
                    idx[id(a)] = i
                    i = i + 1

#------------------------------------------------------------------------------
    def insert_atom(self,index,atom):
        if chempy.feedback['atoms']:
            print " "+str(self.__class__)+': inserting atom "%s" before %d.' % (
                atom.name,index)

        nAtom=self.nAtom
        self.atom.insert(index,atom)
        
# re-index bond table
        for b in self.bond:
            if b.index[0] >= index: 
                b.index[0] = b.index[0] + 1
            if b.index[1] >= index: 
                b.index[1] = b.index[1] + 1

# update index if it exists
        if self.index:
            idx = self.index
            for k in idx.keys():
                if idx[k] >= index:
                    idx[k] = idx[k] + 1
            idx[id(atom)] = index
            
#------------------------------------------------------------------------------
    def index_atom(self,atom):
        c = 0
        id_at = id(atom)
        for a in self.atom:
            if id(a)==id_at:
                return c
            c = c + 1
        return -1
        
#------------------------------------------------------------------------------
    def add_atom(self,atom):
        if chempy.feedback['atoms']:
            print " "+str(self.__class__)+': adding atom "%s".' % atom.name
        self.atom.append(atom)
        index = self.nAtom - 1
        if self.index:
            self.index[id(atom)] = index
        return index

#------------------------------------------------------------------------------
    def add_bond(self,bond):
        if chempy.feedback['bonds']:
            print " "+str(self.__class__)+": adding bond (%d,%d)." % \
                    (bond.index[0],bond.index[1])
        self.bond.append(bond)      

#------------------------------------------------------------------------------
    def remove_bond(self,index):
        if chempy.feedback['bonds']:
            print " "+str(self.__class__)+": removing bond %d." % index
        nBond=len(self.Bond)
        del self.bond[index]
        
    
#------------------------------------------------------------------------------
    def convert_to_connected(self):
        if chempy.feedback['verbose']:
            print " "+str(self.__class__)+": converting to connected model..."
        model = Connected()
        model.molecule = self.molecule
        model.atom = self.atom
        model.bond = []
        model.index = None
        for a in model.atom:
            model.bond.append([])
        for b in self.bond:
            model.bond[b.index[0]].append(b) # note two refs to same object
            model.bond[b.index[1]].append(b) # note two refs to same object 
        self.reset()
        return model
#------------------------------------------------------------------------------
    def from_molobj(self,molobj): 
        self.reset()
        mol = self.molecule
        if len(molobj.title):
            mol.title = molobj.title
        if len(molobj.comments):
            mol.comments = molobj.comments
        mol.chiral = molobj.chiral
        mol.dim_code = molobj.dimcode
        for a in molobj.atom:
            at = chempy.Atom()
            at.symbol = a.symbol
            at.name = a.name
            if a.resn != chempy.Atom.defaults['resn']:
                at.resn = a.resn
            if a.resn_code != chempy.Atom.defaults['resn_code']:
                at.resn_code = a.resn_code
            at.resi = a.resi
            at.resi_number = a.resi_number
            at.b = a.b
            at.q = a.q
            at.alt = a.alt
            at.hetatm = a.hetatm
            if a.segi != chempy.Atom.defaults['segi']:
                at.segi = a.segi
            if a.chain != chempy.Atom.defaults['chain']:
                at.chain = a.chain
            at.color_code = a.color_code
            at.coord = a.coord
            at.formal_charge = a.formal_charge
            at.partial_charge = a.partial_charge
            if a.numeric_type != -99:
                at.numeric_type = a.numeric_type
            if a.text_type != 'UNKNOWN':
                at.text_type = a.text_type
            at.stereo = a.stereo
            if hasattr(a,'flags'):
                at.flags = a.flags
            if hasattr(a,'vdw'):
                at.vdw = a.vdw
            self.atom.append(at)
        for b in molobj.bond:
            bnd = chempy.Bond()
            bnd.index = [b.atom[0],b.atom[1]]
            bnd.order = b.order
            bnd.stereo = b.stereo
            self.bond.append(bnd)
#------------------------------------------------------------------------------
    def sort(self):
        if chempy.feedback['verbose']:
            print " "+__name__+": sorting..."
        if not self.index:
            self.update_index()
        old_index = self.index
        self.atom.sort()      
        self.update_index()
        xref = {}
        new_index = self.index
        for a in new_index.keys():
            xref[old_index[a]] = new_index[a]
        for b in self.bond:
            b.index[0] = xref[b.index[0]]
            b.index[1] = xref[b.index[1]]
        del old_index
        del xref

#------------------------------------------------------------------------------
    def get_internal_tuples(self):
        # generates raw atom sets needed to construct an internal coordinate
        # description of the molecule
        model = self
        # get a connected version too
        cmodel = copy.deepcopy(model).convert_to_connected()
        center = [0.0,0.0,0.0]
        nAtom = model.nAtom
        to_go = nAtom
        done = {}
        if to_go<3:
            z_set = [(0),(1,0)]
        else:
            # get center of molecule
            for a in model.atom:
                center = add(center,a.coord)
            center = scale(center,1.0/nAtom)
            # find most central multivalent atom
            min_a = -1
            c = 0
            for a in model.atom:
                if len(cmodel.bond[c])>1: # must have at least two neighbors
                    d = distance(a.coord,center)
                    if min_a < 0:
                        min_d = d
                        min_a = c
                    elif d<min_d:
                        min_d=d
                        min_a=c
                c = c + 1
            fst = min_a
            # make this our first atom
            z_set = [( fst, )]
            done[fst] = 1
            to_go = to_go - 1
            # for the second atom, try to choose different multivalent neighbor
            nxt = -1
            for b in cmodel.bond[fst]:
                nbr = b.index[0]
                if nbr == fst:
                    nbr = b.index[1]
                if len(cmodel.bond[nbr])>1:
                    nxt = nbr
                    break
            # safety, choose any neighbor
            if nxt<0:
                nbr = b.index[0]
                if nbr == fst:
                    nbr = b.index[1]
                nxt = nbr
            z_set.append((nxt,fst))
            done[nxt] = 1
            to_go = to_go - 1
            # for the third atom, choose a different multivalent neighbor
            trd = -1
            for b in cmodel.bond[fst]:
                nbr = b.index[0]
                if nbr == fst:
                    nbr = b.index[1]
                if len(cmodel.bond[nbr])>1:
                    if not done.has_key(nbr):
                        trd = nbr
                        break
            # safety, choose any unchosen neighbor
            if trd<0:
                for b in cmodel.bond[fst]:
                    nbr = b.index[0]
                    if nbr == fst:
                        nbr = b.index[1]
                    if not done.has_key(nbr):
                        trd = nbr
                        break
            z_set.append((trd,fst,nxt))
            done[trd] = 1
            result = 1
            to_go = to_go - 1
            if to_go:
                # now find all torsions in the molecule
                tors = {}
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
                if len(tors.keys()):
                    # choose remaining atoms based on existing atoms using torsion
                    while to_go:
                        for tor in tors.keys():
                            a0 = tor[0]
                            a1 = tor[1]
                            a2 = tor[2]
                            a3 = tor[3]
                            dh0 = done.has_key(a0)
                            dh1 = done.has_key(a1)
                            dh2 = done.has_key(a2)
                            dh3 = done.has_key(a3)
                            if ( (not dh0) and dh1 and dh2 and dh3 ):
                                z_set.append((a0,a1,a2,a3))
                                done[a0] = 1
                                to_go = to_go - 1
                            elif ( dh0 and dh1 and dh2 and (not dh3) ):
                                z_set.append((a3,a2,a1,a0))
                                done[a3] = 1
                                to_go = to_go - 1
                else: # for molecules with no torsions (dichloromethane, etc.)
                    # we have to generate torsions which include one virtual
                    # bond
                    for b in cmodel.bond[fst]:
                        nbr = b.index[0]
                        if nbr in [fst,nxt,trd]:
                            nbr = b.index[1]
                        if not done.has_key(nbr):
                            z_set.append((nbr,trd,fst,nxt))
                            to_go = to_go - 1
                            done[nbr] = 1
        return z_set
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
class Connected(Base):

    attr_value = {
        'nAtom' : compile('len(self.atom)','Connected','eval'),
        }

    def __getattr__(self,attr):
        if Connected.attr_value.has_key(attr):
            return eval(Connected.attr_value[attr])
        else:
            raise AttributeError(attr)
        
#------------------------------------------------------------------------------
    def __init__(self):
        self.reset()
  
#------------------------------------------------------------------------------
    def reset(self):
        self.index = None
        self.molecule = chempy.Molecule()
        self.atom = []
        self.bond = []
        
#------------------------------------------------------------------------------
    def convert_to_indexed(self):
        if chempy.feedback['verbose']:
            print " "+str(self.__class__)+": converting to indexed model..."
        indexed = Indexed()
        indexed.atom = self.atom
        indexed.molecule = self.molecule
        c = 0 
        for a in self.bond:
            for b in a:
                if b.index[0] == c:
                    indexed.bond.append(b)
            c = c + 1
        self.reset()
        return indexed

#------------------------------------------------------------------------------
    def insert_atom(self,index,atom):
        if chempy.feedback['atoms']:
            print " "+str(self.__class__)+': inserting atom "%s" before %d.' % (
                atom.name,index)

        nAtom=self.nAtom
        self.atom.insert(index,atom)
        
# re-index bond table
        for a in self.bonds:
            for b in a:
                if b.index[0] >= index:
                    b.index[0] = b.index[0] + 1
                if b.index[1] >= index:
                    b.index[1] = b.index[1] + 1

# update index if it exists
        if self.index:
            idx = self.index
            for k in idx.keys():
                if idx[k] >= index:
                    idx[k] = idx[k] + 1
            idx[id(atom)] = index

#------------------------------------------------------------------------------
    def delete_atom(self,index):
        if chempy.feedback['atoms']:
            print " "+str(self.__class__)+": deleting atom %d." % index

        nAtom=self.nAtom

# update index if it exists
        if self.index:
            idx = self.index
            for k in idx.keys():
                if idx[k] > index:
                    idx[k] = idx[k] - 1
            del idx[id(self.atom[index])]

# delete atom
        del self.atom[index]

# delete bonds associated with this atom

        nBond = len(self.bond)

        for a in self.bond:
            i = 0
            templist = []
            for b in a:
                if index in b.index:
                    templist.append(i)
                i = i + 1
            for i in range(len(templist)):
                j = templist[i] - i
                del a[j]

# re-index bond table
        for b in self.bond:
            if b.index[0] > index: 
                b.index[0] = b.index[0] - 1
            if b.index[1] > index: 
                b.index[1] = b.index[1] - 1

#------------------------------------------------------------------------------
    def add_atom(self,atom):
        if chempy.feedback['atoms']:
            print " "+str(self.__class__)+': adding atom "%s".' % atom.name
        self.atom.append(atom)
        self.bond.append([])
        index = self.nAtom - 1
        if self.index:
            self.index[id(atom)] = index
        return index

#------------------------------------------------------------------------------
    def sort(self):
        if chempy.feedback['verbose']:
            print " "+__name__+": sorting..."
        if not self.index:
            self.update_index()
        old_index = self.index
        self.atom.sort()      
        self.update_index()
        xref = {}
        new_index = self.index
        for a in new_index.keys():
            xref[old_index[a]] = new_index[a]
        new_bond = [None] * self.nAtom
        c = 0
        tmp_list = []
        for a in self.bond:
            for b in a:
                if c==b.index[0]:
                    tmp_list.append(b)
            new_bond[xref[c]] = a
            c = c + 1
        for b in tmp_list:
            b.index[0] = xref[b.index[0]]
            b.index[1] = xref[b.index[1]]
        del self.bond
        self.bond = new_bond
        del old_index
        del xref
        



