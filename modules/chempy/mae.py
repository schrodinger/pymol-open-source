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

# ad-hoc maestro file parser

import re
import copy

strip_re = re.compile(r'\#.*\#')
token_re = re.compile(r'"([^"]*)"|([^ ]+)')
array_re = re.compile(r'(.*)\[([0-9]+)\]')

coerce = {
    's' : str,
    'i' : int,
    'r' : float }

class MAEParser:

    def __init__(self,lst=None):
        self.i = 0 # index in list
        self.t = [] # token list
        self.d = [] # hiearchy of data read
        self.lst = lst
        self.lst_len = len(lst)

    def nxt_lin(self):
        if self.lst:
            if self.i<self.lst_len:
                self.i = self.i + 1
                return self.lst[self.i-1]
        return None

    def nxt_tok(self):
        while 1:
            if len(self.t):
                return self.t.pop(0)
            else:
                l = self.nxt_lin()
                if not l:
                    return None
                l = strip_re.sub('',l).strip()
                self.t = token_re.findall(l)
                self.t = [''.join(x) for x in self.t]
        return None

    def push_tok(self,tok):
        self.t.insert(0,tok)

    def parse_top(self):
        dct = {}
        stk = [] # keyword stack
        mode = 0 # 0 = definition, 1 = data
        while 1:
            tok = self.nxt_tok()
            if tok is None:
                break
            if tok==':::':
                mode = 1
            if not len(stk):
                mode = 0
            if mode:
                lab = stk.pop(0)
                dct[lab] = coerce[lab[0]](*(tok,))
            else:
                stk.append(tok)
            if tok=='}':
                break
        return dct

    def parse_array(self,n_rec): # creates a list of homogenous lists
                                          # each containing data for one field
        dct = {}
        data = [] # actual array data
        stk = [] # keyword stack
        coer = [] # coersion functions
        n_fld = 0
        mode = 0 # 0 = definition, 1 = data
        cc = 0
        while 1:
            tok = self.nxt_tok()
            if tok is None:
                break
            if tok=='}':
                break
            if tok==':::':
                if not mode:
                    mode = 1
                    n_fld = len(stk)
                    c = 0
                    for a in stk: # create row index for the array
                        dct[a] = c
                        data.append([]) # add row for each field
                        coer.append(coerce[a[0]])
                        c = c + 1
            elif not mode:
                stk.append(tok)
            else: # here we actually read the array
                self.push_tok(tok)
                for cc in range(n_rec):
                    tok = self.nxt_tok() # chuck index
                    if tok==':::': # truncated/incomplete array
                        break
                    c = 0
                    for c in range(n_fld):
                        tok = self.nxt_tok()
                        if tok is None:
                            break
                        data[c].append(coer[c](*(tok,)))
                        if tok=='}':
                            break
        return (n_rec,dct,data) # return a tuple

    def parse_m_ct(self):
        dct = {}
        stk = [] # keyword stack
        mode = 0 # 0 = definition, 1 = data
        while 1:
            tok = self.nxt_tok()
            if tok is None:
                break
            if tok==':::':
                mode = 1
            elif mode:
                if not len(stk):
                    mode = 0
                    self.push_tok(tok) # go around
                else:
                    dct[stk.pop(0)] = tok
            else:
                arm = array_re.findall(tok)
                if len(arm):
                    arm = arm[0]
                    n_rec=int(arm[1])
                    if arm[0] in ['m_atom','m_bond']:
                        self.nxt_tok() # skip '{'
                        dct[arm[0]]=self.parse_array(n_rec)
                else:
                    stk.append(tok)
            if tok=='}':
                break
        return dct

    def parse(self):
        while 1:
            tok = self.nxt_tok()
            if tok is None:
                break
            if tok=='{':
                self.d.append('top',self.parse_top())
            elif tok in ['f_m_ct','p_m_ct']:
                self.nxt_tok() # skip '{'
                self.d.append(tok,self.parse_m_ct())
        return self.d

class MAE(Storage):

    def _read_m_atom(self,m_atom,model):
        ma = model.atom
        at_ent = m_atom[1]
        at_dat = m_atom[2]

        nAtom = m_atom[0]
        if 'i_m_mmod_type' in at_ent:
            a1 = at_dat[at_ent['i_m_mmod_type']]

            for b in range(nAtom):
                nt = a1[b]
                ma[b].numeric_type = nt
                ma[b].symbol = MMOD_atom_data[nt][1]
                ma[b].text_type = MMOD_atom_data[nt][0]

        if 'r_m_x_coord' in at_ent and \
            'r_m_y_coord' in at_ent and \
            'r_m_z_coord' in at_ent:
            a1 = at_dat[at_ent['r_m_x_coord']]
            a2 = at_dat[at_ent['r_m_y_coord']]
            a3 = at_dat[at_ent['r_m_z_coord']]
            for b in range(nAtom):
                ma[b].coord = [a1[b],a2[b],a3[b] ]

        if 'i_m_residue_number' in at_ent:
            a1 = at_dat[at_ent['i_m_residue_number']]
            for b in range(nAtom):
                resi = a1[b]
                ma[b].resi = str(resi)
                ma[b].resi_number = resi

        if 's_m_mmod_res' in at_ent:
            a1 = at_dat[at_ent['s_m_mmod_res']]
            for b in range(nAtom):
                ma[b].resi_code = a1[b]

        if 's_m_chain_name' in at_ent:
            a1 = at_dat[at_ent['s_m_chain_name']]
            for b in range(nAtom):
                ma[b].chain = a1[b]

        if 'i_m_color' in at_ent:
            a1 = at_dat[at_ent['i_m_color']]
            for b in range(nAtom):
                ma[b].color_code = a1[b]

        if 'r_m_charge1' in at_ent:
            a1 = at_dat[at_ent['r_m_charge1']]
            for b in range(nAtom):
                ma[b].partial_charge = a1[b]

        if 's_m_pdb_residue_name' in at_ent:
            a1 = at_dat[at_ent['s_m_pdb_residue_name']]
            for b in range(nAtom):
                resn = a1[b].strip()
                if len(resn):
                    ma[b].resn = resn

        if 'i_m_formal_charge' in at_ent:
            a1 = at_dat[at_ent['i_m_formal_charge']]
            for b in range(nAtom):
                ma[b].formal_charge = a1[b]

        if 's_m_atom_name' in at_ent:
            a1 = at_dat[at_ent['s_m_atom_name']]
            for b in range(nAtom):
                nam = a1[b].strip()
                if len(nam):
                    ma[b].name = nam

        if 's_m_pdb_atom_name' in at_ent:
            a1 = at_dat[at_ent['s_m_pdb_atom_name']]
            for b in range(nAtom):
                nam = a1[b].strip()
                if len(nam):
                    ma[b].name = nam

    def _read_m_bond(self,m_bond,model):
        bd_ent = m_bond[1]
        bd_dat = m_bond[2]

        nBond = m_bond[0]

        if len(bd_dat[0]): # not empty right?
            if 'i_m_from' in bd_ent and \
                'i_m_to' in bd_ent and \
                'i_m_order' in bd_ent:
                a1 = bd_dat[bd_ent['i_m_from']]
                a2 = bd_dat[bd_ent['i_m_to']]
                a3 = bd_dat[bd_ent['i_m_order']]
                for b in range(nBond):
                    bd1 = a1[b] - 1
                    bd2 = a2[b] - 1
                    bd3 = a3[b]

                    if bd1<bd2:
                        bnd = Bond()
                        bnd.index = [ bd1,bd2 ]
                        bnd.order = bd3
                        model.bond.append(bnd)

#---------------------------------------------------------------------------------
    def fromList(self,MMODList): # returns a list of indexed models

        mp = MAEParser(lst=MMODList)
        mp_rec = mp.parse()

        full_model = None
        result = []

        for mp_ent in mp_rec:
            if mp_ent[0] == 'f_m_ct':
                f_m_ct = mp_ent[1]
                model = Indexed()
                if 's_m_title' in f_m_ct:
                    model.molecule.title = f_m_ct['s_m_title'].strip()

                if 'm_atom' in f_m_ct:
                    m_atom = f_m_ct['m_atom']
                    nAtom = m_atom[0]
                    for a in range(nAtom):
                        model.atom.append(Atom())
                    self._read_m_atom(m_atom,model)

                if 'm_bond' in f_m_ct:
                    m_bond = f_m_ct['m_bond']
                    self._read_m_bond(m_bond,model)
                full_model = model
                result.append(model)

            elif mp_ent[0]=='p_m_ct' and full_model is not None:
                model = copy.deepcopy(full_model)
                f_m_ct = mp_ent[1]
                if 's_m_title' in f_m_ct:
                    model.molecule.title = f_m_ct['s_m_title'].strip()

                if 'm_atom' in f_m_ct:
                    m_atom = f_m_ct['m_atom']
                    nAtom = m_atom[0]
                    self._read_m_atom(m_atom,model)

                if 'm_bond' in f_m_ct:
                    m_bond = f_m_ct['m_bond']
                    self._read_m_bond(m_bond,model)
                full_model = model
                result.append(model)

        return result

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
