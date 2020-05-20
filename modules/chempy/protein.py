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

#
#
#

from . import bond_amber
from . import protein_residues
from . import protein_amber

import chempy.models
from chempy.neighbor import Neighbor
from chempy.models import Connected
from chempy import Bond,place,feedback
from chempy.cpv import *

MAX_BOND_LEN = 2.2
PEPT_CUTOFF = 1.7

N_TERMINAL_ATOMS = ('HT','HT1','HT2','HT3','H1','H2','H3',
                          '1H','2H','3H','1HT','2HT','3HT')

C_TERMINAL_ATOMS = ('OXT','O2','OT1','OT2')

#---------------------------------------------------------------------------------
# NOTE: right now, the only way to get N-terminal residues is to
# submit a structure which contains at least one N_TERMINAL hydrogens

def generate(model, forcefield = protein_amber, histidine = 'HIE',
                 skip_sort=None, bondfield = bond_amber ):

    strip_atom_bonds(model) # remove bonds between non-hetatms (ATOM)
    add_bonds(model,forcefield=forcefield)
    connected = model.convert_to_connected()
    add_hydrogens(connected,forcefield=forcefield,skip_sort=skip_sort)
    place.simple_unknowns(connected,bondfield = bondfield)
    return connected.convert_to_indexed()

#---------------------------------------------------------------------------------
def strip_atom_bonds(model):
    new_bond = []
    matom = model.atom
    for a in model.bond:
        if matom[a.index[0]].hetatm or matom[a.index[1]].hetatm:
            new_bond.append(a)
    model.bond = new_bond

#---------------------------------------------------------------------------------
def assign_types(model, forcefield = protein_amber, histidine = 'HIE' ):
    '''   
assigns types: takes HIS -> HID,HIE,HIP and CYS->CYX where appropriate
but does not add any bonds!
'''
    if feedback['actions']:
        print(" "+str(__name__)+": assigning types...")
    if not isinstance(model, chempy.models.Indexed):
        raise ValueError('model is not an "Indexed" model object')
    if model.nAtom:
        crd = model.get_coord_list()
        nbr = Neighbor(crd,MAX_BOND_LEN)
        res_list = model.get_residues()
        if len(res_list):
            for a in res_list:
                base = model.atom[a[0]]
                if not base.hetatm:
                    resn = base.resn
                    if resn == 'HIS':
                        for c in range(a[0],a[1]): # this residue
                            model.atom[c].resn = histidine
                        resn = histidine
                    if resn == 'N-M': # N-methyl from Insight II,
                        for c in range(a[0],a[1]): # this residue
                            model.atom[c].resn = 'NME'
                        resn = 'NME'
                    # find out if this is n or c terminal residue
                    names = []
                    for b in range(a[0],a[1]):
                        names.append(model.atom[b].name)
                    tmpl = protein_residues.normal
                    if forcefield:
                        ffld = forcefield.normal
                    for b in N_TERMINAL_ATOMS:
                        if b in names:
                            tmpl = protein_residues.n_terminal
                            if forcefield:
                                ffld = forcefield.n_terminal
                            break
                    for b in C_TERMINAL_ATOMS:
                        if b in names:
                            tmpl = protein_residues.c_terminal
                            if forcefield:
                                ffld = forcefield.c_terminal
                            break
                    if resn not in tmpl:
                        raise RuntimeError("unknown residue type '"+resn+"'")
                    else:
                        # reassign atom names and build dictionary
                        dict = {}
                        aliases = tmpl[resn]['aliases']
                        for b in range(a[0],a[1]):
                            at = model.atom[b]
                            if at.name in aliases:
                                at.name = aliases[at.name]
                            dict[at.name] = b
                            if forcefield:
                                k = (resn,at.name)
                                if k in ffld:
                                    at.text_type = ffld[k]['type']
                                    at.partial_charge = ffld[k]['charge']
                                else:
                                    raise RuntimeError("no parameters for '"+str(k)+"'")
                        if 'SG' in dict: # cysteine
                            cur = dict['SG']
                            at = model.atom[cur]
                            lst = nbr.get_neighbors(at.coord)
                            for b in lst:
                                if b>cur: # only do this once (only when b>cur - i.e. this is 1st CYS)
                                    at2 = model.atom[b]
                                    if at2.name=='SG':
                                        if not at2.in_same_residue(at):
                                            dst = distance(at.coord,at2.coord)
                                            if dst<=MAX_BOND_LEN:
                                                if forcefield:
                                                    for c in range(a[0],a[1]): # this residue
                                                        atx = model.atom[c]
                                                        atx.resn = 'CYX'
                                                        resn = atx.resn
                                                        if (c<=b):
                                                            k = ('CYX',atx.name)
                                                            if k in ffld:
                                                                atx.text_type = ffld[k]['type']
                                                                atx.partial_charge = ffld[k]['charge']
                                                            else:
                                                                raise RuntimeError("no parameters for '"+str(k)+"'")
                                                    for d in res_list: # other residue
                                                        if (b>=d[0]) and (b<d[1]):
                                                            for c in range(d[0],d[1]):
                                                                atx = model.atom[c]
                                                                atx.resn = 'CYX'
                                                                # since b>cur, assume assignment later on
                                            break

#---------------------------------------------------------------------------------
def add_bonds(model, forcefield = protein_amber, histidine = 'HIE' ):
    '''
add_bonds(model, forcefield = protein_amber, histidine = 'HIE' )

(1) fixes aliases, assigns types, makes HIS into HIE,HID, or HIP
     and changes cystine to CYX
(2) adds bonds between existing atoms
    '''
    if feedback['actions']:
        print(" "+str(__name__)+": assigning types and bonds...")
    if not isinstance(model, chempy.models.Indexed):
        raise ValueError('model is not an "Indexed" model object')
    if model.nAtom:
        crd = model.get_coord_list()
        nbr = Neighbor(crd,MAX_BOND_LEN)
        res_list = model.get_residues()
        if len(res_list):
            for a in res_list:
                base = model.atom[a[0]]
                if not base.hetatm:
                    resn = base.resn
                    if resn == 'HIS':
                        for c in range(a[0],a[1]): # this residue
                            model.atom[c].resn = histidine
                        resn = histidine
                    if resn == 'N-M': # N-methyl from Insight II,
                        for c in range(a[0],a[1]): # this residue
                            model.atom[c].resn = 'NME'
                        resn = 'NME'
                    # find out if this is n or c terminal residue
                    names = []
                    for b in range(a[0],a[1]):
                        names.append(model.atom[b].name)
                    tmpl = protein_residues.normal
                    if forcefield:
                        ffld = forcefield.normal
                    for b in N_TERMINAL_ATOMS:
                        if b in names:
                            tmpl = protein_residues.n_terminal
                            if forcefield:
                                ffld = forcefield.n_terminal
                            break
                    for b in C_TERMINAL_ATOMS:
                        if b in names:
                            tmpl = protein_residues.c_terminal
                            if forcefield:
                                ffld = forcefield.c_terminal
                            break
                    if resn not in tmpl:
                        raise RuntimeError("unknown residue type '"+resn+"'")
                    else:
                        # reassign atom names and build dictionary
                        dict = {}
                        aliases = tmpl[resn]['aliases']
                        for b in range(a[0],a[1]):
                            at = model.atom[b]
                            if at.name in aliases:
                                at.name = aliases[at.name]
                            dict[at.name] = b
                            if forcefield:
                                k = (resn,at.name)
                                if k in ffld:
                                    at.text_type = ffld[k]['type']
                                    at.partial_charge = ffld[k]['charge']
                                else:
                                    raise RuntimeError("no parameters for '"+str(k)+"'")
                        # now add bonds for atoms which are present
                        bonds = tmpl[resn]['bonds']
                        mbond = model.bond
                        for b in list(bonds.keys()):
                            if b[0] in dict and b[1] in dict:
                                bnd = Bond()
                                bnd.index = [ dict[b[0]], dict[b[1]] ]
                                bnd.order = bonds[b]['order']
                                mbond.append(bnd)
                        if 'N' in dict:  # connect residues N-C based on distance
                            cur_n = dict['N']
                            at = model.atom[cur_n]
                            lst = nbr.get_neighbors(at.coord)
                            for b in lst:
                                at2 = model.atom[b]
                                if at2.name=='C':
                                    if not at2.in_same_residue(at):
                                        dst = distance(at.coord,at2.coord)
                                        if dst<=PEPT_CUTOFF:
                                            bnd=Bond()
                                            bnd.index = [cur_n,b]
                                            bnd.order = 1
                                            mbond.append(bnd)
                                            break
                        if 'SG' in dict: # cysteine
                            cur = dict['SG']
                            at = model.atom[cur]
                            lst = nbr.get_neighbors(at.coord)
                            for b in lst:
                                if b>cur: # only do this once (only when b>cur - i.e. this is 1st CYS)
                                    at2 = model.atom[b]
                                    if at2.name=='SG':
                                        if not at2.in_same_residue(at):
                                            dst = distance(at.coord,at2.coord)
                                            if dst<=MAX_BOND_LEN:
                                                bnd=Bond()
                                                bnd.index = [cur,b]
                                                bnd.order = 1
                                                mbond.append(bnd)
                                                if forcefield:
                                                    for c in range(a[0],a[1]): # this residue
                                                        atx = model.atom[c]
                                                        atx.resn = 'CYX'
                                                        resn = atx.resn
                                                        k = ('CYX',atx.name)
                                                        if k in ffld:
                                                            atx.text_type = ffld[k]['type']
                                                            atx.partial_charge = ffld[k]['charge']
                                                        else:
                                                            raise RuntimeError("no parameters for '"+str(k)+"'")
                                                    for d in res_list:
                                                        if (b>=d[0]) and (b<d[1]): # find other residue
                                                            for c in range(d[0],d[1]):
                                                                atx = model.atom[c]
                                                                atx.resn = 'CYX'
                                                                # since b>cur, assume assignment later on
                                                break

#---------------------------------------------------------------------------------
def add_hydrogens(model,forcefield=protein_amber,skip_sort=None):
    # assumes no bonds between non-hetatms
    if feedback['actions']:
        print(" "+str(__name__)+": adding hydrogens...")
    if not isinstance(model, chempy.models.Connected):
        raise ValueError('model is not a "Connected" model object')
    if model.nAtom:
        if not model.index:
            model.update_index()
        res_list = model.get_residues()
        if len(res_list):
            for a in res_list:
                base = model.atom[a[0]]
                if not base.hetatm:
                    resn = base.resn
                    # find out if this is n or c terminal residue
                    names = []
                    for b in range(a[0],a[1]):
                        names.append(model.atom[b].name)
                    tmpl = protein_residues.normal
                    if forcefield:
                        ffld = forcefield.normal
                    for b in N_TERMINAL_ATOMS:
                        if b in names:
                            tmpl = protein_residues.n_terminal
                            if forcefield:
                                ffld = forcefield.n_terminal
                            break
                    for b in C_TERMINAL_ATOMS:
                        if b in names:
                            tmpl = protein_residues.c_terminal
                            if forcefield:
                                ffld = forcefield.c_terminal
                            break
                    if resn not in tmpl:
                        raise RuntimeError("unknown residue type '"+resn+"'")
                    else:
                        # build dictionary
                        dict = {}
                        for b in range(a[0],a[1]):
                            at = model.atom[b]
                            dict[at.name] = b
                        # find missing bonds with hydrogens
                        bonds = tmpl[resn]['bonds']
                        mbond = model.bond
                        for b in list(bonds.keys()):
                            if b[0] in dict and (b[1] not in dict):
                                at = model.atom[dict[b[0]]]
                                if at.symbol != 'H':
                                    name = b[1]
                                    symbol = tmpl[resn]['atoms'][name]['symbol']
                                    if symbol == 'H':
                                        newat = at.new_in_residue()
                                        newat.name = name
                                        newat.symbol = symbol
                                        k = (resn,newat.name)
                                        newat.text_type = ffld[k]['type']
                                        newat.partial_charge = ffld[k]['charge']
                                        idx1 = model.index[id(at)]
                                        idx2 = model.add_atom(newat)
                                        bnd = Bond()
                                        bnd.index = [ idx1, idx2 ]
                                        bnd.order = bonds[b]['order']
                                        mbond[idx1].append(bnd)
                                        mbond[idx2].append(bnd)
                            if (b[0] not in dict) and b[1] in dict:
                                at = model.atom[dict[b[1]]]
                                if at.symbol != 'H':
                                    name = b[0]
                                    symbol = tmpl[resn]['atoms'][name]['symbol']
                                    if symbol == 'H':
                                        newat = at.new_in_residue()
                                        newat.name = name
                                        newat.symbol = symbol
                                        k = (resn,newat.name)
                                        newat.text_type = ffld[k]['type']
                                        newat.partial_charge = ffld[k]['charge']
                                        idx1 = model.index[id(at)]
                                        idx2 = model.add_atom(newat)
                                        bnd = Bond()
                                        bnd.index = [ idx1, idx2 ]
                                        bnd.order = bonds[b]['order']
                                        mbond[idx1].append(bnd)
                                        mbond[idx2].append(bnd)
        if not skip_sort:
            model.sort()
