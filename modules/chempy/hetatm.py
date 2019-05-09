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

import chempy.models
from chempy.neighbor import Neighbor
from chempy.models import Connected
from chempy import Bond
from chempy import place

MAX_BOND_LEN = 2.2
PEPT_CUTOFF = 1.7

#---------------------------------------------------------------------------------
def generate(model, topology= None, forcefield = None ):

    add_bonds(model,topology=topology,forcefield=forcefield)
    connected = model.convert_to_connected()
    add_hydrogens(connected,topology=topology,forcefield=forcefield)
    place.simple_unknowns(connected)
    return connected.convert_to_indexed()

#---------------------------------------------------------------------------------
def assign_types(model, topology = None, forcefield = None ):
    if not isinstance(model, chempy.models.Indexed):
        raise ValueError('model is not an "Indexed" model object')
    nAtom = model.nAtom
    if nAtom:
        tmpl = topology.normal
        ffld = forcefield.normal
        res_list = model.get_residues()
        if len(res_list):
            for a in res_list:
                base = model.atom[a[0]]
                resn = base.resn
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

#---------------------------------------------------------------------------------
def add_bonds(model, topology = None, forcefield = None ):
    if not isinstance(model, chempy.models.Indexed):
        raise ValueError('model is not an "Indexed" model object')
    nAtom = model.nAtom
    if nAtom:
        tmpl = topology.normal
        ffld = forcefield.normal
        res_list = model.get_residues()
        if len(res_list):
            for a in res_list:
                base = model.atom[a[0]]
                resn = base.resn
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

#---------------------------------------------------------------------------------
def add_hydrogens(model,topology=None,forcefield=None):
    if not isinstance(model, chempy.models.Connected):
        raise ValueError('model is not a "Connected" model object')
    nAtom = model.nAtom
    if nAtom:
        if not model.index:
            model.update_index()
        ffld = forcefield.normal
        tmpl = topology.normal
        res_list = model.get_residues()
        if len(res_list):
            for a in res_list:
                base = model.atom[a[0]]
                resn = base.resn
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
