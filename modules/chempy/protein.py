#
#
#

import protein_residues
import protein_amber

from chempy.neighbor import Neighbor
from chempy.models import Connected
from chempy import Bond
from chempy import place
from chempy.cpv import *

MAX_BOND_LEN = 2.2
PEPT_CUTOFF = 1.7

#---------------------------------------------------------------------------------
def generate(model, forcefield = protein_amber, histidine = 'HIE' ):

   add_bonds(model)
   connected = model.convert_to_connected()
   add_hydrogens(connected)
   place.simple_unknowns(connected)
   return connected.convert_to_indexed()

#---------------------------------------------------------------------------------
def add_bonds(model, forcefield = protein_amber, histidine = 'HIE' ):
   if str(model.__class__) != 'chempy.models.Indexed':
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
               # find out if this is n or c terminal residue
               names = []
               for b in range(a[0],a[1]):
                  names.append(model.atom[b].name)
               tmpl = protein_residues.normal
               if forcefield:
                  ffld = forcefield.normal
               for b in ('HT','HT1','HT2','HT3','H1','H2','H3'):
                  if b in names:
                     tmpl = protein_residues.n_terminal
                     if forcefield:
                        ffld = forcefield.n_terminal
                     break
               for b in ('OXT','O2','OT1','OT2'):
                  if b in names:
                     tmpl = protein_residues.c_terminal
                     if forcefield:
                        ffld = forcefield.c_terminal
                     break
               if not tmpl.has_key(resn):
                  raise RuntimeError("unknown residue type '"+resn+"'")
               else:
                  # reassign atom names and build dictionary
                  dict = {}
                  aliases = tmpl[resn]['aliases']
                  for b in range(a[0],a[1]):
                     at = model.atom[b]
                     if aliases.has_key(at.name):
                        at.name = aliases[at.name]
                     dict[at.name] = b
                     if forcefield:
                        k = (resn,at.name)
                        if ffld.has_key(k):
                           at.text_type = ffld[k]['type']
                           at.partial_charge = ffld[k]['charge']
                        else:
                           raise RuntimeError("no parameters for '"+str(k)+"'")
                  # now add bonds for atoms which are present
                  bonds = tmpl[resn]['bonds']
                  mbond = model.bond
                  for b in bonds.keys():
                     if dict.has_key(b[0]) and dict.has_key(b[1]):
                        bnd = Bond()
                        bnd.index = [ dict[b[0]], dict[b[1]] ]
                        bnd.order = bonds[b]['order']
                        mbond.append(bnd)
                  if dict.has_key('N'):  # connect residues N-C based on distance
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
                  if dict.has_key('SG'): # cysteine
                     cur = dict['SG']
                     at = model.atom[cur]
                     lst = nbr.get_neighbors(at.coord)
                     for b in lst:
                        if b!=cur:
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
                                          k = ('CYX',atx.name)
                                          if ffld.has_key(k):
                                             atx.text_type = ffld[k]['type']
                                             atx.partial_charge = ffld[k]['charge']
                                          else:
                                             raise RuntimeError("no parameters for '"+str(k)+"'")
                                       for d in res_list: # other residue
                                          if (b>=d[0]) and (b<d[1]):
                                             for c in range(d[0],d[1]):
                                                atx = model.atom[c]
                                                atx.resn = 'CYX'
                                                k = ('CYX',atx.name)
                                                if ffld.has_key(k):
                                                   atx.text_type = ffld[k]['type']
                                                   atx.partial_charge = ffld[k]['charge']
                                                else:
                                                   raise RuntimeError("no parameters for '"+str(k)+"'")
                                             
                                    break
                     
#---------------------------------------------------------------------------------
def add_hydrogens(model,forcefield=protein_amber):  # assumes no bonds between non-hetatms
   if str(model.__class__) != 'chempy.models.Connected':
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
               for b in ('HT','HT1','HT2','HT3','H1','H2','H3'):
                  if b in names:
                     tmpl = protein_residues.n_terminal
                     if forcefield:
                        ffld = forcefield.n_terminal
                     break
               for b in ('OXT','O2','OT1','OT2'):
                  if b in names:
                     tmpl = protein_residues.c_terminal
                     if forcefield:
                        ffld = forcefield.c_terminal
                     break
               if not tmpl.has_key(resn):
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
                  for b in bonds.keys():
                     if dict.has_key(b[0]) and (not dict.has_key(b[1])):
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
                     if (not dict.has_key(b[0])) and dict.has_key(b[1]):
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
   model.sort()
