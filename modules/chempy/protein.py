#
#
#

import protein_residues
import protein_amber
import bond_amber

from chempy.neighbor import Neighbor
from chempy.models import Connected
from chempy import Bond
from cpv import *

MAX_BOND_LEN = 2.2
PEPT_CUTOFF = 1.7
TET_TAN = 1.41
TRI_TAN = 1.732
#---------------------------------------------------------------------------------
def generate(model, forcefield = protein_amber, histidine = 'HIE' ):
   connect(model)
   con_model = Connected(model)
   add_hydrogens(con_model)
   place_unknown(con_model)
   return con_model.convert_to_indexed()

#---------------------------------------------------------------------------------
def connect(model, forcefield = protein_amber, histidine = 'HIE' ):
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
#---------------------------------------------------------------------------------
def find_known_secondary(model,anchor,known_list):
   at = model.atom[anchor]
   h_list = []
   for id in known_list:
      for b in model.bond[id]:
         atx2 = b.index[0]
         if atx2 == id:
            atx2 = b.index[1]
         if atx2 != anchor: # another bonded atom, not achor
            at2 = model.atom[atx2]
            if at2.has('coord'):
               if at2.symbol != 'H':
                  return (id,atx2)
               else:
                  h_list.append((id,atx2))
   if len(h_list): # only return hydrogen as a last resort
      return h_list[0]
   return None

#---------------------------------------------------------------------------------
def place_unknown(model,forcefield=protein_amber,bondfield=bond_amber):
   # this can be used to build hydrogens and would robably work for
   # acyclic carbons as well
   if str(model.__class__) != 'chempy.models.Connected':
      raise ValueError('model is not a "Connected" model object')
   if model.nAtom:
      if not model.index:
         model.update_index()
      idx = model.index
      last_count = -1
      while 1:
         need = [ [], [], [], [] ]
         bnd_len = bondfield.length
   # find known atoms with missing neighbors, and keep track of the neighbors
         for a in model.atom:
            if a.has('coord'):
               miss = []
               know = []
               atx1 = idx[id(a)]
               bnd = model.bond[atx1]
               for b in bnd:
                  atx2 = b.index[0]
                  if atx2 == atx1:
                     atx2 = b.index[1]
                  at2 = model.atom[atx2]
                  if not at2.has('coord'):
                     miss.append(atx2)
                  else:
                     know.append(atx2)
               c = len(miss)
               if c:
                  need[c-1].append((atx1,miss,know))

         for a in need[0]: # missing only one atom
            atx1 = a[0]
            at1 = model.atom[atx1]
            atx2 = a[1][0]
            at2 = model.atom[atx2]
            know = a[2]
            if bondfield.nonlinear.has_key(at1.text_type):
               near = find_known_secondary(model,atx1,know)
               if near:
                  at3 = model.atom[near[0]]
                  if bondfield.planer.has_key(at3.text_type): # Phenolic hydrogens, etc.
                     at4 = model.atom[near[1]]
                     d1 = sub(at1.coord,at3.coord)
                     p0 = normalize(d1)
                     d2 = sub(at4.coord,at3.coord)
                     p1 = normalize(cross_product(d2,p0))
                     p2 = normalize(cross_product(p0,p1))                     
                     v = scale(p2,TRI_TAN)
                     v = normalize(add(p0,v))
                     at2.coord = add(at1.coord,scale(v,
                        bnd_len[(at1.text_type,at2.text_type)]))
                  else: # Ser, Cys, Thr hydroxyl hydrogens
                     at4 = model.atom[near[1]]
                     d2 = sub(at3.coord,at4.coord)
                     v = normalize(d2)
                     at2.coord = add(at1.coord,scale(v,
                        bnd_len[(at1.text_type,at2.text_type)]))
               elif len(know):
                  d2 = [1.0,0,0]
                  at3 = model.atom[know[0]]
                  p0 = normalize(sub(at1.coord,at3))
                  p1 = normalize(cross_product(d2,p0))
                  v = scale(p1,TET_TAN)
                  v = normalize(add(p0,v))
                  at2.coord = add(at1.coord,scale(v,
                        bnd_len[(at1.text_type,at2.text_type)]))
               else:
                  at2.coord = random_sphere(at1.coord,
                      bnd_len[(at1.text_type,at2.text_type)])
            elif len(know): # linear sum...amide, tbu, etc
               v = [0.0,0.0,0.0]
               for b in know:
                  d = sub(at1.coord,model.atom[b].coord)
                  v = add(v,normalize(d))
               v = normalize(v)
               at2.coord = add(at1.coord,scale(v,
                   bnd_len[(at1.text_type,at2.text_type)]))
            else:
               at2.coord = random_sphere(at1.coord,
                    bnd_len[(at1.text_type,at2.text_type)])

         for a in need[1]: # missing two atoms
            atx1 = a[0]
            at1 = model.atom[atx1]
            atx2 = a[1][0]
            at2 = model.atom[atx2]
            know = a[2]
            if bondfield.planer.has_key(at1.text_type): # guanido, etc
               near = find_known_secondary(model,atx1,know)
               if near: # 1-4 present
                  at3 = model.atom[near[0]]
                  at4 = model.atom[near[1]]
                  d1 = sub(at1.coord,at3.coord)
                  p0 = normalize(d1)
                  d2 = sub(at4.coord,at3.coord)
                  p1 = normalize(cross_product(d2,p0))
                  p2 = normalize(cross_product(p0,p1))
                  v = scale(p2,TRI_TAN)
                  v = normalize(add(p0,v))
                  at2.coord = add(at1.coord,scale(v,
                    bnd_len[(at1.text_type,at2.text_type)]))                                                         
                  at2 = model.atom[a[1][1]]
                  v = scale(p2,-TRI_TAN)
                  v = normalize(add(p0,v))
                  at2.coord = add(at1.coord,scale(v,
                    bnd_len[(at1.text_type,at2.text_type)]))
               elif len(know): # no 1-4 found
                  d2 = [1.0,0,0]
                  at3 = model.atom[know[0]]
                  p1 = normalize(cross_product(d2,p0))
                  p2 = normalize(cross_product(p0,p1))
                  v = scale(p2,TRI_TAN)
                  v = normalize(add(p0,v))
                  at2.coord = add(at1.coord,scale(v,
                    bnd_len[(at1.text_type,at2.text_type)]))
                  at2 = model.atom[a[1][1]]
                  v = scale(p2,-TRI_TAN)
                  v = normalize(add(p0,v))
                  at2.coord = add(at1.coord,scale(v,
                    bnd_len[(at1.text_type,at2.text_type)]))
               else:
                  at2.coord = random_sphere(at1.coord,
                     bnd_len[(at1.text_type,at2.text_type)])
            elif len(know)>=2: # simple tetrahedral
               at3 = model.atom[know[0]]
               at4 = model.atom[know[1]]
               v = [0.0,0.0,0.0]
               d1 = sub(at1.coord,at3.coord)
               d2 = sub(at1.coord,at4.coord)
               v = add(normalize(d1),normalize(d2))
               p0 = normalize(v)
               p1 = normalize(cross_product(d2,p0))
               v = scale(p1,TET_TAN)
               v = normalize(add(p0,v))
               at2.coord = add(at1.coord,scale(v,
                     bnd_len[(at1.text_type,at2.text_type)]))
               at2 = model.atom[a[1][1]]               
               v = scale(p1,-TET_TAN)
               v = normalize(add(p0,v))
               at2.coord = add(at1.coord,scale(v,
                     bnd_len[(at1.text_type,at2.text_type)]))
            else:
               if len(know): # sulfonamide? 
                  d2 = [1.0,0,0]
                  at3 = model.atom[know[0]]                  
                  p1 = normalize(cross_product(d2,p0))
                  v = scale(p1,TET_TAN)
                  v = normalize(add(p0,v))
                  at2.coord = add(at1.coord,scale(v,
                     bnd_len[(at1.text_type,at2.text_type)]))
               else: # blind
                  at2.coord = random_sphere(at1.coord,
                     bnd_len[(at1.text_type,at2.text_type)])
               at4=at2
               at2=model.atom[a[1][1]]
               v = [0.0,0.0,0.0]
               d1 = sub(at1.coord,at3.coord)
               d2 = sub(at1.coord,at4.coord)
               v = add(normalize(d1),normalize(d2))
               p0 = normalize(v)
               p1 = normalize(cross_product(d2,p0))
               v = scale(p1,TET_TAN)
               v = normalize(add(p0,v))
               at2.coord = add(at1.coord,scale(v,
                  bnd_len[(at1.text_type,at2.text_type)]))

         for a in need[2]: # missing 3 atoms
            atx1 = a[0]
            at1 = model.atom[atx1]
            atx2 = a[1][0]
            at2 = model.atom[atx2]
            know = a[2]
            near = find_known_secondary(model,atx1,know)
            if near: # 1-4 present
               at3 = model.atom[near[0]]
               at4 = model.atom[near[1]]
               d1 = sub(at1.coord,at3.coord)
               p0 = normalize(d1)
               d2 = sub(at4.coord,at3.coord)
               p1 = normalize(cross_product(d2,p0))
               p2 = normalize(cross_product(p0,p1))
               v = scale(p2,-TET_TAN)
               v = normalize(add(p0,v))
               at2.coord = add(at1.coord,scale(v,
                    bnd_len[(at1.text_type,at2.text_type)]))                                                         
               at4 = at2
               at2 = model.atom[a[1][1]]
               d1 = sub(at1.coord,at3.coord)
               d2 = sub(at1.coord,at4.coord)
               v = add(normalize(d1),normalize(d2))
               p0 = normalize(v)
               p1 = normalize(cross_product(d2,p0))
               v = scale(p1,TET_TAN)
               v = normalize(add(p0,v))
               at2.coord = add(at1.coord,scale(v,
                     bnd_len[(at1.text_type,at2.text_type)]))
               at2 = model.atom[a[1][2]]               
               v = scale(p1,-TET_TAN)
               v = normalize(add(p0,v))
               at2.coord = add(at1.coord,scale(v,
                     bnd_len[(at1.text_type,at2.text_type)]))
            elif len(know): # fall-back
               d2 = [1.0,0,0]
               at3 = model.atom[know[0]]                  
               p1 = normalize(cross_product(d2,p0))
               v = scale(p1,TET_TAN)
               v = normalize(add(p0,v))
               at2.coord = add(at1.coord,scale(v,
                  bnd_len[(at1.text_type,at2.text_type)]))
               at4=at2
               at2=model.atom[a[1][1]]
               v = [0.0,0.0,0.0]
               d1 = sub(at1.coord,at3.coord)
               d2 = sub(at1.coord,at4.coord)
               v = add(normalize(d1),normalize(d2))
               p0 = normalize(v)
               p1 = normalize(cross_product(d2,p0))
               v = scale(p1,TET_TAN)
               v = normalize(add(p0,v))
               at2.coord = add(at1.coord,scale(v,
                  bnd_len[(at1.text_type,at2.text_type)]))
               at2=model.atom[a[1][2]]
               v = scale(p1,-TET_TAN)
               v = normalize(add(p0,v))
               at2.coord = add(at1.coord,scale(v,
                  bnd_len[(at1.text_type,at2.text_type)]))
            else: # worst case: add one and get rest next time around
               at2.coord=random_sphere(at2.coord,
                  bnd_len[(at1.text_type,at2.text_type)])

         for a in need[3]: # missing 4 atoms
            atx1 = a[0]
            at1 = model.atom[atx1]
            atx2 = a[1][0]
            at2 = model.atom[atx2]
            # add coordinate and get the rest next time around
            at2.coord=random_sphere(at2.coord,
                bnd_len[(at1.text_type,at2.text_type)])

         c = 0
         for a in model.atom:
            if not a.has('coord'):
               c = c + 1
         if not c:
            break;
         if c==last_count:
            break;
         last_count = c




