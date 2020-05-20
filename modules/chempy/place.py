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

from chempy.cpv import *
from chempy import feedback
import chempy.models


TET_TAN = 1.41
TRI_TAN = 1.732

#------------------------------------------------------------------------------
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

#------------------------------------------------------------------------------
def simple_unknowns(model,bondfield=bond_amber):
    if feedback['actions']:
        print(" "+str(__name__)+": placing unknowns...")
    # this can be used to build hydrogens and would robably work for
    # acyclic carbons as well
    if not isinstance(model, chempy.models.Connected):
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
                if at1.text_type in bondfield.nonlinear:
                    near = find_known_secondary(model,atx1,know)
                    if near:
                        at3 = model.atom[near[0]]
                        if at3.text_type in bondfield.planer: # Phenolic hydrogens, etc.
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
                        p0 = normalize(sub(at1.coord,at3.coord))
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
                if at1.text_type in bondfield.planer: # guanido, etc
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
                        d1 = sub(at1.coord,at3.coord)
                        p0 = normalize(d1)
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
                        d1 = sub(at1.coord,at3.coord)
                        p0 = normalize(d1)
                        p1 = normalize(cross_product(d2,p0))
                        v = scale(p1,TET_TAN)
                        v = normalize(add(p0,v))
                        at2.coord = add(at1.coord,scale(v,
                            bnd_len[(at1.text_type,at2.text_type)]))
                    else: # blind
                        at2.coord = random_sphere(at1.coord,
                            bnd_len[(at1.text_type,at2.text_type)])

                        # 2013-08-14 added by thomas
                        raise NotImplementedError("FIXME: at3 unassigned")
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

                    # 2013-08-14 added by thomas, not sure if this is correct
                    d1 = sub(at1.coord,at3.coord)
                    p0 = normalize(d1)

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


#------------------------------------------------------------------------------
def test_random():
    '''
    This is a simple test function which drops most coordinates from a
    polypeptide and tries to reposition them with simple_unknowns().

    Works fine to position hydrogens, fails to position other atoms.
    '''
    import random
    from pymol import cmd
    from chempy.champ import assign
    cmd.fab('ACDEFGHIKLMNPQRSTVWY', 'm0')
    assign.amber99()
    for i in range(100):
        m = cmd.get_model('m0').convert_to_connected()
        for a in m.atom:
            if random.random() < 0.8:
                del a.coord
        simple_unknowns(m)
        cmd.load_model(m.convert_to_indexed(), 'm' + str(i + 1))
