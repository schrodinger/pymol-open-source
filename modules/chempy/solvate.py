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

from chempy.models import Indexed
from chempy import water_amber,water_residues
from chempy import hetatm,path,io
import os

import copy

def fill_box(area,water=None,box=18.774349):
    if not water:
        test_path = path + "water.pdb"
        if os.path.exists(test_path):
            water= io.pdb.fromFile(test_path)
            water= hetatm.generate(water,forcefield=water_amber,
                                      topology=water_residues)
    if not water:
        raise RunError('Need water structure file.')
    solv_box = Indexed()
    half = box/2
    x = area[0][0]
    c = 0
    while x<area[1][0]:
        y = area[0][1]
        while y<area[1][1]:
            z = area[0][2]
            while z<area[1][2]:
                tmp = copy.deepcopy(water)
                lr = ''
                for a in tmp.atom:
                    if a.resi!=lr:
                        c = c + 1
                        lr = a.resi
                    a.resi = str(c)
                    b = a.coord
                    b[0] = b[0] + x + half
                    b[1] = b[1] + y + half
                    b[2] = b[2] + z + half
                print " "+__name__+": filling box at %8.3f %8.3f %8.3f\n" % (x,y,z)
                solv_box.merge(tmp)
                z = z + box
            y = y + box
        x = x + box
    return solv_box


