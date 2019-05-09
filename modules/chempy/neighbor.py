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

import copy

class Neighbor:

    def __init__(self,vect_list,spacing):
        self.voxel = {}
        self.neighbor = None
        self.spacing = spacing
        c = 0
        voxel = self.voxel
        for v in vect_list:
            k = self.address(v)
            if k not in voxel:
                voxel[k] = [c]
            else:
                voxel[k].append(c)
            c = c + 1

    def optimize(self):
        self.neighbor = {}
        voxel = self.voxel
        for k in voxel.keys():
            lst = self.neighbor[k]
            for a in (k[0]-1,k[0],k[0]+1):
                for b in (k[1]-1,k[1],k[1]+1):
                    for c in (k[2]-1,k[2],k[2]+1):
                        k2 = (a,b,c)
                        if k2 in voxel:
                            lst.extend(voxel[k2])

    def address(self,vect):
        return (int(vect[0]/self.spacing),
                  int(vect[1]/self.spacing),
                  int(vect[2]/self.spacing))

    def get_voxel(self,vect):
        k = self.address(vect)
        if k in self.voxel:
            return self.voxel[k]
        else:
            return []

    def get_neighbors(self,vect):
        if self.neighbor:
            k = self.address(vect)
            if k in self.neighbor:
                return self.neighbor[k]
            else:
                return []
        else:
            voxel = self.voxel
            k = self.address(vect)
            lst = []
            for a in (k[0]-1,k[0],k[0]+1):
                for b in (k[1]-1,k[1],k[1]+1):
                    for c in (k[2]-1,k[2],k[2]+1):
                        k2 = (a,b,c)
                        if k2 in voxel:
                            lst.extend(voxel[k2])
            return lst
