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

import numpy

class Brick(object):
    '''
    Map object to load into PyMOL with

    >>> pymol.importing.load_brick(brickinstance, "name")
    '''

    def __init__(self):
        self.valid = None

    @classmethod
    def from_numpy(cls, data, grid, origin=(0.0, 0.0, 0.0)):
        '''
        @param data: numpy float array with len(data.shape) == 3
        @param range: 3f sequence
        @param origin: 3f sequence
        '''
        data = numpy.asfarray(data)
        assert len(data.shape) == 3

        self = cls()
        self.lvl = data
        self.grid = list(grid)
        self.origin = list(origin)

        # redundant information
        self.dim = list(data.shape)
        self.range = [g * (d - 1) for (g, d) in zip(self.grid, self.dim)]

        return self

    def setup_from_min_max(self,mn,mx,grid,buffer=0.0):
        self.origin = [
            mn[0]-buffer,
            mn[1]-buffer,
            mn[2]-buffer
            ]
        self.range = [
            (mx[0]-mn[0])+2*buffer,
            (mx[1]-mn[1])+2*buffer,
            (mx[2]-mn[2])+2*buffer
            ]
        self.dim = [
            1 + int(self.range[0]/grid[0]),
            1 + int(self.range[1]/grid[1]),
            1 + int(self.range[2]/grid[2])
            ]
        self.grid = list(grid)
        self.lvl = numpy.zeros(self.dim,float)
