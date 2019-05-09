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

class Map:

    def __init__(self):
        self.cell_dim = [ 1.0, 1.0, 1.0]
        self.cell_ang = [ 90.0, 90.0, 90.0 ]
        self.cell_div = [ 10, 10, 10 ] # how many divisions per cell axis
        self.first = [ 0, 0, 0 ] # first slice on each axis
        self.last = [ 9, 9, 9 ] # last slice on each axis
        self.format = 'Empty'

    def from_c_object(self,cobj,format,cell_dim,cell_ang,cell_div,first,last):
        # pack with defined types
        self.format = str(format)
        self.c_object = cobj
        self.cell_ang = [ float(cell_ang[0]),
                                float(cell_ang[1]),
                                float(cell_ang[2])]
        self.cell_dim = [ float(cell_dim[0]),
                                float(cell_dim[1]),
                                float(cell_dim[2])]
        self.cell_div = [ int(cell_div[0]),
                                int(cell_div[1]),
                                int(cell_div[2])]
        self.first = [ int(first[0]),
                            int(first[1]),
                            int(first[2])]
        self.last = [ int(last[0]),
                          int(last[1]),
                          int(last[2])]
