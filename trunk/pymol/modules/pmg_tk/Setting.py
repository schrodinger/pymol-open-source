#A*
#-------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program C*
#copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific.  D*
#-------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.  F*
#-------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information.
#H*
#-------------------------------------------------------------------
#I* Additional authors of this source file include: -* -* NOTE: Based
#on code by John E. Grayson which was in turn -* based on code written
#by Doug Hellmann.  Z*
#-------------------------------------------------------------------

# this section is devoted to making sure that Tkinter variables which
# correspond to Menu-displayed settings are kept synchronized with
# PyMOL

import Tkinter
from Tkinter import IntVar
import time

from pymol import cmd

class PymolVar(Tkinter.Variable, object):
    def __init__(self, index, v):
        super(PymolVar, self).__init__()
        self.set(v)
        self.index = index
        self.skip_w = -1
        self.trace_variable('w', self.trace_w)

    def trace_w(self, *args):
        if self.skip_w > 0:
            self.skip_w -= 1
            return
        cmd.set(self.index, self.get(), log=1)

    def update(self):
        self.skip_w += 1
        if not self.skip_w:
            return
        v = cmd.get_setting_tuple(self.index)[1][0]
        self.set(v)

class ColorVar(PymolVar):
    def set(self, v):
        if isinstance(v, str) and v.strip():
            v = cmd.get_color_index(v)
        super(ColorVar, self).set(v)

class Setting:

    def __init__(self,app):
        self.cmd = app.pymol.cmd
        self.pymol = app.pymol
        
        while not self.cmd.ready(): # make sure PyMOL is ready for action...
            time.sleep(0.1)

        self.active_dict = {}

        self.F=[ None,
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),

                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    ]
        self.SHFTF=[ None,
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),

                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    IntVar(),
                    ]

    def __getattr__(self, name):
        index = self.pymol.setting._get_index(name)
        v_type, v_list = self.cmd.get_setting_tuple(name)

        if v_type < 4: # bool, int, float
            var = PymolVar(index, v_list[0])
        elif v_type == 4: # 3f vector
            raise UserWarning(name)
        elif v_type == 5: # color
            var = ColorVar(index, v_list[0])
        else: # text
            var = PymolVar(index, v_list[0])

        setattr(self, name, var)
        self.active_dict[index] = var
        return var

    def update_scenes(self):
        dict = self.cmd.get_scene_dict()
        if dict != None:
            for x in range(1,13):
                if dict.has_key('F%d'%x):
                    self.F[x].set(1)
                else:
                    self.F[x].set(0)
                if dict.has_key('SHFT-F%d'%x):
                    self.SHFTF[x].set(1)
                else:
                    self.SHFTF[x].set(0)
                
    def refresh(self): # get any settings changes from PyMOL and update menus
        lst = self.cmd.get_setting_updates()
        if lst!=None:
            for a in lst:
                if a in self.active_dict:
                    self.active_dict[a].update()
