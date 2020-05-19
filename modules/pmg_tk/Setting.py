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

import sys
if True:
    import tkinter as Tkinter
    from tkinter import IntVar
import time

from pymol import cmd, Scratch_Storage

class PymolVar(Tkinter.Variable, object):
    '''Tk variable synced with PyMOL setting'''

    _items = False

    def __init__(self, sinst, index, v):
        self.s = sinst     # pmg_tk.Setting.Setting instance
        self.index = index # setting index
        self.skip_w = -1   # prevent trace_w/update loops

        super(PymolVar, self).__init__()
        self.set(v)
        self.trace_variable('w', self.trace_w)

    def trace_w(self, *args):
        '''Set the PyMOL setting to this var's value'''
        if self.skip_w > 0:
            self.skip_w -= 1
            return
        self.s.set(self.index, self.get())

    def update(self):
        '''Set this var's value to the PyMOL setting's value'''
        self.skip_w += 1
        if not self.skip_w:
            return
        vt = self.s.get_setting_tuple(self.index)[1]
        if self._items:
            self.set(str(vt))
            for item, v in zip(self._items, vt):
                item.skip_w += 1
                item.set(v)
        else:
            self.set(vt[0])

    def __getitem__(self, i):
        '''For 3f settings: provide sub-variables for items'''
        if not self._items:
            self._items = []
        N = len(self._items)
        if N <= i:
            v = self.s.get_setting_tuple(self.index)[1]
            for j in range(N, i+1):
                self._items.append(ListVarItem(j,
                    self.s, self.index, round(v[j], 5)))
        return self._items[i]

class ListVarItem(PymolVar):
    '''Tk variable for 3f setting items'''
    def __init__(self, i, *args, **kwargs):
        self.i = i
        super(ListVarItem, self).__init__(*args, **kwargs)

    def trace_w(self, *args):
        if self.skip_w > 0:
            self.skip_w -= 1
            return
        v = list(self.s.get_setting_tuple(self.index)[1])
        try:
            v[self.i] = float(self.get())
        except ValueError:
            return
        self.s.set(self.index, v)

    def update(self):
        raise NotImplementedError

class ColorVar(PymolVar):
    '''Numeric color variable'''
    def set(self, v):
        if isinstance(v, str) and v.strip():
            v = cmd.get_color_index(v)
        super(ColorVar, self).set(v)

class Setting:
    '''Proxy to global or object level PyMOL settings'''

    def __init__(self,app, sele='', state=0):
        self.cmd = app.pymol.cmd
        self.pymol = app.pymol
        self.sele = sele
        self.state = state
        
        while not self.cmd.ready(): # make sure PyMOL is ready for action...
            time.sleep(0.1)

        self.active_dict = {}

    def set(self, name, value):
        self.cmd.set(name, value, self.sele, self.state, log=1, quiet=0)

    def get(self, name):
        return self.cmd.get(name, self.sele, self.state)

    def get_setting_tuple(self, name):
        r = self.cmd.get_setting_tuple(name, self.sele, self.state)
        if r[0] == 0:
            r = self.cmd.get_setting_tuple(name)
        return r

    def __getattr__(self, name):
        index = self.pymol.setting.index_dict[name]
        v_type, v_list = self.get_setting_tuple(name)

        if v_type < 4: # bool, int, float
            var = PymolVar(self, index, v_list[0])
        elif v_type == 4: # 3f vector
            var = PymolVar(self, index, self.get(name))
        elif v_type == 5: # color
            var = ColorVar(self, index, v_list[0])
        else: # text
            var = PymolVar(self, index, v_list[0])

        setattr(self, name, var)
        self.active_dict[index] = var.update
        return var

    def refresh(self): # get any settings changes from PyMOL and update menus
        for a in self.cmd.get_setting_updates(self.sele, self.state) or ():
            try:
                self.active_dict[a]()
            except KeyError:
                pass
