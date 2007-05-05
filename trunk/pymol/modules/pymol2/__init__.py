#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2007 by Warren Lyford Delano of DeLano Scientific.
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

import __main__
__main__.pymol_launch = 5 

import pymol
from pymol import _cmd
__main__.pymol = pymol

from pymol import selector
import pymol.menu
import pymol.povray

import threading
import traceback

pymol2_lock = threading.RLock() 

from cmd2 import Cmd

class PyMOL:

    def __init__(self): # initialize a PyMOL instance

        pymol2_lock.acquire(1)
        try:

            pymol._init_internals(self)
            
            self._COb = _cmd._new(self)

            # initialize the cmd API

            self.cmd = Cmd(self._COb)
            self.cmd._pymol = self

            # begin assembling the instance member by member

            # key instance methods

            self.exec_str = pymol.exec_str
            self.adapt_to_hardware = pymol.adapt_to_hardware
            self.exec_deferred = pymol.exec_deferred

            self.util = pymol.util

            # Python components

        except:
            traceback.print_exc()
            pymol2_lock.release()
        
    def __del__(self):
        _cmd._del(self._COb)

    def start(self):
        pymol2_lock.acquire()
        try:
            _cmd._start(self._COb, self.cmd)

            # add additional properties from the pymol module

            self.menu = pymol.menu
            self.setting = pymol.setting
            self.povray = pymol.povray
            self.preset = pymol.preset

            self.chempy = pymol.chempy
            self.bonds = pymol.bonds
            self.models = pymol.models
        except:
            traceback.print_exc()            
            pymol2_lock.release()
        
    def stop(self):
        _cmd._stop(self._COb)



        
    
