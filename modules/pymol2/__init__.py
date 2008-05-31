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
from copy import deepcopy

import threading
import traceback
import sys

pymol2_lock = threading.RLock() 

from cmd2 import Cmd

class PyMOL:

    def __init__(self,scheme=None): # initialize a PyMOL instance
        pymol2_lock.acquire(1)
        try:

            pymol._init_internals(self)

            self.invocation = self._invocation

            options = self.invocation.options

            if scheme!=None: #
                if scheme == 'presentation':
                    options.quiet = 0
                    options.show_splash = 0
                    options.external_gui = 0
                    options.internal_feedback = 0
                    options.no_quit = 1
                    options.internal_gui = 0
                    options.presentation = 1
                elif scheme == 'widget': # An embedded widget of some type
                    options.quiet = 0
                    options.show_splash = 0
                    options.external_gui = 0
                    options.internal_feedback = 1
                    options.no_quit = 1
            else:
                options.show_splash = 0 # suppress this annoyance by default
                    
            self._COb = _cmd._new(self,self.invocation.options)

            # initialize the cmd API

            self.cmd = Cmd(self,self._COb)
            
            # begin assembling the instance member by member

            # key instance methods

            self.exec_str = pymol.exec_str
            self.adapt_to_hardware = pymol.adapt_to_hardware
            self.exec_deferred = pymol.exec_deferred

            # Python components

            self.util = pymol.util
            self.menu = pymol.menu
            self.setting = pymol.setting
            self.povray = pymol.povray
            self.preset = pymol.preset
            
        except:
            traceback.print_exc()
            pymol2_lock.release()
        
    def __del__(self):
        _cmd._del(self._COb)
        
    def start(self):
        pymol2_lock.acquire()
        try:

            # fire off the C code
            
            _cmd._start(self._COb, self.cmd)

            # add in some additional Python modules
            
            self.chempy = pymol.chempy
            self.bonds = pymol.bonds
            self.models = pymol.models
                
        except:
            traceback.print_exc()            
            pymol2_lock.release()
        
    def startWithTclTk(self, gui = None, skin=None):
        self.start()
        if gui == None:
            gui = self.invocation.options.gui
        if skin == None:
            skin = self.invocation.options.skin
        poll = 0
        __import__(gui)
        sys.modules[gui].__init__(self,poll,skin)
        
    def stop(self):
        _cmd._stop(self._COb)

    def idle(self):
        return _cmd._idle(self._COb)

    def reshape(self, width, height, force=0):
        _cmd._reshape(self._COb,width,height,force)

    def getRedisplay(self,reset):
        return _cmd._getRedisplay(self._COb,reset)

    def draw(self):
        _cmd._draw(self._COb)
    
    def button(self,button,state,x,y,modifiers):
        _cmd._button(self._COb,button,state,x,y,modifiers)

    def drag(self,x,y,modifiers):
        _cmd._drag(self._COb,x,y,modifiers)

        
    
