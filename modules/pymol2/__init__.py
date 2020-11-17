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

import pymol

from pymol import _cmd

import threading
import sys

pymol2_lock = threading.RLock()

##
## FIXME: The PyMOL and SingletonPyMOL classes are partly redundant with the
## instance tracking of the "cmd" module (and the pymol2.cmd2.Cmd class),
## which also holds the _COb pointer.
##

class SingletonPyMOL:
    '''
    Start an exclusive PyMOL instance, only one instance allowed
    '''
    def idle(self):
        return _cmd._idle(self._COb)

    def getRedisplay(self, reset=True):
        return _cmd._getRedisplay(self._COb, reset)

    def reshape(self, width, height, force=0):
        _cmd._reshape(self._COb, width, height, force)

    def draw(self):
        _cmd._draw(self._COb)

    def button(self, button, state, x, y, modifiers):
        _cmd._button(self._COb, button, state, x, y, modifiers)

    def drag(self, x, y, modifiers):
        _cmd._drag(self._COb, x, y, modifiers)

    def start(self):
        cmd = pymol.cmd
        if cmd._COb is not None:
            raise RuntimeError('can only start SingletonPyMOL once')

        with pymol2_lock:
            cmd._COb = _cmd._new(pymol, pymol.invocation.options, True)
            _cmd._start(cmd._COb, cmd)

        # this instance tracking is redundant with the "cmd" module itself
        self._COb = cmd._COb
        self.cmd = cmd

    def stop(self):
        with pymol2_lock:
            _cmd._stop(self._COb)

        pymol.cmd._COb = None

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exc_type, exc_value, tb):
        self.stop()


class PyMOL(SingletonPyMOL):
    '''
    Start a non-exclusive PyMOL instance, multiple instances are possible
    '''

    def __getattr__(self, key):
        # Make this a proxy to the "pymol" module.
        return getattr(pymol, key)

    def __init__(self,scheme=None): # initialize a PyMOL instance
        from .cmd2 import Cmd

        with pymol2_lock:
            pymol._init_internals(self)

            self.invocation = self._invocation

            options = self.invocation.options

            if scheme is not None: #
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

            self.glutThread = None

    def __del__(self):
        self.cmd.__dict__.clear()

    def start(self):
        with pymol2_lock:
            _cmd._start(self._COb, self.cmd)

    def stop(self):
        _cmd._stop(self._COb)
