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
#-* NOTE: Based on code by John E. Grayson which was in turn 
#-* based on code written by Doug Hellmann. 
#Z* -------------------------------------------------------------------

from Tkinter import *
import Pmw
from pymol import cmd

pm = cmd

class Setting:

	
	def __init__(self):
	
		self.ray_trace_frames = IntVar()
		self.ray_trace_frames.set(0)

		self.cache_frames = IntVar()
		self.cache_frames.set(0)

		self.ortho = IntVar()
		self.ortho.set(0)

		self.antialias = IntVar()
		self.antialias.set(0)

		self.all_states = IntVar()
		self.all_states.set(0)

		self.overlay = IntVar()
		self.overlay.set(0)

		self.normal_workaround = IntVar()
		self.normal_workaround.set(0)

		self.valence = IntVar()
		self.valence.set(0)

		self.xref = { 
'ray_trace_frames' : (lambda s,a: (
	cmd.set(a,("%1.0f" % s.ray_trace_frames.get())),
	s.cache_frames.set(s.ray_trace_frames.get()),
	s.update('cache_frames'))),
'cache_frames'  : (lambda s,a: (cmd.set(a,("%1.0f" % s.cache_frames.get())))),
'ortho'         : (lambda s,a: (cmd.set(a,("%1.0f" % s.ortho.get())))),
'antialias'     : (lambda s,a: (cmd.set(a,("%1.0f" % s.antialias.get())))),
'valence'       : (lambda s,a: (cmd.set(a,("%1.2f" % (float(s.valence.get())/20.0))))),
'all_states'    : (lambda s,a: (cmd.set(a,("%1.0f" % s.all_states.get())))),
'normal_workaround': (lambda s,a: (cmd.set(a,("%1.0f" % s.normal_workaround.get())))),
'overlay'       : (lambda s,a: (cmd.set(a,("%1.0f" % s.overlay.get()))))
			}
		
	def update(self,sttng):
		set_fn = self.xref[sttng]
		set_fn(self,sttng)


