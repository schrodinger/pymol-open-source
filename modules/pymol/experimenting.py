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

import selector
from cmd import _cmd,lock,unlock,Shortcut,QuietException
import cmd
import threading

def expfit(a,b): # Huh?
   try:
      lock()   
      r = _cmd.fit(a,b,2)
   finally:
      unlock()
   return r

def focus():  # BROKEN
   try:
      lock()
      r = _cmd.focus()
   finally:
      unlock()
   return r

def spheroid(object=""):  # EXPERIMENTAL
   try:
      print "Warning: 'spheroid' is experimental, incomplete, and unstable."
      lock()
      r = _cmd.spheroid(str(object))
   finally:
      unlock()
   return r

def mem():
   '''
DESCRIPTION

   "mem" Dumps current memory state to standard output. This is a
   debugging feature, not an official part of the API.

'''
   try:
      lock()
      r = _cmd.mem()
   finally:
      unlock()
   return r

   
def check(selection=None,preserve=0):
# UNSUPPORTED
# This function relies on code that is not currently part of PyMOL/ChemPy
   # NOTE: the realtime module relies on code that is not yet part of PyMOL/ChemPy
   from chempy.tinker import realtime
   if selection==None:
      arg = cmd.get_names("objects")
      arg = arg[0:1]
      if arg:
         if len(arg):
            selection = arg
   if selection!=None:
      selection = selector.process(selection)
      realtime.assign("("+selection+")",int(preserve))
      realtime.setup("("+selection+")")

def fast_minimize(*arg):
# OBSOLETE, TO BE REMOVED
   from chempy.tinker import realtime  
   grad  = 0.01
   iter = 500
   interval = 50
   la = len(arg)
   if not la:
      arg = cmd.get_names("objects")
      arg = arg[0:1]
      la = len(arg)
   if la:
      sele  = "("+arg[0]+")"
      if la>1:
         iter = int(arg[1])
      if la>2:
         grad = float(arg[2])
      if la>3:
         interval = int(arg[3])
      t = threading.Thread(target=realtime.mini,args=(iter,grad,interval,arg[0]))
      t.setDaemon(1)
      t.start()
   
def minimize(*arg):
# OBSOLETE, TO BE REMOVED
   # NOTE: the realtime module relies on code that is not yet part of PyMOL/ChemPy
   from chempy.tinker import realtime  
   grad  = 0.01
   iter = 500
   interval = 50
   la = len(arg)
   if not la:
      arg = cmd.get_names("objects")
      arg = arg[0:1]
      la = len(arg)
   if la:
      sele  = "("+arg[0]+")"
      if la>1:
         iter = int(arg[1])
      if la>2:
         grad = float(arg[2])
      if la>3:
         interval = int(arg[3])
      if realtime.setup(sele):
         t = threading.Thread(target=realtime.mini,args=(iter,grad,interval,arg[0]))
         t.setDaemon(1)
         t.start()
      else:
         print " minimize: missing parameters, can't continue"


def dump(fnam,obj):
   try:
      lock()
      r = _cmd.dump(str(fnam),obj)
   finally:
      unlock()
   return r


def dummy(*arg):
   return None

def test(object,action): # generic test routine for development
   try:
      lock()   
      r=_cmd.test(str(object),int(action))
   finally:
      unlock()
   return r

def import_coords(obj,state,mechio): # experimental
   r = None
   try:
      lock()   
      r = _cmd.import_coords(str(obj),int(state)-1,mechio)
   finally:
      unlock()
   return r

def load_coords(*arg): # UNSUPPORTED
   r = 1
   try:
      lock()
      ok = 1
      ftype = loadable.model
      state = 0
      model = arg[0];
      if len(arg)<2:
         ok=0
      if len(arg)>=2:
         oname = string.strip(arg[1])
      if len(arg)>=3:
         state = int(arg[2])-1
      if ok:
         r = _cmd.load_coords(str(oname),model,
                              int(state)-1,int(ftype))
      else:
         print "Error: invalid arguments."
   finally:
      unlock()
   return r
