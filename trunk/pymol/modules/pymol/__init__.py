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

import thread 
import threading 
import os
import sys
import re
import string 
import time
import invocation
import traceback
import _cmd

# PyMOL __init__.py

sys.path.append(os.environ['PYMOL_PATH']+'/modules')

sys.setcheckinterval(1)
lock_api = threading.RLock()

def start_pymol():
	global glutThread
	glutThread = thread.get_ident()
	_cmd.runpymol()
	
def exec_str(s):
   try:
      exec s in globals(),globals()
   except StandardError:
      traceback.print_exc()
   return None
   
def exec_deferred():
#   pm.do("@t.pml")
   for a in invocation.options.deferred:
      if re.search(r"pymol\.py$",a):
         pass
      elif re.search(r"\.py$",a):
         pm.do("run %s" % a)
      elif re.search(r"\.pdb$|\.mol$|\.mmod$|\.mmd$|\.xplor$|\.pkl$|\.sdf$",a):
         pm.load(a)
      elif re.search(r"\.pml$",a):
         pm.do("@%s" % a)

def launch_gui():
   if invocation.options.external_gui:
      __import__(invocation.options.gui)

import cmd

if os.environ.has_key('DISPLAY'):
   from xwin import *
   
