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
import math
import threading

# PyMOL __init__.py


# Create a temporary object "stored" in the PyMOL global namespace
# for usage with evaluate based-commands such as alter

class Scratch_Storage:
   pass

stored = Scratch_Storage()

#

sys.path.append(os.environ['PYMOL_PATH']+'/modules')

sys.setcheckinterval(1)
lock_api = threading.RLock() # mutex for API 
lock_api_c = threading.RLock() # mutex for C management of python threads

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

def stdin_reader(): # dedicated thread for reading standard input
	import sys
	from pymol import cmd
	while 1:
		cmd.do(sys.stdin.readline())
			   
def exec_deferred():
   cmd.config_mouse(quiet=1)
   for a in invocation.options.deferred:
      if a[0:4]=="_do_":
         cmd.do(a[4:])
      elif re.search(r"pymol\.py$",a):
         pass
      elif re.search(r"\.py$|\.pym|\.pyc$",a,re.I):
         cmd.do("run %s" % a)
      elif cmd.file_ext_re.search(a):
         cmd.load(a)
      elif re.search(r"\.pml$",a,re.I):
         cmd.do("@%s" % a)
      else:
         cmd.load(a)
   if invocation.options.read_stdin:
      t = threading.Thread(target=stdin_reader)
      t.setDaemon(1)
      t.start()

def launch_gui():
   if invocation.options.external_gui:
      __import__(invocation.options.gui)

import cmd

if os.environ.has_key('DISPLAY'):
   from xwin import *
