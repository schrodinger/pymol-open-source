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

# invocation.py
#
# This module unifies argument handling for embedded and modular PyMOL
#

import copy
import re
import os
import glob

pattern = '.pymolrc*'


class generic:
   pass

options = generic();

options.stereo_capable = 0
options.deferred = []
options.no_gui = 0
options.internal_gui = 1
options.internal_feedback = 1
options.external_gui = 1
options.gui = 'pmg_tk'
options.show_splash = 1
options.read_stdin = 0

py_re = re.compile(r"\.py$|\.pym$\.PY$|\.PYM$")

def get_user_config():
   lst = glob.glob(pattern)
   if not len(lst):
      if os.environ.has_key("HOME"):
         lst = glob.glob(os.environ['HOME']+"/"+pattern)
   if not len(lst):
      if os.environ.has_key("PYMOL_PATH"):
         lst = glob.glob(os.environ['PYMOL_PATH']+"/"+pattern)
   first = []
   second = []
   for a in lst:
      if py_re.search(a):
         first.append("_do__ run "+a) # preceeding "_ " cloaks 
      else:
         second.append("_do__ @"+a) # preceeding "_ " cloaks 
   first.sort()
   second.sort()
   return first+second

def parse_args(argv):
   av = copy.deepcopy(argv)
   av = av[1:] # throw out the executable path
   av.reverse()
   global options
   options.deferred = []
   # append user settings file as an option
   options.deferred.extend(get_user_config())
   while 1:
      if not len(av):
         break
      a = av.pop()
      a = re.sub(r'''^"|"$|^'|'$''','',a) # strip extra quotes
      if a[0:1]=='-':
         if "c" in a:
            options.no_gui=1
            options.external_gui=0
         if "s" in a:
            pass # stereo now autodetected
         if "q" in a:
            options.show_splash = 0
         if "i" in a:
            options.internal_gui = 0
         if "f" in a:
            options.internal_feedback = int(av.pop())
         if "x" in a:
            options.external_gui = 0
         if "t" in a:
            options.gui = 'pmg_tk'
         if "w" in a:
            options.gui = 'pmg_wx'
         if "d" in a:
            options.deferred.append("_do_%s"%av.pop())
         if "l" in a:
            options.deferred.append("_do_spawn %s"%av.pop())
         if "r" in a:
            options.deferred.append("_do_run %s"%av.pop())
         if "p" in a:
            options.read_stdin = 1 
      else:
         options.deferred.append(a)


