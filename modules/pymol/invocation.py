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

class generic:
   pass

options = generic();

options.stereo_capable = 0
options.deferred = []
options.no_gui = 0
options.internal_gui = 1
options.external_gui = 1
options.gui = 'pmg_tk'
options.show_splash = 1
options.read_stdin = 0

def parse_args(argv):
   av = copy.deepcopy(argv)
   av.reverse()
   global options
   options.deferred = []
   while 1:
      if not len(av):
         break
      a = av.pop()
      if a[0:1]=='-':
         if "c" in a:
            options.no_gui=1
            options.external_gui=0
         if "s" in a:
            options.stereo_capable = 2
         if "q" in a:
            options.show_splash = 0
         if "i" in a:
            options.internal_gui = 0
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


