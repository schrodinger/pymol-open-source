# NOTE: this file is now obsolete, however it has been left intact and
# functional in order to promote backwards compatibility with existing
# PyMOL installs.

import thread 
import threading 
import os
import sys
import time
import __main__

# let pymol/__init__.py known that we're launching using the old way

__main__.pymol_launch = 0 

if hasattr(__main__,"pymol_argv"):
   pymol_argv = __main__.pymol_argv
else:
   pymol_argv = sys.argv

modules_path = os.environ['PYMOL_PATH']+'/modules'

if modules_path not in sys.path:
   sys.path.append(modules_path)

import pymol

pymol.invocation.parse_args(pymol_argv)

pymol.start_pymol()





