import thread 
import threading 
import os
import sys
import time
import __main__


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





