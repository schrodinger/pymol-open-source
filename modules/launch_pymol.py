import thread 
import threading 
import os
import sys
import time

modules_path = os.environ['PYMOL_PATH']+'/modules'

if modules_path not in sys.path:
   sys.path.append(modules_path)

import pymol

pymol.invocation.parse_args(sys.argv)

pymol.start_pymol()





