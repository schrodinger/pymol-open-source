import thread 
import threading 
import os
import sys
import time

sys.path.append(os.environ['PYMOL_PATH']+'/modules')

import pymol

pymol.invocation.parse_args(sys.argv)

pymol.start_pymol()





