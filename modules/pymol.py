import thread 
import threading 
import os
import sys
import re
import string 
import time

sys.path.append(os.environ['PYMOL_PATH']+'/modules')

import _pm

sys.setcheckinterval(1)

lock_oq = threading.RLock()
lock_iq = threading.RLock()
lock_api = threading.RLock()

glutThread = 0

if "-c" in sys.argv:
   gui=0
else:
   gui=1

if "-s" in sys.argv:
   stereo=2
else:
   stereo=0

def start_pymol():
	global glutThread
	glutThread = thread.get_ident()
	_pm.runpymol(gui|stereo)
	
threading.Thread(target=start_pymol,args=()).start()

import pm

while not pm.ready():
	time.sleep(0.1)
time.sleep(0.2)

for a in sys.argv:
	if re.search(r"pymol\.py$",a):
		pass
	elif re.search(r"\.py$",a):
		pm.do("run %s" % a)
	elif re.search(r"\.pdb$|\.mol$|\.mmod",a):
		pm.do("load %s" % a)
	elif re.search(r"\.pml$",a):
		pm.do("@%s" % a)

if gui:
   execfile(os.environ['PYMOL_PATH']+'/modules/pmg.py',globals(),locals())

#t2=threading.Thread(target=execfile,args=(os.environ['PYMOL_PATH']+'/modules/pmg.py',globals(),locals()))
#t2.setDaemon(1)
#t2.start()


