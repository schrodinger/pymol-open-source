import thread 
import threading 
import os
import sys
import re
import _pm
import string 
import time

sys.path.append(os.environ['PYMOL_PATH']+'/modules')
sys.setcheckinterval(1)

lock_oq = threading.RLock()
lock_iq = threading.RLock()
lock_api = threading.RLock()

glutThread = 0

def start_pymol():
	global glutThread
	glutThread = thread.get_ident()
	_pm.runpymol()
	
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
	elif re.search(r"\.pdb$|\.mol$",a):
		pm.do("load %s" % a)
	elif re.search(r"\.pml$",a):
		pm.do("@%s" % a)

execfile(os.environ['PYMOL_PATH']+'/modules/pmg.py',globals(),locals())

#t2=threading.Thread(target=execfile,args=(os.environ['PYMOL_PATH']+'/modules/pmg.py',globals(),locals()))
#t2.setDaemon(1)
#t2.start()


