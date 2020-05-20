import os
import shutil
import glob
import re
import string
import sys
import time
import socket  # for gethostname()
import getpass # for getuser()
import re

from chempy import feedback

# "do" is the preferred command for running tinker

def do(run_prefix):
    if feedback['bmin']:
        print(" "+str(__name__)+': launching %s...'%bmin_path)
        if hasattr(sys.stdout,"flush"):
            sys.stdout.flush()
    for a in [ ".m1", ".m2", ".log", ".out" ]:
        if os.path.exists(run_prefix + a):
            os.unlink(run_prefix + a)
    pth = "."
#   pth = os.getcwd()
#   pth = re.sub(r".*\/"+getpass.getuser()+"\/",'',pth)
#   cmd = "rsh "+socket.gethostname()+" "+bmin_path+" "+pth+"/"+run_prefix
    cmd = bmin_path+" "+pth+"/"+run_prefix
    print(cmd)
    os.system(cmd)
    while 1:
        if os.path.exists(run_prefix+".out"): break
        time.sleep(0.1)
    if feedback['bmin']:
        os.system("cat bmintmp.log")
        print(" "+str(__name__)+': bmin job complete. ')
        if hasattr(sys.stdout,"flush"):
            sys.stdout.flush()

if 'SCHRODINGER' in os.environ:
    base = os.environ['SCHRODINGER']
    bmin_path = base + '/bmin'
#   os.environ['MMSHARE_EXEC'] = '/apps/schrodinger/mmshare-v10028/bin/Linux-x86'
#   os.environ['MMOD_EXEC'] = '/apps/schrodinger/macromodel-v71008/bin/Linux-x86'
#   bmin_path = os.environ['MMOD_EXEC']+"/bmin"
else:
    bmin_path = "bmin"
