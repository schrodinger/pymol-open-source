# 
#
      
from glob import glob

import threading
import time
from pymol import cmd
import sys, os, os.path
from chempy.champ import Champ
from pymol import preset

ent_dir = "pdb"

cmd.set("auto_zoom","off")

def load():
   cmd.set("valence")
   r = 0
   list = glob("pdb/*/*")

   start = time.time()
   count = 0
   for file in list:
      cmd.load(file,str(count),quiet=1)
      cmd.refresh()
      cmd.disable(str(count))
      count = count + 1
      atoms = cmd.count_atoms()
      passed = time.time()-start
      print "%3d structures/%5.1f sec = %8.1f atom/sec over %6d atoms"%(count,passed,atoms/passed,atoms)
      if count>500: break

#cmd.feedback('disable','symmetry objectmolecule executive','everything')
load()


