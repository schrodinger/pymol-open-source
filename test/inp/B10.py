# 
# fastest run through the PDB
#
      
from glob import glob

import threading
import time
from pymol import cmd
import sys, os, os.path
import traceback

ent_dir = "pdb"


def load():
   cmd.set("valence")
   cmd.unset("auto_zoom")
   cmd.zoom("center",100)
   r = 0
   list = glob("pdb/*/*")
#   while list[0]!="pdb/f8/pdb1f8u":
#      list.pop(0)
   for file in list:
      try:
         cmd.delete('pdb')
         cmd.load(file,'pdb')
         cmd.set_title('pdb',1,os.path.split(file)[-1])
         cmd.refresh()
      except:
         traceback.print_exc()
cmd.feedback('disable','symmetry objectmolecule executive','everything')
load()


