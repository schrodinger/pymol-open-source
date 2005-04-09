# 
# basic run through the PDB with ray-tracing
#
      
from glob import glob

import threading
import time

from pymol import cmd
import sys, os, os.path
import pymol


ent_dir = "pdb"

pymol.run = 1

cmd.set_key('F1',lambda p=pymol:(setattr(p,'run',0)))
cmd.set_key('F2',lambda p=pymol:(setattr(p,'run',1)))

def load():
   r = 0
   rep = [ "simple", "technical", "ligands", "ligand_sites","ligand_sites_trans","ligand_sites_mesh","pretty", "publication" ]
   list = glob("pdb/*/*")
#   while list[0]!="pdb/f8/pdb1f8u":
#      list.pop(0)
   for file in list:
      print file
      cmd.delete('pdb')
      cmd.load(file,'pdb')
      cmd.orient('pdb')
      cur_rep = rep.pop(0)
      rep = rep + [cur_rep]
      getattr(pymol.preset,cur_rep)('pdb')
      cmd.refresh()
# give PyMOL a chance
      time.sleep(0.02)
      time.sleep(0.02)
      time.sleep(0.02)
      time.sleep(0.02)
      cmd.refresh()
      time.sleep(0.02)
      time.sleep(0.02)
      time.sleep(0.02)
      time.sleep(0.02)

      while(pymol.run==0):
         time.sleep(0.1)
load()


