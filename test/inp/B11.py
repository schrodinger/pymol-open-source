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
   (x,y,z)=(0,0,0)
   start = time.time()
   count = 1
   scale = 100.0
   for file in list:
      cmd.load(file,str(count),quiet=1)
      cmd.translate([x*scale,y*scale,z*scale],object=str(count))
#      cmd.disable(str(count))
      atoms = cmd.count_atoms()
      passed = time.time()-start
      print "%3d structures/%5.1f sec = %8.1f atom/sec over %6d atoms"%(count,passed,atoms/passed,atoms)
      count = count + 1
      if count>100: break
      x = x + 1
      if x>7:
         x = 0
         y = y + 1
         if y>7:
            y = 0
            z = z + 1
            if z>9:
               y = 0
               z = z + 1
   cmd.zoom()
   cmd.set("sphere_mode",1)
   cmd.as("spheres")
   cmd.rebuild()
   cmd.set("hash_max",250)
   
#cmd.feedback('disable','symmetry objectmolecule executive','everything')
load()
cmd.set("matrix_mode",1)


