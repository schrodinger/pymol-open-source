# 
# basic run through the PDB, test hydrogen addition and SMILES conversion
#
      
from glob import glob

import threading
import time
from pymol import cmd
import sys, os, os.path
from chempy.champ import Champ

ent_dir = "pdb"

def load():
   cmd.set("valence")
   r = 0
   list = glob("pdb/*/*")

#   list = ["pdb/06/pdb406d"]
#   while list[0]!="pdb/f8/pdb1f8u":
#      list.pop(0)
   for file in list:
      cmd.delete('pdb')
      cmd.load(file,'pdb')
      cmd.remove("not alt ''+A")
      cmd.fix_chemistry("all","all")
      cmd.h_add()
      sys.__stderr__.write(file)
      sys.__stderr__.flush()
      ch = Champ()
      model = cmd.get_model('pdb')
      idx = ch.insert_model(model)

      ch.pattern_orient_bonds(idx)
      print " %5d"%cmd.count_atoms(),"%5d"%len(ch.pattern_get_string(idx)),
      ch.pattern_detect_chirality(idx)
      pat = ch.pattern_get_string(idx)
      print "%5d"%len(pat),pat[0:22]+"..."+pat[-10:]

cmd.feedback('disable','symmetry objectmolecule executive','everything')
load()


