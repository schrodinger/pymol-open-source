
from pymol.wizard import Wizard
from pymol import cmd
from chempy import io

import pymol
import os

sele_name = "_mutate"

obj_name = "mutate"

class Mutagenesis(Wizard):

   count = 0
   cutoff = 3.5
   
   def __init__(self):

      Wizard.__init__(self)

      self.library = io.pkl.fromFile(os.environ['PYMOL_PATH']+"modules/chempy/sidechains/sc_library.pkl")
      
      self.status = 0 # 0 no selection, 1 mutagenizing
      self.error = None
      self.object_name = None
      self.mode = 'current'
      
   def get_panel(self):
      return [
         [ 1, 'Mutagenesis',''],
         [ 2, 'Apply' , 'cmd.get_wizard().apply()'],         
         [ 2, 'Clear' , 'cmd.get_wizard().clear()'],
         [ 2, 'Done','cmd.set_wizard()'],
         ]

   def cleanup(self):
      self.clear()
      
   def clear(self):
      cmd.delete(sele_name)
      cmd.delete(obj_name)

   def apply(self):
      pass

   def get_prompt(self):
      self.prompt = None
      if self.status==0:
         self.prompt = [ 'Pick a residue to modify...']
      elif self.status==1:
         self.prompt = [ 'Select a conformation...' ]
      return self.prompt

   def do_library(self):
      pymol.stored.resn=""
      cmd.iterate("(%s and n;ca)"%sele_name,"stored.resn=resn")
      res_type = pymol.stored.resn
      if self.library.has_key(res_type):
         lib = self.library[res_type]
         state = 1
         for a in lib:
            auto_zoom = cmd.get_setting_text('auto_zoom')
            cmd.set('auto_zoom',"0",quiet=1)
            cmd.create(obj_name,sele_name,1,state)
            for b in a.keys():
               if b!='FREQ':
                  cmd.set_dihedral("(%s & n;%s)"%(obj_name,b[0]),
                                   "(%s & n;%s)"%(obj_name,b[1]),
                                   "(%s & n;%s)"%(obj_name,b[2]),
                                   "(%s & n;%s)"%(obj_name,b[3]),
                                   a[b],state=state)
               else:
                  print a[b]
                  cmd.set_title(obj_name,state,"%4.3f"%a[b])
#               print state,b,a[b]
            state = state + 1
            cmd.set('auto_zoom',auto_zoom,quiet=1)
            pymol.util.cbaw(obj_name)
         print " Mutagenesis: %d conformations loaded."%len(self.library[res_type])
         cmd.frame(0)
         cmd.unpick()
         
   def do_pick(self,bondFlag):
      if bondFlag:
         self.error = "Error: please select an atom, not a bond."
         print self.error
      else:
         if self.status!=0:
            cmd.delete(obj_name)
         cmd.select(sele_name,"(byres pk1)")
         cmd.unpick()
         cmd.enable(sele_name)
         self.status = 1
         self.error = None
         self.do_library()
      cmd.refresh_wizard()



