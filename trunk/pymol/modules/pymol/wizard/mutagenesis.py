
from pymol.wizard import Wizard
from pymol import cmd
from chempy import io

import pymol
import os
import string

sele_name = "_mutate"

obj_name = "mutate"
tmp_name = "tmp_mut"

class Mutagenesis(Wizard):

   count = 0
   cutoff = 3.5
   
   def __init__(self):

      Wizard.__init__(self)

      self.library = io.pkl.fromFile(os.environ['PYMOL_PATH']+"/modules/chempy/sidechains/sc_library.pkl")
      
      self.status = 0 # 0 no selection, 1 mutagenizing
      self.error = None
      self.object_name = None
      self.modes = [
         'Current'
         ]
      self.mode = self.modes[0]
      residues = self.library.keys()
      residues.sort()
      self.modes.extend(residues)
      self.mode_name={}
      for a in self.modes:
         print "[",a,"]"
         self.mode_name[a] = a
      
      smm = []
      smm.append([ 2, 'Residue', '' ])
      for a in self.modes:
         smm.append([ 1, self.mode_name[a], 'cmd.get_wizard().set_mode("'+a+'")'])
      self.menu['mode']=smm

   def set_mode(self,mode):
      if mode in self.modes:
         self.mode = mode
      cmd.refresh_wizard()
      
   def get_panel(self):
      return [
         [ 1, 'Mutagenesis',''],
         [ 3, self.mode_name[self.mode],'mode'],         
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
      auto_zoom = cmd.get_setting_text('auto_zoom')
      cmd.set('auto_zoom',"0",quiet=1)
      cmd.frame(0)
      cmd.delete(tmp_name)      
      if self.mode=="Current":
         pymol.stored.resn=""
         cmd.iterate("(%s and n;ca)"%sele_name,"stored.resn=resn")
         res_type = pymol.stored.resn
         cmd.create(tmp_name,sele_name,1,1)         
      else:
         res_type = self.mode
         cmd.fragment(string.lower(self.mode),tmp_name)
         cmd.remove("("+tmp_name+" and hydro)")
         # copy identifying information
         cmd.iterate("(%s and n;ca)"%sele_name,"stored.chain=chain")
         cmd.alter("(%s)"%tmp_name,"chain=stored.chain")
         cmd.iterate("(%s and n;ca)"%sele_name,"stored.resi=resi")
         cmd.alter("(%s)"%tmp_name,"resi=stored.resi")
         cmd.iterate("(%s and n;ca)"%sele_name,"stored.segi=segi")
         cmd.alter("(%s)"%tmp_name,"segi=stored.segi")
         # move the fragment
         cmd.pair_fit("(%s and n;ca)"%tmp_name,
                       "(%s and n;ca)"%sele_name,
                       "(%s and n;c)"%tmp_name,
                       "(%s and n;c)"%sele_name,
                       "(%s and n;n)"%tmp_name,
                       "(%s and n;n)"%sele_name)
         # fix the carbonyl position...
         cmd.iterate_state(1,"(%s and n;o)"%sele_name,"stored.list=[x,y,z]")
         cmd.alter_state(1,"(%s and n;o)"%tmp_name,"(x,y,z)=stored.list")
         
      if self.library.has_key(res_type):
         lib = self.library[res_type]
         state = 1
         for a in lib:
            cmd.create(obj_name,tmp_name,1,state)
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
            pymol.util.cbaw(obj_name)
#         cmd.delete(tmp_name)
         print " Mutagenesis: %d conformations loaded."%len(self.library[res_type])
         cmd.set('auto_zoom',auto_zoom,quiet=1)
         cmd.delete(tmp_name)
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



