
from pymol.wizard import Wizard
from pymol import cmd
from chempy import io

import pymol
import os
import string

sele_name = "_mutate_sel"

obj_name = "mutation"
tmp_name = "_tmp_mut"
tmp_obj1 = "_tmp_obj1"
tmp_obj2 = "_tmp_obj2"
tmp_obj3 = "_tmp_obj3"

default_mode = "current"
default_rep = "lines"

class Mutagenesis(Wizard):

   count = 0
   cutoff = 3.5
   
   def __init__(self):

      cmd.unpick()
      
      Wizard.__init__(self)

      self.library = io.pkl.fromFile(os.environ['PYMOL_PATH']+
                                     "/data/chempy/sidechains/sc_library.pkl")
      
      self.status = 0 # 0 no selection, 1 mutagenizing
      self.error = None
      self.object_name = None
      self.modes = [
         'current'
         ]
      self.mode = default_mode
      self.rep = default_rep
      residues = self.library.keys()
      residues.extend(['GLY','PRO','ALA'])
      residues.sort()
      self.modes.extend(residues)
      self.mode_name={}
      for a in self.modes:
         self.mode_name[a] = "-> "+a
      self.mode_name['current']="No Mutation"

      smm = []
      smm.append([ 2, 'Substitution', '' ])
      for a in self.modes:
         smm.append([ 1, self.mode_name[a], 'cmd.get_wizard().set_mode("'+a+'")'])
      self.menu['mode']=smm

      self.reps = [
         'lines',
         'sticks',
         'spheres',
         'dots'
         ]

      self.rep_name = {
         'lines' : "Show Lines",
         'sticks' : "Show Sticks",
         'spheres' : "Show Spheres",
         'dots' : "Show Dots",
         }

      smm = []
      smm.append([ 2, 'Representation', '' ])
      for a in self.reps:
         smm.append([ 1, self.rep_name[a], 'cmd.get_wizard().set_rep("'+a+'")'])
      self.menu['rep']=smm

      if 'pk1' in cmd.get_names('selections'):
         cmd.select(sele_name,"(byres pk1)")
         cmd.unpick()
         cmd.enable(sele_name)
         self.status = 1
         self.error = None
         self.do_library()
         cmd.refresh_wizard()

   def set_mode(self,mode):
      if mode in self.modes:
         self.mode = mode
      if self.status==1:
         self.do_library()
      cmd.refresh_wizard()
      
   def set_rep(self,rep):
      if rep in self.reps:
         self.rep=rep
      cmd.hide("("+obj_name+")")
      cmd.show('lines',obj_name) # always show lines      
      cmd.show(self.rep,obj_name)
      cmd.refresh_wizard()
      
   def get_panel(self):
      return [
         [ 1, 'Mutagenesis',''],
         [ 3, self.mode_name[self.mode],'mode'],
         [ 3, self.rep_name[self.rep],'rep'],                  
         [ 2, 'Apply' , 'cmd.get_wizard().apply()'],         
         [ 2, 'Clear' , 'cmd.get_wizard().clear()'],
         [ 2, 'Done','cmd.set_wizard()'],
         ]

   def cleanup(self):
      global default_mode,default_rep
      default_mode = self.mode
      default_rep = self.rep
      self.clear()
      
   def clear(self):
      self.status=0
      cmd.delete(tmp_obj1)
      cmd.delete(tmp_obj2)
      cmd.delete(tmp_obj3)
      cmd.delete(sele_name)
      cmd.delete(obj_name)
      cmd.refresh_wizard()
      
   def apply(self):
      if self.status==1:
         # find the name of the object which contains the selection
         new_name = None
         obj_list = cmd.get_names('objects')
         for a in obj_list:
            if cmd.get_type(a)=="object:molecule":
               if cmd.count_atoms("(%s and %s)"%(a,sele_name)):
                  new_name = a
                  break
         src_frame = cmd.get_state()
         if new_name==None:
            print " Mutagenesis: object not found."
         else:
            auto_zoom = cmd.get_setting_text('auto_zoom')
            cmd.set('auto_zoom',"0",quiet=1)
            if self.mode!="current":
               # create copy w/o residue
               cmd.create(tmp_obj1,"(%s and not %s)"%(new_name,sele_name))
               # save copy for bonded atom reference
               cmd.create(tmp_obj3,new_name)
               # transfer the selection to copy
               cmd.select(sele_name,"(%s in %s)"%(tmp_obj3,sele_name))
               # create copy with mutant in correct frame
               cmd.create(tmp_obj2,obj_name,src_frame,1)
               cmd.delete(new_name)
               # create the merged molecule
               cmd.create(new_name,"(%s or %s)"%(tmp_obj1,tmp_obj2),1) # only one state in merged object...
               # now connect them
               cmd.bond("(%s in %s and n;N)"%(new_name,tmp_obj2),
                        "(name C and (%s in (neighbor %s)))"%
                        (new_name,sele_name))
               cmd.bond("(%s in %s and n;C)"%(new_name,tmp_obj2),
                        "(name N and (%s in (neighbor %s)))"%
                        (new_name,sele_name))
               # now transfer selection back to the modified object
               cmd.delete(tmp_obj1)
               cmd.delete(tmp_obj2)
               cmd.delete(tmp_obj3)
               self.clear()
               # and return to frame 1
               cmd.frame(1)
               cmd.refresh_wizard()               
            else:
               # create copy with conformation in correct state
               cmd.create(tmp_obj2,obj_name,src_frame,1)
               # save existing conformation on undo stack
#               cmd.edit("((%s in %s) and name ca)"%(new_name,sele_name))
               cmd.push_undo("("+sele_name+")")
               # modify the conformation
               cmd.update(new_name,tmp_obj2)
#               cmd.unpick()
               cmd.delete(tmp_obj2)
               self.clear()
               # and return to frame 1
               cmd.frame(1)
               cmd.refresh_wizard()                              
            cmd.set('auto_zoom',auto_zoom,quiet=1)
               
   def get_prompt(self):
      self.prompt = None
      if self.status==0:
         self.prompt = [ 'Pick a residue...']
      elif self.status==1:
         self.prompt = [ 'Select a conformational state, or pick a new residue...' ]
      return self.prompt

   def do_library(self):
      cmd.feedback("push")
      cmd.feedback("disable","selector","everythin")
      cmd.feedback("disable","editor","actions")

      self.prompt = [ 'Loading rotamers...']
      
      auto_zoom = cmd.get_setting_text('auto_zoom')
      cmd.set('auto_zoom',"0",quiet=1)
      cmd.frame(0)
      cmd.delete(tmp_name)
      if self.mode=="current":
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
         cmd.iterate("(%s and n;ca)"%sele_name,"stored.ss=ss")
         cmd.alter("(%s)"%tmp_name,"ss=stored.ss")
         # move the fragment
         if ((cmd.count_atoms("(%s and n;cb)"%tmp_name)>0) and
             (cmd.count_atoms("(%s and n;cb)"%sele_name)>0)):
            cmd.pair_fit("(%s and n;ca)"%tmp_name,
                         "(%s and n;ca)"%sele_name,
                         "(%s and n;cb)"%tmp_name,
                         "(%s and n;cb)"%sele_name,
                         "(%s and n;c)"%tmp_name,
                         "(%s and n;c)"%sele_name,
                         "(%s and n;n)"%tmp_name,
                         "(%s and n;n)"%sele_name)
         else:
            cmd.pair_fit("(%s and n;ca)"%tmp_name,
                         "(%s and n;ca)"%sele_name,
                         "(%s and n;c)"%tmp_name,
                         "(%s and n;c)"%sele_name,
                         "(%s and n;n)"%tmp_name,
                         "(%s and n;n)"%sele_name)
         # fix the carbonyl position...
         cmd.iterate_state(1,"(%s and n;o)"%sele_name,"stored.list=[x,y,z]")
         cmd.alter_state(1,"(%s and n;o)"%tmp_name,"(x,y,z)=stored.list")

      cartoon = (cmd.count_atoms("(%s and n;ca and rep cartoon)"%sele_name)>0)
      sticks = (cmd.count_atoms("(%s and n;ca and rep sticks)"%sele_name)>0)
         
      cmd.delete(obj_name)
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
                  cmd.set_title(obj_name,state,"%1.1f%%"%(a[b]*100))
            state = state + 1
         cmd.delete(tmp_name)
         cmd.set("seq_view",0,obj_name,quiet=1)

         print " Mutagenesis: %d conformations loaded."%len(self.library[res_type])
      else:
         cmd.create(obj_name,tmp_name,1,1)
         print " Mutagenesis: no additional conformations in library."
      pymol.util.cbaw(obj_name)
      cmd.hide("("+obj_name+")")
      cmd.show(self.rep,obj_name)
      cmd.show('lines',obj_name) # always show lines
      if cartoon:
         cmd.show("cartoon",obj_name)
      if sticks:
         cmd.show("sticks",obj_name)
      cmd.set('auto_zoom',auto_zoom,quiet=1)
      cmd.delete(tmp_name)
      cmd.frame(0)
      cmd.unpick()
      cmd.feedback("pop")

   def do_select(self,selection):
      if self.status!=0:
         cmd.delete(obj_name)
      cmd.select(sele_name,"(byres first (%s))"%selection)
      cmd.unpick()
      cmd.enable(sele_name)
      self.status = 1
      self.error = None
      self.do_library()
      cmd.refresh_wizard()
      
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



