
from pymol.wizard import Wizard
from pymol import cmd
import pymol

sele_prefix = "_dw"
sele_prefix_len = len(sele_prefix)

dist_prefix = "wdist"

class Distance(Wizard):

   count = 0
   cutoff = 3.5
   
   def __init__(self):

      Wizard.__init__(self)
      
      self.status = 0 # 0 no atoms selections, 1 atom selected
      self.error = None
      self.object_name = None

      # mode selection subsystem
      
      self.mode = 'pairs'
      self.modes = [
         'polar',
         'heavy',
         'neigh',
         'pairs',
         ]
      self.mode_name = {
         'polar':'Polar Neighbors',
         'heavy':'Heavy Neighbors',
         'neigh':'Neighbors',
         'pairs':'Pairwise Distances',
         }

      smm = []
      smm.append([ 2, 'Measurement Mode', '' ])
      for a in self.modes:
         smm.append([ 1, self.mode_name[a], 'cmd.get_wizard().set_mode("'+a+'")'])
      self.menu['mode']=smm

      # overwrite mode selection subsystem
      
      self.object_mode='overwr'
      self.object_modes = [
         'overwr',
         'append',
         ]
      self.object_mode_name = {
         'overwr':'Overwrite',
         'append':'Append',         
         }

      smm = []
      smm.append([ 2, 'New Distances?', '' ])
      for a in self.object_modes:
         smm.append([ 1, self.object_mode_name[a], 'cmd.get_wizard().set_object_mode("'+a+'")'])
      self.menu['object_mode']=smm

# generic set routines

   def set_mode(self,mode):
      if mode in self.modes:
         self.mode = mode
      self.status = 0
      self.clear()
      cmd.refresh_wizard()

   def set_object_mode(self,mode):
      if mode in self.object_modes:
         self.object_mode = mode
      self.status = 0
      cmd.refresh_wizard()

      
   def get_panel(self):
      return [
         [ 1, 'Distance Measurement',''],
         [ 3, self.mode_name[self.mode],'mode'],
         [ 3, self.object_mode_name[self.object_mode],'object_mode'],
         [ 2, 'Delete Last' , 'cmd.get_wizard().delete_last()'],
         [ 2, 'Done','cmd.set_wizard()'],
         ]

   def cleanup(self):
      self.clear()
      
   def clear(self):
      cmd.delete(sele_prefix+"*")
      
   def get_prompt(self):
      self.prompt = None
      if self.mode == 'pairs':
         if self.status==0:
            self.prompt = [ 'Pick the first atom...']
         elif self.status==1:
            self.prompt = [ 'Pick the second atom...' ]
      elif self.mode in [ 'polar', 'neigh' ]:
         self.prompt = [ 'Pick an atom...']
      if self.error!=None:
         self.prompt.append(self.error)
      return self.prompt
   
   def delete_last(self):
      if self.status==0:
         if self.__class__.count>0:
            cmd.delete(dist_prefix+"%02d"%self.__class__.count)
            self.__class__.count = self.__class__.count - 1
      self.status=0
      self.error = None
      self.clear()
      cmd.refresh_wizard()
      
   def do_pick(self,bondFlag):
      if bondFlag:
         self.error = "Error: please select an atom, not a bond."
         print self.error
      else:
         if self.mode == 'pairs':
            if self.status==0:
               name = sele_prefix 
               cmd.select(name,"(pk1)")
               cmd.unpick()
               cmd.enable(name)
               self.status = 1
               self.error = None
            elif self.status==1:
               if ((self.object_mode=='append') or (not self.__class__.count)):
                  self.__class__.count = self.__class__.count + 1
               else:
                  cmd.delete(dist_prefix+"%2d"%self.__class__.count)
               name = dist_prefix + "%02d"%self.__class__.count
               cmd.dist(name,sele_prefix,"(pk1)")
               cmd.delete(sele_prefix)
               cmd.unpick()
               cmd.enable(name)
               self.status = 0
         elif self.mode in ['neigh','polar','heavy']:
            if ((self.object_mode=='append') or (not self.__class__.count)):
               self.__class__.count = self.__class__.count + 1
            else:
               cmd.delete(dist_prefix+"%2d"%self.__class__.count)
            name = dist_prefix + "%02d"%self.__class__.count
            if self.mode == 'neigh':
               cmd.dist(name,"(pk1)","((pk1 a; %f) and (not (neighbor pk1)) and (not (neighbor (neighbor pk1))) and (not (neighbor (neighbor (neighbor pk1)))))"%self.__class__.cutoff)
            elif self.mode == 'polar':
               cmd.dist(name,"(pk1)","((pk1 a; %f) and (e;n,o) and (not (neighbor pk1)) and (not (neighbor (neighbor pk1))) and (not (neighbor (neighbor (neighbor pk1)))))"%self.__class__.cutoff)            
            elif self.mode == 'heavy':
               cmd.dist(name,"(pk1)","((pk1 a; %f) and (not hydro) and (not (neighbor pk1)) and (not (neighbor (neighbor pk1))) and (not (neighbor (neighbor (neighbor pk1)))))"%self.__class__.cutoff)            
            cmd.unpick()
            cmd.enable(name)
      cmd.refresh_wizard()
