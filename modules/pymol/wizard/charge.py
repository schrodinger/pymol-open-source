
from pymol.wizard import Wizard
from pymol import cmd
import pymol

class Charge(Wizard):

   def __init__(self):

      Wizard.__init__(self)
      
      self.modes = [
         'labchg',
         'cpychg',
         'addchg',
         'zrochg',
         'sumchg',
         ]

      self.mode = 'labchg'
      self.status = 0
      
      self.mode_name = {
         'labchg':'Show Charges',
         'cpychg':'Copy Charges',
         'addchg':'Add Charges',
         'zrochg':'Zero Charges',
         'sumchg':'Get Total Charge',
         }
      
      # initialize mode menu
      
      smm = []
      smm.append([ 2, 'Atom Charge Mode', '' ])
      for a in self.modes:
         smm.append([ 1, self.mode_name[a], 'cmd.get_wizard().set_mode("'+a+'")'])

      self.menu['mode']=smm
      
      self.memory = 0
         

   def get_panel(self):
      return [
         [ 1, 'Charge Wizard',''],
         [ 3, self.mode_name[self.mode],'mode'],
         [ 2, 'Clear','cmd.get_wizard().clear()'],
         [ 2, 'Done','cmd.set_wizard()'],
         ]

   def cleanup(self):
      self.clear()
      
   def clear(self):
      self.set_status(0)
      lst = cmd.get_names('selections')
      if '_charge' in lst:
         cmd.edit("_charge")
         cmd.label("pkchain",'') # fastest clear command
         cmd.delete("_charge")
         cmd.unpick()
      cmd.unpick()
      
   def get_prompt(self):
      self.prompt = None
      if self.mode == 'cpychg':
         if self.status==0:
            self.prompt = [ 'Pick source atom...' ]
         elif self.status==1:
            self.prompt = [ 'Pick destination atom on which to assign charge %6.4f'%self.partial_charge ]
      if self.mode == 'addchg':
         if self.status==0:
            self.prompt = [ 'Pick source atom...' ]
         elif self.status==1:
            self.prompt = [ 'Pick destination atom on which to add charge %6.4f'%self.partial_charge ]
      if self.mode == 'sumchg':
         if self.status==0:
            self.prompt = [ 'Pick an atom on the chain...' ]
         if self.status==1:
            self.prompt = [ 'Total charge on the chain is %6.4f'%self.partial_charge,
                            'Pick an atom on the chain...' ]
            
      if self.mode == 'zrochg':
            self.prompt = [ 'Pick atom on which to zero charge...' ]


      if self.mode == 'labchg':
            self.prompt = [ 'Pick atom on which to show charge...' ]

         
      return self.prompt
   
   def set_mode(self,mode):
      if mode in self.modes:
         self.mode = mode
      self.status = 0
      cmd.refresh_wizard()
      
   def set_status(self,status):
      self.status = status
      cmd.refresh_wizard()

   def do_pick(self,bondFlag):
      if bondFlag:
         print " Error: please select a single atom"
         
      if self.mode == 'cpychg':
         # picking up
         if self.status==0:
            if cmd.iterate("(pk1)","stored.charge = partial_charge"):
               self.partial_charge = pymol.stored.charge
               self.status = 1
               cmd.label("(pk1)","'%6.4f'%partial_charge")
               cmd.select("_charge","(pk1)")
               cmd.unpick()
               cmd.enable("_charge")
               
         # dropping off
         elif self.status==1:
            pymol.stored.charge=self.partial_charge
            if cmd.alter("(pk1)","partial_charge = stored.charge"):
               self.status = 0
               cmd.label("(pk1)","'%6.4f'%partial_charge")
               cmd.select("_charge","(pk1)")
               cmd.unpick()
               cmd.enable("_charge")

      if self.mode == 'addchg':
         # picking up
         if self.status==0:
            if cmd.iterate("(pk1)","stored.charge = partial_charge"):
               self.partial_charge = pymol.stored.charge
               self.status = 1
               cmd.label("(pk1)","'%6.4f'%partial_charge")
               cmd.select("_charge","(pk1)")
               cmd.unpick()
               cmd.enable("_charge")
               
         # dropping off
         elif self.status==1:
            pymol.stored.charge=self.partial_charge
            if cmd.alter("(pk1)","partial_charge = partial_charge + stored.charge"):
               self.status = 0
               cmd.label("(pk1)","'%6.4f'%partial_charge")
               cmd.select("_charge","(pk1)")
               cmd.unpick()
               cmd.enable("_charge")

      if self.mode == 'zrochg':
         if cmd.alter("(pk1)","partial_charge = 0.0"):
            cmd.label("(pk1)","'%6.4f'%partial_charge")
            cmd.select("_charge","(pk1)")
            cmd.unpick()
            cmd.enable("_charge")
               
      if self.mode == 'labchg':
         cmd.label("(pk1)","'%6.4f'%partial_charge")
         cmd.select("_charge","(pk1)")
         cmd.unpick()
         cmd.enable("_charge")
               
      if self.mode == 'sumchg':
         pymol.stored.charge = 0.0
         if cmd.iterate("(pkchain)","stored.charge = stored.charge + partial_charge"):
            self.partial_charge = pymol.stored.charge
            self.status = 1
            cmd.select("_charge","(pkchain)")
            cmd.unpick()
            cmd.enable("_charge")
         
      cmd.refresh_wizard()
