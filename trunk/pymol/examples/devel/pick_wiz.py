
# demonstrate how to use PyMOL's atom pick "events" in a Wizard

# Run this file as:
#  DOS/Unix> pymol pick_wiz.py
#     or
#  PyMOL> run pick_wiz.py

from pymol.wizard import Wizard
from pymol import cmd
import pymol

class PickWizard(Wizard):

   def reset(self):
      self.pk1_st = None
      self.pk2_st = None
      self.pk1_xyz = None
      cmd.refresh_wizard()

   def __init__(self):
      Wizard.__init__(self)
      self.reset()
      
   def get_prompt(self):

      if self.pk2_st!=None:
         return ["You picked the bond between %s and %s"%(
            self.pk1_st, self.pk2_st)]
      elif self.pk1_st!=None:
         return ["You picked atom %s"%(self.pk1_st),
                 "At X=%1.2f Y=%1.2f Z=%1.2f"%self.pk1_xyz]
      else:
         return ["Please pick an atom or a bond..."]

   def do_pick(self,picked_bond):

      self.reset()

      cmd.iterate("pk1","setattr(cmd.get_wizard(),'pk1_st',"
                  "'%s/%s/%s/%s/%s'%(model,segi,chain,resi,name))")
      if picked_bond:
         cmd.iterate("pk1","setattr(cmd.get_wizard(),'pk2_st',"
                     "'%s/%s/%s/%s/%s'%(model,segi,chain,resi,name))")
      else:

         # for single atom, also get 3D coordinates (EXAMPLE)
         
         cmd.iterate_state( cmd.get_state(),
                           "pk1","setattr(cmd.get_wizard(),'pk1_xyz',(x,y,z))")
         
      cmd.unpick() 
      cmd.refresh_wizard()

   def get_panel(self):
      return [
         [ 1, 'Example Wizard',''],         
         [ 2, 'Reset','cmd.get_wizard().reset()'],
         [ 2, 'Done','cmd.set_wizard()'],
         ]

# create an instane

wiz = PickWizard()

# make this the active wizard

cmd.set_wizard(wiz)



