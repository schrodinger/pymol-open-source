
# demonstrate how to capture PyMOL's atom selection events

from pymol.wizard import Wizard
from pymol import cmd
import pymol

sele_name = "lb" # must be set to "lb" for now...

class PickWizard(Wizard):

   def set_buttons(self):
      
      # just use selections (disable atom and bond picking)
      
      cmd.button('m','ctrl','+lb')
      cmd.button('r','ctrl','none')
      cmd.button('r','ctsh','none')
      
   def get_prompt(self):

      # returns prompt for the viewer window (optional)
      
      if sele_name in cmd.get_names('selections'):
         n_atom = cmd.count_atoms(sele_name)
      else:
         n_atom = 0
      if n_atom:
         list = cmd.identify(sele_name)
         return ["%d atoms selected..."%n_atom,str(list)]
      else:
         return ["Please select some atoms..."]

   def do_select(self,name):

      # handle mouse selection callback
      
      if not sele_name in cmd.get_names('selections'):
         cmd.select(sele_name,'none')
      cmd.enable(sele_name)
      cmd.refresh_wizard()

   def get_panel(self):
      return [
         [ 1, 'Example Wizard',''],         
         [ 2, 'Clear Selection',
           'cmd.delete("'+sele_name+'");cmd.refresh_wizard()'],
         ]

wiz = PickWizard()

# make this the active wizard

cmd.set_wizard(wiz)

# reconfigure the mouse panel

wiz.set_buttons()

