from pymol.wizard import Wizard
from pymol import cmd
import pymol
import types

class Message(Wizard):

   def __init__(self,*arg):
      self.message = []
      for a in arg:
         if not isinstance(a,types.ListType):
            self.message.append(a)
         else:
            self.message.extend(a)
      for a in self.message:
         print a
         
   def get_prompt(self):
      self.prompt = self.message
      return self.prompt


   def get_panel(self):
      return [
         [ 1, 'Message', '' ],
         [ 2, 'Dismiss', 'cmd.set_wizard()' ]
         ]



