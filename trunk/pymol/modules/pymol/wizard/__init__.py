
class Wizard:

   def __init__(self):
      self.menu = {}
      self.prompt = None
      self.panel = None
      
   def get_prompt(self):
      return self.prompt

   def get_panel(self):
      return self.panel

   def do_pick(self,bondFlag):
      pass

   def cleanup(self):
      pass

   def get_menu(self,tag):
      result = None
      if self.menu.has_key(tag):
         result = self.menu[tag]
      return result

