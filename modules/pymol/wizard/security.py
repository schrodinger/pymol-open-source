# This wizard contributed by Ezequiel "Zac" Panepucci 011114
# modified by Warren L. DeLano

from pymol.wizard import Wizard
from pymol import cmd
import pymol

class Security(Wizard):

   def __init__(self):
      for a in self.get_prompt():
         print a
         
   def get_prompt(self):
      self.prompt = [ '=== PyMOL SECURITY WARNING ===',
                      '',
                      'Do you trust the person who created the file you just loaded?',
                      '',
                      'It contains movie commands which could be used maliciously',
                      'to control your computer or damage your files.',
                      '',
                      'Click or enter "accept" to assume the risks of running the commands.',
                      '',
                      'Click or enter "decline" to decline the risks and delete the commands.',
                      '',
                      'Click or enter "mdump" to print out the commands in the movie.',
                      '',
                      'To avoid this message in the future, "set security,off"',
                      'or launch pymol with the -o option.',
                      ]
      return self.prompt


   def get_panel(self):
      return [
         [ 1, 'Assume Movie Risks?',''],
         [ 2, 'accept','cmd.accept()'],
         [ 1, '',''],
         [ 2, 'decline','cmd.decline()'],
         [ 1, '',''],
         [ 2, 'mdump','cmd.mdump()'],
         ]


