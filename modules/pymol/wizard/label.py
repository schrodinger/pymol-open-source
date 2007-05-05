# This wizard contributed by Ezequiel "Zac" Panepucci 011114
# modified by Warren L. DeLano

from pymol.wizard import Wizard
from pymol import cmd
import pymol

class Label(Wizard):

    atom=None
    messages=1
    labeling=1
    obj_name=None

    def __init__(self,_self=cmd):
        Wizard.__init__(self,_self)
        self.cmd.unpick()

        
    def get_prompt(self):
        self.prompt = []
        if (not self.messages):
            return None

        if (self.atom == None):
            self.prompt = ['Click atoms...']
        else:
            if self.atom.chain == '':
                self.prompt.append( '%s %s %s %s B = %.2f  XYZ = %.3f %.3f %.3f' %
                                     (self.obj_name,
                                      self.atom.resn,
                                      self.atom.resi,
                                      self.atom.name,
                                      self.atom.b,
                                      self.atom.coord[0],
                                      self.atom.coord[1],
                                      self.atom.coord[2]) )
            else:
                self.prompt.append('%s %s %s%s %s B = %.2f  XYZ = %.3f %.3f %.3f' %
                                     (self.obj_name,
                                      self.atom.resn,
                                      self.atom.chain,
                                      self.atom.resi,
                                      self.atom.name,
                                      self.atom.b,
                                      self.atom.coord[0],
                                      self.atom.coord[1],
                                      self.atom.coord[2]) )

        return self.prompt

    def toggle_messages(self):
        self.messages = not self.messages

    def toggle_labeling(self):
        self.labeling = not self.labeling

    def get_panel(self):
        return [
            [ 1, 'Labeling',''],
            [ 2, 'Toggle add/erase','cmd.get_wizard().toggle_labeling()'],
            [ 2, 'Toggle messages','cmd.get_wizard().toggle_messages()'],
            [ 2, 'Clear All','cmd.label()'],
            [ 2, 'Done','cmd.set_wizard()'],
            ]

    def do_pick(self,bondFlag):
        self.obj_name = None

#      if 'pk1' in cmd.get_names('selections'):
        if cmd.count_atoms('pk1',1):
            self.obj_name = cmd.identify('pk1',1)[0][0]

        model = cmd.get_model("(pk1)")
        self.atom = model.atom.pop()
        if not self.labeling:
            cmd.label("(pk1)", '""')
        elif self.atom.name == 'CA':
            cmd.label("(pk1)", '" %s %s" % (resn,resi)')
        else:
            cmd.label("(pk1)", '" %s %s" % (name,resi)')
        cmd.unpick()
        cmd.refresh_wizard()


