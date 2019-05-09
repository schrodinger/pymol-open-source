# This wizard contributed by Ezequiel "Zac" Panepucci 011114
# modified by Warren L. DeLano
# modified by Thomas Holder

from pymol.wizard import Wizard
from pymol.exporting import _resn_to_aa as one_letter
from pymol import cmd

class Label(Wizard):

    atom=None
    messages=1
    mode = 0
    mode_names = [
        '{resn}-{resi}',
        '{onelettercode}{resi}',
        '{chain}/{resn}`{resi}',
        '{chain}/{resn}`{resi}/{name}`{alt}',
        '/{model}/{segi}/{chain}/{resn}`{resi}/{name}`{alt}',
        '{chain}',
        '{resn}',
        '{resi}',
        '{name}',
    ]

    field_names = ('model', 'segi', 'chain', 'resn', 'resi', 'name', 'alt',
            'b', 'x', 'y', 'z', 'label')
    fields_str_tuple = '(' + ','.join(field_names) + ')'

    def __init__(self,_self=cmd):
        Wizard.__init__(self,_self)
        self.cmd.unpick()
        self.cmd.deselect()
        self.cmd.set('mouse_selection_mode', 0)

        self.menu = {
            'mode': [
                [2, 'Mode', ''],
            ] + [
                [1, name, 'cmd.get_wizard().set_mode(%d)' % i]
                for (i, name) in enumerate(self.mode_names)
            ]
        }

    def get_event_mask(self):
        return Wizard.event_mask_pick + Wizard.event_mask_select

    def set_mode(self, i):
        self.mode = i
        self.cmd.refresh_wizard()

    def get_prompt(self):
        self.prompt = []
        if (not self.messages):
            return None

        if (self.atom is None):
            self.prompt = ['Click atoms...']
        else:
            self.prompt = [
                '/%s/%s/%s/%s`%s/%s`%s  B = %.2f  XYZ = %.3f %.3f %.3f' % self.atom[:-1]
            ]

        return self.prompt

    def toggle_messages(self):
        self.messages = not self.messages
        self.cmd.refresh_wizard()

    def get_panel(self):
        messages_str = 'On' if self.messages else 'Off'
        return [
            [ 1, 'Labeling',''],
            [ 3, 'Mode: ' + self.mode_names[self.mode], 'mode' ],
            [ 2, 'Messages: ' + messages_str, 'cmd.get_wizard().toggle_messages()'],
            [ 2, 'Done','cmd.set_wizard()'],
            ]

    def do_pick(self,bondFlag):
        self.do_select('pk1')
        self.cmd.unpick()

    def do_select(self, sele):
        self.atom = None
        n = self.cmd.iterate_state(-1, 'first ?' + sele,
                'self.atom = ' + self.fields_str_tuple,
                space=locals())

        if n == 0:
            return

        if self.atom[-1]:
            label = ''
        else:
            d = dict(zip(self.field_names, self.atom))
            d['onelettercode'] = one_letter.get(d['resn'], d['resn'])
            label = repr(self.mode_names[self.mode].format(**d))

        self.cmd.label('first ?' + sele, label)
        self.cmd.deselect()
        self.cmd.refresh_wizard()
