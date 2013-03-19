'''
PYMOL-317
can't select by atom type
'''

import os
from pymol import cmd, testing, stored

@testing.requires('incentive')
class TestPYMOL317(testing.PyMOLTestCase):

    def testSelectByAtomType(self):
        cmd.load(self.datafile('1oky.pdb'))
        cmd.remove("alt 'B'")
        cmd.label('resn STU', 'text_type')
        cmd.select('foo', "resn STU and text_type 'C.3'")
        self.assertEquals(cmd.count_atoms('foo'), 8)

    def testSelectByAtomTypeNoLabel(self):
        cmd.load(self.datafile('1oky.pdb'))
        cmd.remove("alt 'B'")
        cmd.select('foo', "resn STU and text_type 'C.3'")
        self.assertEquals(cmd.count_atoms('foo'), 8)

    def _testSelectByAtomTypeNoLabelDiscrete(self):
        cmd.load(self.datafile('ligs3d.sdf'))
        cmd.select('foo', "text_type 'C.3'")
#        self.assertEquals(cmd.count_atoms('foo'), 8)

