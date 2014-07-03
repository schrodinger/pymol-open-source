'''
PYMOL-317
can't select by atom type
'''

import os
from pymol import cmd, testing, stored

@testing.requires('incentive')
class TestPYMOL1567(testing.PyMOLTestCase):

    @testing.foreach(1, 0)
    def testUndoAfterRemoveAtomOnDiscrete(self, discr):
        cmd.load(self.datafile('ligs3d.sdf'), discrete=discr)
        natms = cmd.count_atoms()
        cmd.remove("index 1")
        cmd.undo()
        self.assertEqual(natms, cmd.count_atoms())

