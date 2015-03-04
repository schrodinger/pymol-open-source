'''
undo atom type change
'''

import unittest
from pymol import cmd, testing, stored

@testing.requires('incentive')
@testing.requires('no_edu')
class Test1276(testing.PyMOLTestCase):

    def test(self):
        cmd.fragment('gly', 'm1')
        cmd.create('m2', 'm1')
        cmd.edit('m1 & name N')
        cmd.replace('I', 1, 1)
        self.assertEqual(cmd.count_atoms('m1 in m2'), 5)
        cmd.undo()
        self.assertEqual(cmd.count_atoms('m1 in m2'), 7)
