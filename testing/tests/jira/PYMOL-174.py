'''
RigiMOL fails with fewer than 4 atoms
'''

import unittest
from pymol import cmd, testing, stored

@testing.requires('incentive')
class Test174(testing.PyMOLTestCase):
    def test(self):
        cmd.set('suspend_undo')
        cmd.fragment('gly', 'm1')
        cmd.remove('not ID 0+1')
        cmd.create('m1', 'm1', 1, 2)
        cmd.rotate('x', 90, 'm1', 2)
        steps = 5
        cmd.morph('mout', 'm1', refinement=0, steps=steps, method='rigimol')
        self.assertEqual(steps, cmd.count_states('mout'))
        self.assertEqual(cmd.count_atoms('m1'),
                         cmd.count_atoms('mout'))
