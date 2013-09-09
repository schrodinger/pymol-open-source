'''
Testing selection by properties
'''

import unittest
from pymol import cmd, testing

@testing.requires('properties')
class TestProperties(testing.PyMOLTestCase):
    def test(self):
        cmd.fragment('ala')
        cmd.alter_state(1, 'all', '(p.x_lt_0, p.y_gt_0) = (x < 0, y > 0)')
        cmd.select('sele_p_x', 'p. x_lt_0')
        cmd.select('sele_p_y', 'p. y_gt_0')
        cmd.select('sele_x', 'x < 0.0')
        cmd.select('sele_y', 'y > 0.0')
        counts = [
            cmd.count_atoms('sele_p_x'),
            cmd.count_atoms('sele_x'),
            cmd.count_atoms('sele_x & sele_p_x'),
        ]
        self.assertEqual(counts[0], 8)
        self.assertEqual(counts[0], counts[1])
        self.assertEqual(counts[0], counts[2])
        counts = [
            cmd.count_atoms('sele_p_y'),
            cmd.count_atoms('sele_y'),
            cmd.count_atoms('sele_y & sele_p_y'),
        ]
        self.assertEqual(counts[0], 6)
        self.assertEqual(counts[0], counts[1])
        self.assertEqual(counts[0], counts[2])
