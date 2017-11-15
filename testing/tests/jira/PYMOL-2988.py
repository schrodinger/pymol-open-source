'''
some elements assigned wrong from GRO files
'''

import os
from pymol import cmd, CmdException, testing, stored

@testing.requires_version('1.7.4')
class TestElem(testing.PyMOLTestCase):

    def test(self):
        cmd.set('retain_order')
        cmd.load(self.datafile('MCRKY-Mg-Ca-Cl.gro'))

        elem_set = set()
        cmd.iterate('index 1-98', 'elem_set.add(elem)', space=locals())
        self.assertEqual(elem_set, set(['C', 'H', 'N', 'O', 'S']))

        elem_set = set()
        cmd.iterate('index 99-101', 'elem_set.add(elem)', space=locals())
        self.assertEqual(elem_set, set(['Ca', 'Cl', 'Mg']))
