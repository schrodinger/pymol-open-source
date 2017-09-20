'''
Object visibility with hidden parent in scenes
'''

import os
from pymol import cmd, CmdException, testing, stored

@testing.requires_version('1.8.7')
class TestSceneObjVis(testing.PyMOLTestCase):
    def test(self):
        cmd.pseudoatom('m1')
        cmd.pseudoatom('m2')
        cmd.pseudoatom('m3')
        cmd.group('g1', 'm1')
        cmd.disable('g1')
        cmd.disable('m2')
        cmd.ray(1, 1)  # force scene update
        cmd.scene('s1', 'store')

        self._test_recall()
        cmd.disable('*')
        self._test_recall()
        cmd.enable('*')
        self._test_recall()

    def _test_recall(self):
        cmd.scene('s1', 'recall')
        cmd.ray(1, 1)  # force scene update
        self.assertEqual(['m3'], cmd.get_object_list('(visible)'))
        self.assertEqual(['m3'], cmd.get_object_list('(enabled)'))
        cmd.enable('g1')
        cmd.ray(1, 1)  # force scene update
        self.assertEqual(['m1', 'm3'], cmd.get_object_list('(visible)'))
        self.assertEqual(['m1', 'm3'], cmd.get_object_list('(enabled)'))
