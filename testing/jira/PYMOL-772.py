'''
Regression test for PYMOL-772
partial load of PSE file duplicates groups with identical name 
'''

import unittest
from pymol import cmd, testing

class Test772(testing.PyMOLTestCase):

    v_ref = ['m1', 'm2', 'g1']

    def _load_session_twice(self):
        cmd.load('PYMOL-772-example.pse.gz', partial=1)
        cmd.load('PYMOL-772-example.pse.gz', partial=1)

    def testDoNotDuplicateGroup(self):
        cmd.pseudoatom('m3')
        cmd.group('g1', 'm3')
        cmd.load('PYMOL-772-example.pse.gz', partial=1)

        v = cmd.get_names()
        self.assertEqual(v, ['m3', 'g1', 'm1', 'm2'])

    @unittest.skip("alter not implemented yet")
    def testAppend(self):
        cmd.set('auto_rename_duplicate_objects', 0)

        self._load_session_twice()

        v = cmd.get_names()
        self.assertEqual(v, self.v_ref)
 
    def testRename(self):
        cmd.set('auto_rename_duplicate_objects', 1)

        self._load_session_twice()

        v = cmd.get_names()
        v_ref = self.v_ref + [n + '_2' for n in self.v_ref
                if not n.startswith('g')]
        self.assertEqual(v, v_ref)

# vi:nowrap
