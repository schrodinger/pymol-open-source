'''
PYMOL-317
can't select by atom type
'''

import os
from pymol import cmd, testing, stored
import unittest

@testing.requires('mmlibs')
@testing.requires('incentive')
@testing.requires('no_edu')
class TestPYMOL317(testing.PyMOLTestCase):

    def _testLabelByAtomType(self):
        cmd.load(self.datafile('1oky.pdb.gz'))
        cmd.remove("alt 'B'")
        cmd.label('resn STU', 'text_type')
        c3labels = []
        cmd.iterate("resn STU", "c3labels.append(label == 'C.3')", space=locals())
        self.assertEquals(sum(c3labels), 8)

    def testSelectByAtomType(self):
        cmd.load(self.datafile('1oky.pdb.gz'))
        cmd.remove("alt 'B'")
        cmd.select('foo', "resn STU and text_type 'C.3'")
        self.assertEquals(cmd.count_atoms('foo'), 8)

        # discrete object
        cmd.delete('*')
        cmd.load(self.datafile('ligs3d.sdf'))
        self.assertTrue(cmd.count_discrete('*') > 0)
        n = cmd.select('foo', "text_type 'C.3'")
        self.assertEquals(n, 115)

        # only run label test if selection tests passed. Don't wnat
        # more than one failing tests for builds with NO_MMLIBS
        cmd.delete('*')
        self._testLabelByAtomType()
