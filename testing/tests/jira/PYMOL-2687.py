'''
PYMOL-2687 byring selection operator
'''

from pymol import cmd, CmdException, testing, stored

class TestPYMOL2687(testing.PyMOLTestCase):

    @testing.requires_version('1.8.1.0')
    def testByRing(self):
        cmd.fragment('trp')
        self.assertEqual(cmd.count_atoms('byring name CB'), 0)
        self.assertEqual(cmd.count_atoms('byring name CG'), 5)
        self.assertEqual(cmd.count_atoms('byring name CD2'), 9)
        self.assertEqual(cmd.count_atoms('byring name CE3'), 6)
        self.assertEqual(cmd.count_atoms('byring all'), 9)
