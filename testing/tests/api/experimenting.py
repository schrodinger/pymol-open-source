import pymol
from pymol import cmd, testing, stored

class TestExperimenting(testing.PyMOLTestCase):

    def testGetBondPrint(self):
        cmd.get_bond_print
        self.skipTest("TODO")

    def testSpheroid(self):
        cmd.fragment('gly', 'm1')
        cmd.create('m1', 'm1', 1, 2)
        cmd.spheroid('m1')
        self.assertEqual(1, cmd.count_states('m1'))

    def testCheck(self):
        cmd.check
        self.skipTest("TODO")

    def testMinimize(self):
        cmd.minimize
        self.skipTest("TODO")

    def testFastMinimize(self):
        cmd.fast_minimize
        self.skipTest("TODO")

    def testDump(self):
        cmd.dump
        self.skipTest("TODO")

    def testImportCoords(self):
        cmd.import_coords
        self.skipTest("TODO")

    @testing.requires('incentive')
    def testFocalblur(self):
        cmd.viewport(100, 100)
        cmd.fragment('gly', 'm1')
        cmd.focal_blur(4.0, 3)
