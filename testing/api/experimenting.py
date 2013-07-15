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

    def testLoadCoords(self):
        cmd.fragment('gly', 'm1')
        coords = cmd.get_model('m1').get_coord_list()
        # buggy state argument
        cmd.load_coords(coords, 'm1', state=2+1)
        self.assertEqual(2, cmd.count_states('m1'))

    @testing.requires('incentive')
    def testFocalblur(self):
        cmd.viewport(100, 100)
        cmd.fragment('gly', 'm1')
        cmd.focalblur(4.0, 3)
