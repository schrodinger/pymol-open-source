from pymol import cmd, testing, stored

@testing.requires_version('1.7.1')
class TestCIF(testing.PyMOLTestCase):

    def testPYMOL1737(self):
        cmd.load(self.datafile('ice_IV.cif'))
        self.assertTrue('ice_IV' in cmd.get_object_list())

    def testPYMOL1533(self):
        cmd.load(self.datafile('1v5a-3models.cif'))
        self.assertEqual(cmd.count_states(), 3)
