from pymol import cmd, testing

# fixed in released version 1.8.4.2
@testing.requires_version('1.8.4.1')
@testing.requires('incentive')
class TestLoadMaeTableSubBlocks(testing.PyMOLTestCase):

    def test(self):
        cmd.load(self.datafile('subblock.mae'))
        self.assertEqual(5, cmd.count_atoms())
