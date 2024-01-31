from pymol import cmd, testing

# fixed in 1.8.2.4 and 1.8.4.2
@testing.requires_version('1.8.2.3')
@testing.requires('incentive')
class TestLoadMaeTableSubBlocks(testing.PyMOLTestCase):

    def test(self):
        cmd.load(self.datafile('subblock.mae'))
        self.assertEqual(5, cmd.count_atoms())
