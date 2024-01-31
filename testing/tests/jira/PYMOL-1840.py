from pymol import cmd, testing, stored

@testing.requires_version('1.8.3.1')
class TestPyMOL1840(testing.PyMOLTestCase):

    def testSdfV3000(self):
        cmd.load(self.datafile('1rx1.pdb'))
        with testing.mktemp('.sdf') as filename:
            cmd.save(filename)
            cmd.delete('*')
            cmd.load(filename)
        self.assertEqual(cmd.count_atoms(), 1435)

        # check for capitalized element symbols
        cmd.set('ignore_case', 0)
        self.assertEqual(cmd.count_atoms('elem Ca'), 1)
