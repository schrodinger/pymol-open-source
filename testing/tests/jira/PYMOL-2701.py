'''
mol2 formal charge assigment (import)
'''

from pymol import cmd, CmdException, testing, stored

class TestPYMOL2701(testing.PyMOLTestCase):

    @testing.requires_version('1.8.0.4')
    def testMol2FormalChargeAssignment(self):
        cmd.fragment('lys')
        cmd.fragment('lys')

        with testing.mktemp('.mol2') as filename:
            cmd.save(filename, state=0)
            cmd.delete('*')
            cmd.load(filename)

        self.assertEqual(cmd.get_model('name NZ', state=1).atom[0].formal_charge, 1)
        self.assertEqual(cmd.get_model('name NZ', state=2).atom[0].formal_charge, 1)
