'''
ignore_case_chain
'''

from pymol import cmd, CmdException, testing, stored

class TestPYMOL2706(testing.PyMOLTestCase):

    @testing.requires_version('1.8.1.0')
    def testIgnoreCase(self):
        # check defaults
        self.assertEqual(cmd.get_setting_int('ignore_case'), 1)
        self.assertEqual(cmd.get_setting_int('ignore_case_chain'), 0)

        # data
        natoms = 10
        cmd.fragment('ala', 'm1')
        cmd.copy('m2', 'm1')
        cmd.alter('m1', 'chain, segi = "C", "S"')
        cmd.alter('m2', 'chain, segi = "c", "s"')
        cmd.alter('m2', 'resn, name = resn.lower(), name.lower()')

        self.assertEqual(cmd.count_atoms('chain C'), natoms)
        self.assertEqual(cmd.count_atoms('segi  S'), natoms)
        self.assertEqual(cmd.count_atoms('chain c'), natoms)
        self.assertEqual(cmd.count_atoms('segi  s'), natoms)
        self.assertEqual(cmd.count_atoms('resn ALA'), natoms * 2)
        self.assertEqual(cmd.count_atoms('resn ala'), natoms * 2)
        self.assertEqual(cmd.count_atoms('name CA'), 2)
        self.assertEqual(cmd.count_atoms('name ca'), 2)

        cmd.set('ignore_case_chain')

        self.assertEqual(cmd.count_atoms('chain C'), natoms * 2)
        self.assertEqual(cmd.count_atoms('segi  S'), natoms * 2)
        self.assertEqual(cmd.count_atoms('chain c'), natoms * 2)
        self.assertEqual(cmd.count_atoms('segi  s'), natoms * 2)

        cmd.set('ignore_case', 0)

        self.assertEqual(cmd.count_atoms('resn ALA'), natoms)
        self.assertEqual(cmd.count_atoms('resn ala'), natoms)
        self.assertEqual(cmd.count_atoms('name CA'), 1)
        self.assertEqual(cmd.count_atoms('name ca'), 1)
