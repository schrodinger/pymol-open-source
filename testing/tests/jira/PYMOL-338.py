'''
auto_show_classified
'''

from pymol import cmd, CmdException, testing, stored

n_total      = 1435
n_polymer    = 1268
n_organic    = 48
n_inorganic  = 1
n_solvent    = 118
n_org_inorg  = n_organic + n_inorganic

@testing.requires_version('1.8.1.0')
class TestPYMOL338(testing.PyMOLTestCase):

    @testing.foreach(
            (0, 1, 0, 0, 0, 0, 0, n_total, n_total),
            (1, 1, n_polymer, 0, n_organic, n_organic, n_inorganic, n_solvent, n_solvent),
            (2, 1, n_polymer, 0, n_organic, n_organic, n_inorganic, n_total, n_total),
            (3, 0, 0, n_polymer, 0, 0, 0, n_org_inorg, n_org_inorg),
        )
    def testAutoShowClassified(self, asc, asl,
            n_cartoon,
            n_ribbon,
            n_sticks,
            n_nb_spheres,
            n_spheres,
            n_lines,
            n_nonbonded):
        cmd.set('auto_show_classified', asc)
        cmd.set('auto_show_lines', asl)
        cmd.set('auto_show_nonbonded', asl)

        cmd.load(self.datafile('1rx1.pdb'))

        self.assertEqual(cmd.count_atoms('rep lines'), n_lines)
        self.assertEqual(cmd.count_atoms('rep nonbonded'), n_nonbonded)

        self.assertEqual(cmd.count_atoms('rep cartoon'), n_cartoon)
        self.assertEqual(cmd.count_atoms('rep sticks'), n_sticks)
        self.assertEqual(cmd.count_atoms('rep nb_spheres'), n_nb_spheres)
        self.assertEqual(cmd.count_atoms('rep spheres'), n_spheres)
