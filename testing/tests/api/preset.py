
from pymol import cmd, testing, stored, preset

class TestPreset(testing.PyMOLTestCase):

    def _assertSimpleReps(self):
        self.assertEqual(cmd.count_atoms('rep cartoon'), 0)
        self.assertEqual(cmd.count_atoms('rep ribbon & guide'), cmd.count_atoms('guide'))
        self.assertEqual(cmd.count_atoms('rep sticks & organic'), cmd.count_atoms('organic'))
        self.assertEqual(cmd.count_atoms('rep nonbonded & solvent'), cmd.count_atoms('solvent'))

    def _assertPrettyReps(self):
        self.assertEqual(cmd.count_atoms('rep cartoon & guide'), cmd.count_atoms('guide'))
        self.assertEqual(cmd.count_atoms('rep ribbon'), 0)
        self.assertEqual(cmd.count_atoms('rep sticks & organic'), cmd.count_atoms('organic'))
        self.assertEqual(cmd.count_atoms('rep nonbonded'), 0)
        self.assertEqual(cmd.count_atoms('rep nonbonded'), 0)

    def testSimple(self):
        cmd.load(self.datafile('1oky.pdb.gz'))
        preset.simple('*')
        self._assertSimpleReps()

    def testTechnical(self):
        cmd.load(self.datafile('1oky.pdb.gz'))
        preset.technical('*')
        self._assertSimpleReps()
        self.assertTrue('1oky_pol_conts' in cmd.get_names('all'))

    def testPublication(self):
        cmd.load(self.datafile('1oky.pdb.gz'))
        preset.publication('*')
        self._assertPrettyReps()
        self.assertTrue(cmd.get_setting_boolean('cartoon_smooth_loops', '1oky'))
        self.assertTrue(cmd.get_setting_boolean('cartoon_fancy_helices', '1oky'))

    def testPretty(self):
        cmd.load(self.datafile('1oky.pdb.gz'))
        preset.pretty('*')
        self._assertPrettyReps()
        self.assertFalse(cmd.get_setting_boolean('cartoon_smooth_loops', '1oky'))
        self.assertFalse(cmd.get_setting_boolean('cartoon_fancy_helices', '1oky'))

    def testLigands(self):
        cmd.load(self.datafile('1oky.pdb.gz'))
        preset.ligands('*')
        self.assertEqual(cmd.count_atoms('rep cartoon'), 0)
        self.assertGreater(cmd.count_atoms('rep ribbon'), 0)
        self.assertEqual(cmd.count_atoms('rep sticks & organic'), cmd.count_atoms('organic'))
        solvent_nb_count = cmd.count_atoms('rep nonbonded & solvent')
        self.assertGreater(solvent_nb_count, 0)
        self.assertLess(solvent_nb_count, cmd.count_atoms('solvent'))
        self.assertGreater(cmd.count_atoms('rep lines & polymer'), 0)
        self.assertEqual(cmd.count_atoms('(polymer be. 10 of hetatm) and rep lines'), 0)
        self.assertTrue('1oky_pol_conts' in cmd.get_names('all'))

    def testDefault(self):
        cmd.load(self.datafile('1oky.pdb.gz'))
        cmd.color("yellow", "1oky") # color index 6
        preset.default('1oky')
        self.assertEqual(cmd.count_atoms('rep cartoon'), 0)
        self.assertEqual(cmd.count_atoms('rep ribbon'), 0)
        self.assertEqual(cmd.count_atoms('rep sticks'), 0)
        self.assertEqual(cmd.count_atoms('rep spheres'), 0)
        self.assertEqual(cmd.count_atoms('rep lines'), cmd.count_atoms())
        self.assertEqual(cmd.count_atoms('rep nonbonded & solvent'), cmd.count_atoms('solvent'))
        # object color must be preserved for carbons
        color_set = set()
        cmd.iterate('elem C', 'color_set.add(color)', space=locals())
        self.assertItemsEqual(color_set, [6])
