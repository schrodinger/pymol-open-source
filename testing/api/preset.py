
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
