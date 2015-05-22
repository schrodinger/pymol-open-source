
from pymol import cmd, testing, stored

class TestSetting(testing.PyMOLTestCase):

    @testing.requires_version('1.7.5')
    def test_indices(self):
        # setting indices must not change, since they are used in session files
        self.assertEqual(742, cmd.pymol.setting._get_index('collada_geometry_mode'))

    def testGet(self):
        name_bool = 'cartoon_fancy_helices'
        name_obj = 'ala'
        cmd.fragment('ala', name_obj)
        cmd.set(name_bool, 1, name_obj)
        self.assertEqual('off', cmd.get(name_bool))
        self.assertEqual('on', cmd.get(name_bool, name_obj))
        cmd.unset(name_bool, '*')
        cmd.set(name_bool)
        self.assertEqual('on', cmd.get(name_bool))
        self.assertEqual('on', cmd.get(name_bool, name_obj))
        cmd.set(name_bool, 0, name_obj)
        self.assertEqual('on', cmd.get(name_bool))
        self.assertEqual('off', cmd.get(name_bool, name_obj))

    def testGetBond(self):
        # see testSetBond
        pass

    def testGetSettingBoolean(self):
        cmd.get_setting_boolean
        self.skipTest("TODO")

    def testGetSettingFloat(self):
        cmd.get_setting_float
        self.skipTest("TODO")

    def testGetSettingInt(self):
        cmd.get_setting_int
        self.skipTest("TODO")

    def testGetSettingText(self):
        cmd.get_setting_text
        self.skipTest("TODO")

    def testGetSettingTuple(self):
        cmd.get_setting_tuple
        self.skipTest("TODO")

    def testGetSettingUpdates(self):
        cmd.get_setting_updates
        self.skipTest("TODO")

    @testing.foreach.product(
            ('', 'ala', '(elem O)'),
            ('sphere_scale',),
            (2.3,),
            (1.0,),
            )
    @testing.requires('incentive')
    def testSet(self, sele, name, value, defaultvalue):
        cmd.fragment('ala')
        cmd.set(name, value, sele)
        n = cmd.iterate('first (%s)' % (sele or 'all'), 'stored.v = s.' + name)
        self.assertEqual(n, 1)
        self.assertAlmostEqual(stored.v, value)

        if sele:
            cmd.unset(name, sele)
            n = cmd.iterate('first (%s)' % (sele or 'all'), 'stored.v = s.' + name)
            self.assertEqual(n, 1)
            self.assertEqual(stored.v, defaultvalue)

    def testSetBond(self):
        value = 2.3
        cmd.fragment('ala')
        cmd.set_bond('stick_radius', value, '*', '*')
        v_list = cmd.get_bond('stick_radius', 'first *', '*')
        self.assertAlmostEqual(v_list[0][1][0][2], value)

        # unset
        # TODO get_bond reports [0,0,0] after unset_bond
        # cmd.unset_bond('stick_radius', '*', '*')
        # v_list = cmd.get_bond('stick_radius', 'first *', '*')
        # self.assertEqual(v_list, [])

    def testUnset(self):
        # see testSet
        pass

    def testUnsetBond(self):
        # see testSetBond
        pass
