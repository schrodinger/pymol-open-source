
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
        for v_ref in (0, 1):
            cmd.set('orthoscopic', v_ref)
            v = cmd.get_setting_boolean('orthoscopic')
            self.assertTrue(isinstance(v, int))
            self.assertEqual(v, v_ref)

    def testGetSettingFloat(self):
        for v_ref in (0.0, 1.0):
            # bool as float
            cmd.set('orthoscopic', v_ref)
            v = cmd.get_setting_float('orthoscopic')
            self.assertTrue(isinstance(v, float))
            self.assertEqual(v, v_ref)

        # float with 6 significant digits
        cmd.set('sphere_scale', 1.234565)
        v = cmd.get_setting_float('sphere_scale')
        self.assertEqual(int(v * 1e5), 123456)

    @testing.requires_version('2.4')
    def testGetSettingFloatExact(self):
        # PyMOL 2.4 uses pymol::pretty_f2d which rounds 7 significant digits
        for v_ref in (1.234567, 1234567.e6, 9999999.0):
            cmd.set('sphere_scale', v_ref)
            v = cmd.get_setting_float('sphere_scale')
            self.assertEqual(v, v_ref)

    def testGetSettingInt(self):
        cmd.set('light_count', 4)
        v = cmd.get_setting_int('light_count')
        self.assertTrue(isinstance(v, int))
        self.assertEqual(v, 4)

        # float as int gets floored
        cmd.set('sphere_scale', 3.7)
        v = cmd.get_setting_int('sphere_scale')
        self.assertTrue(isinstance(v, int))
        self.assertEqual(v, 3)

    def testGetSettingText(self):
        # int
        cmd.set('light_count', 4)
        v = cmd.get_setting_text('light_count')
        self.assertTrue(isinstance(v, str))
        self.assertEqual(v, '4')

        # bool
        cmd.set('orthoscopic', 0)
        v = cmd.get_setting_text('orthoscopic')
        self.assertEqual(v, 'off')

        # float gets rounded to 5 digits
        cmd.set('sphere_scale', 3.7)
        v = cmd.get_setting_text('sphere_scale')
        self.assertEqual(v, '3.70000')

        # float3 gets nicely formatted
        cmd.set('label_position', (0.100009, 2.3, 4.56789))
        v = cmd.get_setting_text('label_position')
        self.assertEqual(v, '[ 0.10001, 2.30000, 4.56789 ]')

        # string
        cmd.set('pdb_echo_tags', 'Hello World')
        v = cmd.get_setting_text('pdb_echo_tags')
        self.assertEqual(v, 'Hello World')

    def testGetSettingTuple(self):
        cmd.set('label_position', (1, 2, 4))
        v = cmd.get_setting_tuple('label_position')
        self.assertEqual(v, (4, (1.0, 2.0, 4.0)))

    def testGetSettingUpdates(self):
        v = cmd.get_setting_updates() # consume whatever
        cmd.set('orthoscopic', 1)
        v = cmd.get_setting_updates()
        self.assertEqual(v, [23])

    @testing.foreach.product(
            ('', 'ala', '(elem O)'),
            ('sphere_scale',),
            (2.3,),
            (1.0,),
            )
    @testing.requires_version('1.7')
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

    @testing.requires_version('1.8.4')
    def testSetBond(self):
        value = 2.3
        cmd.fragment('ala')
        cmd.set_bond('stick_radius', value, '*', '*')
        v_list = cmd.get_bond('stick_radius', 'first *', '*')
        self.assertAlmostEqual(v_list[0][1][0][2], value)

        # unset
        # before 1.8.4, get_bond reported [0,0,0] after unset_bond
        cmd.unset_bond('stick_radius', '*', '*')
        v_list = cmd.get_bond('stick_radius', 'first *', '*')
        self.assertEqual(v_list[0][1][0][2], None)

    @testing.requires_version('2.5')
    def testUnset(self):
        # unset with selection: see testSet

        # unset global settings since PyMOL 2.5:
        for name, value in [
            ('light_count', 5),  # int, non-zero-default
            ('sphere_scale', 2.0),  # float, non-zero default
            ('ray_shadow', 'off'),  # bool, non-zero-default
            ('orthoscopic', 'on'),  # bool, zero-default
            ('cartoon_highlight_color', 'blue'),  # color
            ('fetch_path', '/some/path'),  # string
        ]:
            old_value = cmd.get(name)
            cmd.set(name, value)
            assert old_value != cmd.get(name)
            cmd.unset(name)
            self.assertEqual(old_value, cmd.get(name))

    def testUnsetBond(self):
        # see testSetBond
        pass

    @testing.requires_version('1.8.3')
    def testUnsetDeep(self):
        cmd.fragment('ala')
        cmd.fragment('gly')
        cmd.fragment('his')

        # global
        sphere_scale_global = 0.8
        cmd.set('sphere_scale', sphere_scale_global)
        cmd.set('stick_color', 'blue')

        # object-level
        cmd.set('sphere_scale', .5, 'ala')
        cmd.set('stick_radius', .4, 'gly')

        # atom-level
        cmd.set('sphere_scale', .3, 'index 1-3')
        cmd.set('sphere_color', 'yellow', 'index 2-4')

        # bond-level
        cmd.set('stick_radius', .6, 'elem C')
        cmd.set('stick_color', 'red', 'index 1-5')

        cmd.unset_deep()
        names = cmd.get_object_list()

        for oname in names:
            # object-level check
            for sname in ['sphere_scale', 'sphere_color']:
                self.assertEqual(cmd.get(sname), cmd.get(sname, oname))

            # atom-level check (1)
            # type float
            sname = 'sphere_scale'
            a_level_values = set()
            cmd.iterate(oname, 'a_level_values.add(s.' + sname + ')', space=locals())
            self.assertEqual(len(a_level_values), 1)
            self.assertAlmostEqual(sphere_scale_global, list(a_level_values)[0], delta=1e-4)

            # atom-level check (2)
            # type color ('default' -> None)
            sname = 'sphere_color'
            a_level_values = set()
            cmd.iterate(oname, 'a_level_values.add(s.' + sname + ')', space=locals())
            self.assertEqual(len(a_level_values), 1)
            self.assertEqual(None, list(a_level_values)[0])

        # bond-level check (None for unset settings)
        for sname in ['stick_radius', 'stick_color']:
            b_level_values = cmd.get_bond(sname, '*')
            for o_set in b_level_values:
                for b_set in o_set[1]:
                    self.assertEqual(None, b_set[2], msg=sname + ' ' + o_set[0] + ' ' + str(b_set))
