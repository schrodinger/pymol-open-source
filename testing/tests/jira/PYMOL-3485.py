'''
PYMOL-3485 Partial load of unique color settings
PYMOL-3486 Partial load of ramps
'''

from pymol import cmd, testing

@testing.requires_version('2.5')
class Test3485(testing.PyMOLTestCase):

    def test(self):
        session1 = cmd.get_session()

        cmd.fragment('gly', 'm1')
        cmd.ramp_new('ramp1', 'none')
        cmd.set('surface_color', 'ramp1', '(m1)')
        session2 = cmd.get_session()

        cColorExtCutoff = cmd.get_color_index('ramp1')

        cmd.set_session(session1)
        cmd.ramp_new('ramp2', 'none')

        self.assertEqual(cmd.get_color_index('ramp1'), -1)
        self.assertEqual(cmd.get_color_index('ramp2'), cColorExtCutoff)

        cmd.set_session(session2, partial=1)

        self.assertEqual(cmd.get_color_index('ramp1'), cColorExtCutoff - 1)
        self.assertEqual(cmd.get_color_index('ramp2'), cColorExtCutoff)

        colors = set()
        cmd.iterate('m1', 'colors.add(s.surface_color)', space=locals())
        self.assertEqual(list(colors), [cColorExtCutoff - 1])
