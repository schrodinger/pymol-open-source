'''
PYMOL-2667
Atom-level unique ids PSE loading
'''

from pymol import cmd, CmdException, testing, stored

class TestPYMOL2667(testing.PyMOLTestCase):

    @testing.requires_version('1.7.6.7')
    def test2667(self):
        cmd.fragment('ala', 'm1')
        cmd.alter_state(1, '%m1', 's.label_placement_offset = [1., 0., 0.]')

        with testing.mktemp('.pse') as filename:
            cmd.save(filename)
            cmd.set_name('m1', 'm2')
            cmd.load(filename, partial=1)

        # now we have two copies (m1 m2). If unique ids are converted upon
        # loading, then there will be no settings cross-leaking. Otherwise
        # changing m1 settings will affect m2.

        cmd.alter_state(1, '%m1', 's.label_placement_offset = [0., 1., 0.]')

        m2_setting = []
        cmd.iterate_state(1,
                'first %m2', 'm2_setting.append(s.label_placement_offset)',
                space=locals())

        self.assertArrayEqual([1., 0., 0.], m2_setting[0], 0.001)
