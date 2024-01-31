
from pymol import cmd, testing, stored

class TestControlling(testing.PyMOLTestCase):

    def testButton(self):
        cmd.button
        self.skipTest("TODO")

    def testConfigMouse(self):
        cmd.config_mouse
        self.skipTest("TODO")

    def testEditMode(self):
        cmd.edit_mode
        self.skipTest("TODO")

    def testMask(self):
        cmd.fragment('ala')
        cmd.mask('elem C')
        self.assertEqual(3, cmd.count_atoms('masked'))
        cmd.unmask('name CA')
        self.assertEqual(2, cmd.count_atoms('masked'))

    def testMouse(self):
        cmd.mouse
        self.skipTest("TODO")

    def testOrder(self):
        # proper group ordering only in incentive, see jira/PYMOL-1382.py

        for i in range(3):
            cmd.pseudoatom('m%i' % i)
        self.assertEqual(cmd.get_names(), ['m0', 'm1', 'm2'])

        cmd.order('m2 m1')
        self.assertEqual(cmd.get_names(), ['m0', 'm2', 'm1'])

        cmd.order('m2 m1', location="top")
        self.assertEqual(cmd.get_names(), ['m2', 'm1', 'm0'])

        cmd.order('m2', location="bottom")
        self.assertEqual(cmd.get_names(), ['m1', 'm0', 'm2'])

        cmd.order('m2 m0', location="upper")
        self.assertEqual(cmd.get_names(), ['m1', 'm2', 'm0'])

        cmd.order('*', 'yes')
        self.assertEqual(cmd.get_names(), ['m0', 'm1', 'm2'])

    def testSetKey(self):
        if testing.PYMOL_VERSION[1] > 1.84:
            cmd.set('auto_show_classified', 0)

        cmd.fragment('gly')
        N = cmd.count_atoms()
        self.assertEqual(0, cmd.count_atoms('rep sticks | rep spheres'))
        cmd.set_key('F3', cmd.show, ('sticks',))
        cmd._special(3, 0, 0)
        self.assertEqual(N, cmd.count_atoms('rep sticks'))
        cmd.do('set_key F3, hide sticks; show spheres')
        cmd._special(3, 0, 0)
        self.assertEqual(0, cmd.count_atoms('rep sticks'))
        self.assertEqual(N, cmd.count_atoms('rep spheres'))

    def testUnmask(self):
        # see testMask
        pass

