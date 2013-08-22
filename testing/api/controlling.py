
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
        cmd.mask
        self.skipTest("TODO")

    def testMouse(self):
        cmd.mouse
        self.skipTest("TODO")

    def testOrder(self):
        cmd.order
        self.skipTest("TODO")

    def testSetKey(self):
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
        cmd.unmask
        self.skipTest("TODO")

