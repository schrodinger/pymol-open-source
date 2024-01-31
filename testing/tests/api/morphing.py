from pymol import cmd, testing, stored

@testing.requires('incentive')
@testing.requires('no_edu')
class TestMorphing(testing.PyMOLTestCase):

    def testMorphRigimol(self):
        cmd.set('suspend_undo')
        import epymol.rigimol
        cmd.fab('ACD', 'm1')
        cmd.create('m1', 'm1', 1, 2)
        cmd.rotate('x', 90, 'm1', 2)
        steps = 5
        cmd.morph('mout', 'm1', refinement=1, steps=steps, method='rigimol')
        self.assertEqual(steps, cmd.count_states('mout'))
        self.assertEqual(cmd.count_atoms('m1'),
                         cmd.count_atoms('mout'))

    def testMorphMulti(self):
        cmd.set('suspend_undo')
        cmd.fab('ACD', 'm1')
        cmd.create('m1', 'm1', 1, 2)
        cmd.create('m1', 'm1', 1, 3)
        cmd.rotate('x', 90, 'm1', 2)
        cmd.rotate('x', 45, 'm1', 3)
        steps = 5
        cmd.morph('mout', 'm1', None, 0, 0, 0, steps, 'linear')
        self.assertEqual(steps * 2, cmd.count_states('mout'))
        self.assertEqual(cmd.count_atoms('m1'),
                         cmd.count_atoms('mout'))
        cmd.morph('mo2', 'm1', None, 0, -1, 0, steps, 'linear')
        self.assertEqual(steps * 3, cmd.count_states('mo2'))
        self.assertEqual(cmd.count_atoms('m1'),
                         cmd.count_atoms('mo2'))
