import os
from pymol import cmd, testing, stored

@testing.requires('incentive')
@testing.requires('no_edu')
class TestPYMOL1697(testing.PyMOLTestCase):

    @testing.foreach(1, 0)
    @testing.requires_version('1.7.0')
    def testUndoAfterRemoveAtomOnDiscrete(self, discr):
        cmd.set('suspend_undo', 0)
        cmd.load(self.datafile('ligs3d.sdf'), discrete=discr)
        natms = cmd.count_atoms()
        cmd.remove("index 1")
        cmd.undo()
        self.assertEqual(natms, cmd.count_atoms())

    @testing.requires_version('1.7.3')
    def testRegression2207(self):
        '''
        The initial fix for PYMOL-1697 crashed this example
        (PYMOL-2207)
        '''
        cmd.load(self.datafile('1oky-frag.pdb'))
        cmd.create('1oky-frag', 'resi 80-100', 1, 2)
        cmd.remove('resi 85-95')
        # keep ray efford low
        cmd.viewport(50, 50)
        # caused segfault
        img = self.get_imagearray(ray=1)
