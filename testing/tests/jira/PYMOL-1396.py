'''
PYMOL-1356
spheroid command does not work
'''

from pymol import cmd, testing, stored

@testing.requires('gui')
class Test1356(testing.PyMOLTestCase):

    @testing.foreach(True, False)
    def testSpheresForSpheroidTrajectory(self, use_shader):
        pdbfile = self.datafile("sampletrajectory.pdb")
        dcdfile = self.datafile("sampletrajectory.dcd")
        cmd.load(pdbfile)
        cmd.load_traj(dcdfile)
        cmd.set('use_shaders', use_shader)
        cmd.spheroid("sampletrajectory", 3)
        self.ambientOnly()
        cmd.color('blue')
        cmd.show_as('spheres')
        self.assertImageHasColor('blue')
