'''
PYMOL-1356
spheroid command does not work
'''

from pymol import cmd, testing, stored

class Test1356(testing.PyMOLTestCase):

    @testing.requires('gui')
    @testing.foreach(True, False)
    def testSpheresForSpheroidTrajectory(self, use_shader):
        pdbfile = self.datafile("sampletrajectory.pdb")
        dcdfile = self.datafile("sampletrajectory.dcd")
        cmd.load(pdbfile)
        cmd.load_traj(dcdfile)
        cmd.set('use_shaders', use_shader)
        cmd.spheroid("SampleTrajectory", 100)
        cmd.set('ambient', 1)
        cmd.set('specular', 0)
        cmd.color('blue')
        cmd.show_as('spheres')
        self.assertImageHasColor('blue')
