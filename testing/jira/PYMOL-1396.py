'''
PYMOL-1356
rotate ANISOU when aligning molecules
'''

from pymol import cmd, testing, stored

class Test1356(testing.PyMOLTestCase):

    @testing.foreach(True, False)
    def test(self, use_shader):

        pdbfile = self.datafile("sampletrajectory.pdb")
        dcdfile = self.datafile("sampletrajectory.dcd")
        cmd.load(pdbfile)
        cmd.load_traj(dcdfile)
        cmd.set('use_shaders', use_shader)
        cmd.spheroid("SampleTrajectory", 100)
        cmd.set('ambient', 1)
        cmd.set('specular', 0)
        cmd.show_as('spheres')
        cmd.color('blue')
        self.assertImageHasColor('blue')
