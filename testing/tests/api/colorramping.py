
from pymol import cmd, testing, stored, colorramping

@testing.requires_version('1.7.1')
class TestColorramping(testing.PyMOLTestCase):

    def _sample_data(self):
        cmd.fragment('his')
        cmd.alter('all', 'b = 20')
        cmd.map_new('map')

    def testVolume(self):
        '''
        cmd.volume
        cmd.volume_color
        cmd.volume_ramp_new
        '''
        self._sample_data()

        cmd.volume('vol', 'map', '2fofc')
        ramp = cmd.volume_color('vol')
        self.assertTrue(len(ramp) > 5)
        self.assertAlmostEqual(ramp[5], 1.0, delta=0.01)

        cmd.volume_ramp_new('named', [
            2.0, 'blue', .1,
            3.0, 'yellow', .2])
        cmd.volume_color('vol', 'named')
        ramp = cmd.volume_color('vol')
        self.assertEqual(len(ramp), 10)
        self.assertArrayEqual(ramp, [2., 0., 0., 1., .1, 3., 1., 1., 0., .2], delta=0.01)

    @testing.requires('gui')
    def testVolumePanel(self):
        self.skipTest("no tk")
        self._sample_data()
        cmd.volume('vol', 'map')
        cmd.volume_panel('vol')
