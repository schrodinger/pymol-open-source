'''
ccp4 map type 0
'''

from pymol import cmd, CmdException, testing, stored

@testing.requires_version('1.7.1')
class TestPYMOL1887(testing.PyMOLTestCase):

    @testing.foreach(
        (1, [-10.729, 9.836, 0.0, 1.0]),
        (0, [-128.0, 127.0, 5.039, 12.4]),
    )
    def test(self, normalize, ref):
        cmd.set('volume_data_range', 0)
        cmd.set('normalize_ccp4_maps', normalize)
        cmd.load(self.datafile('emd_1155.ccp4'), 'map')
        mmms = cmd.get_volume_histogram('map', 0)
        self.assertArrayEqual(mmms, ref, delta=0.05)
