'''
mmpymolx specific functionality

See also: tests/jira/PYMOL-317.py (text_type test)
'''

from pymol import cmd, testing, stored
import unittest

@testing.requires('incentive')
@testing.requires('no_edu')
class TestMMPyMOLx(testing.PyMOLTestCase):

    def testStereo(self):
        cmd.fragment('ala')
        cmd.remove('hydro')

        # default: S configuration
        labels = []
        cmd.iterate('name CA', 'labels.append(stereo)', space=locals())
        self.assertEqual(labels, ['S'])

        # load R configuration (moves CB)
        ala_conf_R = {
            'N'  : (-0.67690, -1.23030, -0.49050),
            'CA' : (-0.00090,  0.06370, -0.49050),
            'C'  : ( 1.49910, -0.11030, -0.49050),
            'O'  : ( 2.03010, -1.22730, -0.50150),
        #   'CB' : (-0.50890,  0.85570,  0.72650), # S configuration
            'CB' : (-0.33784,  0.82664, -1.78310), # R configuration
        }
        cmd.alter_state(1, 'ala', '(x,y,z) = ala_conf_R.get(name)', space=locals())
        labels = []
        cmd.iterate('name CA', 'labels.append(stereo)', space=locals())
        self.assertEqual(labels, ['R'])

    def testStereoPYMOL2782(self):
        # L-Cysteine has R configuration
        cmd.fragment('cys')
        labels = []
        cmd.iterate('name CA', 'labels.append(stereo)', space=locals())
        self.assertEqual(labels, ['R'])

    def testLoadMTZ(self):
        cmd.load_mtz(self.datafile('4rwb.mtz'), 'foo',
                'cryst_1/data_1/FP',
                'cryst_1/data_1/PHIC')

        self.assertEqual(cmd.get_names(), ['foo'])

        extent = cmd.get_extent('foo')
        self.assertArrayEqual(extent, [
            [   0.000,   0.000,   0.000],
            [  37.585,  40.586,  27.214]], delta=1e-3, msg='extend wrong')

        # min, max, mean, stdev
        mmms = cmd.get_volume_histogram('foo', 0)
        self.assertArrayEqual(mmms, [-2.6767, 4.9998, 0.0, 1.0], delta=1e-4,
                msg='histogram wrong')
