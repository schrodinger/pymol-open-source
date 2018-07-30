import os
import sys
from pymol import cmd, testing, stored

@testing.requires('incentive')
class TestPh4(testing.PyMOLTestCase):

    @testing.requires_version('1.7.4')
    def test_load_hypothesis_xyz(self):
        from epymol import ph4
        ph4.load_hypothesis_xyz(
                self.datafile('phase/hypothesis.xyz'),
                ref=self.datafile('phase/ALL.xyz'))
        extent = cmd.get_extent('hypothesis')
        self.assertArrayEqual(extent,
                [[-4.218, -2.742, -3.702], [10.536, 6.152, 3.536]],
                delta=1e-2)

    @testing.requires_version('1.7.4')
    def test_load_hypothesis_xvol(self):
        from epymol import ph4
        ph4.load_hypothesis_xvol(
                self.datafile('phase/1M7Q_hypo.xvol'))
        extent = cmd.get_extent('1M7Q_hypo')
        self.assertArrayEqual(extent,
                [[33.622, 23.055, 23.741], [54.344, 38.530, 37.948]],
                delta=1e-2)

    @testing.requires_version('1.8.6')
    def test_load_hypothesis_phypo(self):
        from epymol import ph4
        ph4.load_phypo(
                self.datafile('phase/ADHHRR_1_trans.phypo'),
                'g1', zoom=-1, mimic=1, atom_props='', _self=cmd)
        self.assertEqual(cmd.count_atoms('g1.ADHHRR_1 and name A2+H6'), 2)
        self.assertEqual(cmd.count_atoms('g1.endo-1'), 35)
        self.assertEqual(cmd.count_atoms('g1.endo-2'), 37)
        self.assertEqual(cmd.count_atoms('g1.endo-36'), 28)
        self.assertEqual(cmd.count_atoms('g1.ADHHRR_1_xvol'), 936)
        self.assertTrue('g1.ADHHRR_1_cgo' in cmd.get_names_of_type('object:cgo'))
