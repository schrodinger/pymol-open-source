import os
from pymol import cmd, testing, stored, headering

mtzfiles = [
    '/schrodinger/mmshare/mmlibs/mmreflect/10_warpNtrace.mtz',
    '/schrodinger/tasks/other_maps/ec_1.1.mtz',
    '/schrodinger/tasks/other_maps/refine.mtz',
    '/schrodinger/tasks/other_maps/refine_without_eIF1A.mtz',
    '/schrodinger/tmp/1atx4.mtz',
]

class TestHeadering(testing.PyMOLTestCase):
    def test(self):
        for filename in mtzfiles:
            if not os.path.exists(filename):
                continue
            header = headering.MTZHeader(filename)
            self.assertTrue(isinstance(header, headering.MTZHeader))
            self.assertTrue(header.datasets)

    def test2(self):
        filename = self.datafile('4rwb.mtz')
        header = headering.MTZHeader(filename)

        self.assertTrue(isinstance(header, headering.MTZHeader))

        ds0 = header.datasets['0']
        self.assertEqual(ds0['name'], 'HKL_base')
        self.assertEqual(ds0['project'], 'HKL_base')
        self.assertAlmostEqual(float(ds0['alpha']), 90.0)
        self.assertAlmostEqual(float(ds0['beta']), 95.677)
        self.assertAlmostEqual(float(ds0['gamma']), 90.0)

        ds1 = header.datasets['1']
        self.assertEqual(ds1['name'], 'data_1')
        self.assertEqual(ds1['project'], 'sf_convert')
        self.assertAlmostEqual(float(ds1['x']), 40.93)
        self.assertAlmostEqual(float(ds1['y']), 41.23)
        self.assertAlmostEqual(float(ds1['z']), 27.93)

        self.assertEqual(header.getColumns(), [
            'HKL_base/HKL_base/H', 'HKL_base/HKL_base/K',
            'HKL_base/HKL_base/L', 'cryst_1/data_1/FREE', 'cryst_1/data_1/FP',
            'cryst_1/data_1/SIGFP', 'cryst_1/data_1/FC', 'cryst_1/data_1/PHIC',
            'cryst_1/data_1/FOM'
        ])

        self.assertEqual(header.getColumnsOfType("W"), ['cryst_1/data_1/FOM'])
