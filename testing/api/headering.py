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
