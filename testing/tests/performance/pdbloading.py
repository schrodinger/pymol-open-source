

'''
Stress testing for pdb loading
'''

import time
from pymol import cmd, testing, stored

@testing.requires('no_run_all')
class StressPDBLoading(testing.PyMOLTestCase):
    @testing.foreach('1rx1.pdb', '1aon.pdb')
    def testPDBLoad(self, dfile):
        tm = []
        for i in range(1, 100):
            start = time.time()
            cmd.load(self.datafile(dfile))
            tm.append(time.time() - start)
        print "min time('", dfile, "')=", min(tm)
