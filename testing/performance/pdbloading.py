

'''
Stress testing for pdb loading
'''

from pymol import cmd, testing, stored

class StressPDBLoading(testing.PyMOLTestCase):
    @testing.foreach('1rx1.pdb', '1aon.pdb')
    def testPDBLoad(self, dfile):
        for i in range(1, 100):
            with self.timing():
                cmd.load(self.datafile(dfile))
        tm = map(lambda x: x[1], self.timing().test.timings)
        self.timing().test.timings = []
        print "min time('", dfile, "')=", min(tm)
