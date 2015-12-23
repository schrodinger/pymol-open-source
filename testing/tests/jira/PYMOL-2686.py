'''
PYMOL-2686 new feature: distance mode=4 (centroids)
'''

from pymol import cmd, CmdException, testing, stored

class TestPYMOL2686(testing.PyMOLTestCase):

    @testing.requires_version('1.8.1.0')
    def testDistanceMode4(self):
        # create two 2-state objects
        cmd.fragment('gly', 'm1', origin=0)
        cmd.create('m1', 'm1', 1, 2)
        cmd.create('m2', 'm1', 1, 1)
        cmd.create('m2', 'm1', 1, 2)

        # shift one state by 5 and the other by 9 angstrom
        cs = cmd.get_coordset('m1', state=1, copy=0)
        cs += [5., 0., 0.]
        cs = cmd.get_coordset('m2', state=2, copy=0)
        cs += [9., 0., 0.]

        # creates one distance measure per state
        d = cmd.distance('d1', '?m1', '?m2', mode=4)

        # expecting 7 = (5 + 9) / 2
        self.assertAlmostEqual(d, 7.0, delta=1e-4)

        # one distance measurement has 2 * 3 coordinates
        d1states = cmd.get_session('d1', 1, 1, 0, 0)['names'][0][5][2]
        self.assertEqual(len(d1states[0][1]), 6)
        self.assertEqual(len(d1states[1][1]), 6)
