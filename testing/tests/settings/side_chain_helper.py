'''
test side chain helper settings

cartoon_side_chain_helper
cartoon_nucleic_acid_mode
ribbon_side_chain_helper
ribbon_nucleic_acid_mode
'''

from pymol import cmd, testing

class TestSideChainHelper(testing.PyMOLTestCase):

    @testing.foreach.product(
            ['cartoon', 'ribbon'],  # backbone rep
            ['sticks', 'lines'],    # sidechain rep
            ['cartoon', 'ribbon'],  # for setting X_side_chain_helper and X_nucleic_acid_mode
            [0, 1])                 # X_nucleic_acid_mode
    @testing.requires_version('1.7.7')
    def test(self, bb_rep, sc_rep, bb_set, n_a_m):
        cmd.viewport(100, 100)

        # color classes
        color = {'cartoon': 'blue', 'ribbon': 'red'}
        sc_color = 'white'
        op_color = 'green'
        p_color = 'yellow'

        # lighting setup for color testing
        self.ambientOnly()
        cmd.set('dynamic_width', 0)
        cmd.set('line_width', 4)
        cmd.set('ribbon_width', 4)

        # settings setup
        cmd.set(bb_set + '_side_chain_helper')
        cmd.set(bb_set + '_nucleic_acid_mode', n_a_m)
        cmd.set('ribbon_color', color['ribbon'])
        cmd.set('cartoon_color', color['cartoon'])

        # data
        cmd.load(self.datafile('1ehz-5.pdb'))
        cmd.orient()

        # atom colors
        cmd.color(sc_color)
        cmd.color(op_color, "name OP1+OP2")         # not visible with SCH
        cmd.color(p_color, "name P+O3'+C5'+O5'")    # not visible with SCH and NAM=1

        # need to check for OP1 and O1P naming scheme, so alter some atoms
        cmd.alter('name OP1 & resi 1-2', 'name="O1P"')
        cmd.alter('name OP2 & resi 1-2', 'name="O2P"')

        # reps
        cmd.show_as(bb_rep)
        cmd.show(sc_rep)

        # test
        img = self.get_imagearray()
        self.assertImageHasColor(sc_color, img) # always visible
        for bb_test in ['cartoon', 'ribbon']:
            self._assertImageHasColor(bb_rep == bb_test, color[bb_test], img, 0, bb_test + ' wrong')
        self._assertImageHasColor(not (bb_rep == bb_set), op_color, img, 0, 'OP wrong')
        self._assertImageHasColor(not (bb_rep == bb_set and n_a_m == 1), p_color, img, 0, 'NAM=1 wrong')
