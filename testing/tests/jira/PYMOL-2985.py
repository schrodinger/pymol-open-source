'''
"clean" depends on atom sorting
'''

from pymol import cmd, testing

@testing.requires_version('2.1')
@testing.requires('incentive', 'no_edu')
class TestClean(testing.PyMOLTestCase):

    def test(self):
        cmd.load(self.datafile('PYMOL-2985.sdf'), 'm1')
        cset_before = cmd.get_coordset('m1')
        cmd.clean('m1')
        cset_after = cmd.get_coordset('m1')
        self.assertArrayNotEqual(cset_before, cset_after)
