'''
order with groups
'''

from pymol import cmd, testing

@testing.requires('incentive')
class Test1382(testing.PyMOLTestCase):

    def test(self):
        names_ungrouped = ['a0', 'a1', 'b0', 'b1', 'c0', 'c1', 'g1']
        names_grouped   = ['a0', 'a1', 'c0', 'c1', 'b0', 'g1', 'b1']

        for name in names_ungrouped[:-1]:
            cmd.pseudoatom(name)
        cmd.group(names_ungrouped[-1])

        self.assertEqual(cmd.get_names(), names_ungrouped)

        cmd.group('g1', 'b*')
        cmd.order('b0 g1', location="upper")

        self.assertEqual(cmd.get_names(), names_grouped)
