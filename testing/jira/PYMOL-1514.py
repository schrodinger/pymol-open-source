'''
negative b-factor in MAE file.
'''

from pymol import cmd, testing, stored

class Test1514(testing.PyMOLTestCase):

    def test(self):
        cmd.load("PYMOL-1514.mae")
        b_list = []
        cmd.iterate("all", "b_list.append(b)", space=locals())
        self.assertArrayEqual(b_list, [-1.0, -1.0, -2.3, 2.86], 0.001)
