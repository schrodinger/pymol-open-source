'''
Testing template for PyMOL testing
'''

from pymol import cmd, testing

class TestVersion(testing.PyMOLTestCase):
    def _get_num(self, x):
        retlist = []
        for ele in x:
            if ele.isdigit() or ele == '.':
                retlist.append(ele)
            else:
                break
        return ''.join(retlist)

    def test_version(self):
        v = cmd.get_version()[:3]
        # check types
        self.assertTrue(isinstance(v[0], str))
        self.assertTrue(isinstance(v[1], float))
        self.assertTrue(isinstance(v[2], int))
        # check that major of str and float match
        self.assertEqual(int(v[0].split('.')[0]), int(v[1]))
