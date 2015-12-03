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
        # check that all elements of v are the same, need to strip number off beginning,
        # remove "." and strip trailing 0's
        s = set([self._get_num(str(s)).replace('.', '').rstrip('0') for s in v])
        self.assertEqual(len(s), 1)
