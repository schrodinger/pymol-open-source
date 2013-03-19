'''
Testing template for PyMOL testing
'''

from pymol import cmd, testing

class TestVersion(testing.PyMOLTestCase):

    def test_version(self):
        v = cmd.get_version()
        # check that all elements of v are the same, need to remove "." and strip trailing 0's
        s = set(map(lambda s: str(s).replace('.', '').rstrip('0'), v))
        self.assertEqual(len(s), 1)
