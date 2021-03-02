import pymol
from pymol import cmd, testing

@testing.requires_version('2.3')
class TestCTest(testing.PyMOLTestCase):

    def test2(self):
        try:
            status = pymol._cmd.test2()
        except NotImplementedError:
            self.skipTest("not compiled with --testing")
        self.assertTrue(status in (0, None))
