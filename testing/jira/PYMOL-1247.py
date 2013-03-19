'''
PYMOL-1247 commands with LITERAL parsing mode and leading whitespace fail
'''

from pymol import cmd, testing, stored, parsing

def myfunc(a='adefault', b='bdefault'):
    stored.v = (a, b)

class TestPYMOL1247(testing.PyMOLTestCase):

    @testing.foreach('STRICT', 'LITERAL2')
    def testLeadingSpace(self, mode):
        cmd.keyword['myfunc'] = [myfunc, 0, 0, '', getattr(parsing, mode)]
        cmd.do(' myfunc x', 0, 0, 1)
        self.assertEqual(stored.v, ('x', 'bdefault'))
