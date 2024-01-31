'''
Crash when making object discrete
'''

from pymol import cmd, testing, stored

@testing.requires_version('1.7.1')
class Test1917(testing.PyMOLTestCase):

    def test(self):
        cmd.load(self.datafile('small01.mol'), discrete=0)
        cmd.load(self.datafile('small01.mol'), discrete=0)
        cmd.load(self.datafile('small01.mol'), discrete=0)
        cmd.load(self.datafile('small01.mol'), discrete=0)
        cmd.load(self.datafile('small01.mol'), discrete=1)
