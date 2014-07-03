import os
from pymol import cmd, testing, stored

def touch(filename):
    open(filename, 'ab').close()

class TestExterning(testing.PyMOLTestCase):

    def testCdLsPwd(self):
        with testing.mkdtemp() as path:
            cmd.cd(path)
            self.assertEqual(os.getcwd(),
                    os.path.realpath(path))

            touch('foo1.txt')
            touch('foo2.txt')
            touch('foo3.bin')

            cmd.feedback("disable", "python", "output")

            cmd.pwd()
            # no test of output possible

            cmd.ls('*.txt')
            # no test of output possible

    def testPaste(self):
        cmd.paste
        # partly functional command (?)

    def testSystem(self):
        cmd.system
        self.skipTest("TODO")

