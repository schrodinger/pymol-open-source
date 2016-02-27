'''
selection macros
'''

from pymol import cmd, CmdException, testing, stored

class TestPYMOL2720(testing.PyMOLTestCase):

    def testSeleMacro(self):
        cmd.fab('GGGGG',  'foo-bar-com')
        cmd.fragment('gly', 'foo-com-bar')

        self.assertTrue( cmd.count_atoms('/foo-bar-com'))
        self.assertTrue( cmd.count_atoms('/foo-bar-com////c*'))
        self.assertTrue( cmd.count_atoms('/foo-bar/'))
        self.assertTrue( cmd.count_atoms('/foo-bar*/'))
        self.assertTrue( cmd.count_atoms('/foo*/'))
        self.assertEqual(cmd.count_atoms('/foo-bar-com///2:4/ca'), 3)
        self.assertEqual(cmd.count_atoms('/foo-bar-com///2-4/ca'), 3)
        self.assertTrue( cmd.count_atoms('foo-bar-com////'))

        # test ambiguous prefix
        # -1 on failure
        cmd.set('raise_exceptions', 0)
        self.assertEqual(cmd.count_atoms('/foo/'), -1)
