import pymol
from pymol import cmd, testing, stored

foos = ['foo']
ba_s = ['bar', 'baz']
coms = ['com', 'com_bla', 'com_xxx']
words = foos + ba_s + coms

class TestShortcut(testing.PyMOLTestCase):

    def testShortcut(self):
        # build shortcut
        sc = cmd.Shortcut(words)
        
        # get all keywords
        self.assertItemsEqual(words, sc.interpret(''))

        # full/prefix hits
        self.assertEqual('foo', sc.interpret('f'))
        self.assertEqual('foo', sc.interpret('fo'))
        self.assertEqual('foo', sc.interpret('foo'))

        self.assertItemsEqual(ba_s,  sc.interpret('b'))
        self.assertItemsEqual(ba_s,  sc.interpret('ba'))
        self.assertEqual('bar', sc.interpret('bar'))

        self.assertItemsEqual(coms, sc.interpret('c'))
        self.assertItemsEqual(coms, sc.interpret('co'))
        self.assertEqual('com', sc.interpret('com'))

        # add one
        sc.append('foo_new')
        self.assertItemsEqual(['foo', 'foo_new'], sc.interpret('f'))
        self.assertEqual('foo', sc.interpret('foo'))
        self.assertEqual('foo_new', sc.interpret('foo_'))
        
        self.assertEqual(False, sc.has_key(''))


        # abbreviations
        self.assertEqual('foo_new', sc.interpret('f_'))
        self.assertEqual('foo_new', sc.interpret('f_new'))
        self.assertEqual('foo_new', sc.interpret('fo_'))
        self.assertEqual('com_xxx', sc.interpret('c_x'))
        self.assertEqual('com_xxx', sc.interpret('c_xxx'))
        self.assertEqual('com_xxx', sc.interpret('co_x'))

        # missing key
        self.assertEqual(None, sc.interpret('missing_key'))

        # auto error
        self.assertEqual(None, sc.auto_err(''))
        self.assertEqual(None, sc.auto_err('missing_key'))
        self.assertItemsEqual(coms, sc.auto_err('co'))
        self.assertEqual('com', sc.auto_err('com'))

    def testShortcutMode1(self):
        # build shortcut
        sc = cmd.Shortcut(words)

        # full/prefix hits
        self.assertEqual('foo', sc.interpret('f', 1))
        self.assertItemsEqual(coms, sc.interpret('com', 1))

        # add one
        sc.append('foo_new')
        self.assertItemsEqual(['foo', 'foo_new'], sc.interpret('foo', 1))

    def testShortcutRebuild(self):
        sc = cmd.Shortcut(words)
        sc.rebuild(coms)

        self.assertEqual(None, sc.interpret('f'))
        self.assertEqual(None, sc.interpret('foo'))

        self.assertItemsEqual(coms, sc.interpret('c'))
        self.assertItemsEqual(coms, sc.interpret('com', 1))
        self.assertEqual('com', sc.interpret('com'))
        self.assertEqual('com_xxx', sc.interpret('c_x'))

