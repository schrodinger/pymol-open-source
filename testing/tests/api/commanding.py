
import sys
import pymol
import __main__
from pymol import cmd, testing, stored

class TestCommanding(testing.PyMOLTestCase):

    def testAlias(self):
        stored.v = None
        cmd.alias('foo', '/stored.v = 123')
        cmd.do('_ foo', echo=0)
        self.assertEqual(stored.v, 123)

    @testing.requires('gui', 'no_run_all')
    def testCls(self):
        cmd.set('internal_prompt', 0)
        cmd.set('text', 1)
        cmd.cls()

        # check on black screen
        img = self.get_imagearray()
        self.assertFalse(img[...,:3].any())

    def testDelete(self):
        cmd.pseudoatom('m1')
        cmd.pseudoatom('m2')
        cmd.pseudoatom('m3')
        cmd.delete('m1 m2')
        self.assertEqual(cmd.get_names(), ['m3'])

    def testDo(self):
        # tested with other methods
        pass

    def testExtend(self):
        def check(v):
            stored.v = None
            cmd.do('foo', echo=0)
            self.assertEqual(stored.v, v)

        cmd.extend('foo', lambda: setattr(stored, 'v', 123))
        check(123)

        @cmd.extend
        def foo():
            stored.v = 456
        check(456)

    def testLog(self):
        with testing.mktemp('.pml') as logfile:
            cmd.log_open(logfile)
            cmd.do('_ color blue')
            cmd.log('hello world')
            cmd.log_close()
            lines = filter(None, map(str.strip, open(logfile)))
            self.assertEqual(lines, ['color blue', 'hello world'])

    def testLogClose(self):
        # see testLog
        pass

    def testLogOpen(self):
        # see testLog
        pass

    def testQuit(self):
        cmd.quit
        self.skipTest("cannot test quit")

    def testReinitialize(self):
        def check(v, names):
            self.assertEqual(v, cmd.get('pdb_conect_all'))
            self.assertEqual(cmd.get_names(), names)

        cmd.pseudoatom('m1')
        v = cmd.set('pdb_conect_all')
        cmd.reinitialize('store_defaults')

        cmd.reinitialize('original_settings')
        check('off', ['m1'])

        cmd.reinitialize('settings')
        check('on', ['m1'])

        cmd.reinitialize('purge_defaults')
        cmd.reinitialize('settings')
        check('off', ['m1'])

        cmd.reinitialize('everything')
        check('off', [])

    def testResume(self):
        with testing.mktemp('.pml') as logfile:
            with open(logfile, 'w') as f:
                print >> f, 'bg yellow'
            cmd.resume(logfile)
            self.assertEqual('yellow', cmd.get('bg_rgb'))
            cmd.log('hello world')
            cmd.log_close()
            lines = filter(None, map(str.strip, open(logfile)))
            self.assertEqual(lines, ['bg yellow', 'hello world'])

    def testSplash(self):
        cmd.feedback('disable', 'all', 'output')
        cmd.splash()

    def testSync(self):
        cmd.sync()

    @testing.foreach(
        ('local', 'pymol', False),
        ('global', 'pymol', True),
        ('module', '', False),
        ('main', '__main__', True),
        ('private', '__main__', False),
    )
    def testRun(self, namespace, mod, rw):
        stored.tmp = False
        with testing.mktemp('.py') as filename:
            varname = '_tmp_' + namespace
            with open(filename, 'w') as handle:
                print >> handle, 'from pymol import stored'
                if mod:
                    print >> handle, 'stored.tmp = (__name__ == "%s")' % (mod)
                else:
                    print >> handle, 'stored.tmp = True'
                print >> handle, varname + ' = True'
            cmd.do('run %s, %s' % (filename, namespace), 0, 0)
            self.assertTrue(stored.tmp)
            if mod:
                self.assertEqual(rw, hasattr(sys.modules[mod], varname))
