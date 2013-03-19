
from pymol import cmd, testing, stored

class TestCommanding(testing.PyMOLTestCase):

    def testAlias(self):
        stored.v = None
        cmd.alias('foo', '/stored.v = 123')
        cmd.do('foo', echo=0)
        self.assertEqual(stored.v, 123)

    def testCls(self):
        cmd.cls
        self.skipTest("TODO")

    def testDelete(self):
        cmd.pseudoatom('m1')
        cmd.pseudoatom('m2')
        cmd.pseudoatom('m3')
        cmd.delete('m1 m2')
        self.assertEqual(cmd.get_names(), ['m3'])

    def testDo(self):
        # tested with other methods
        pass

    def testDummy(self):
        v = cmd.dummy(1, 2, x=3, _self=cmd)
        self.assertEqual(v, None)

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
        cmd.log
        self.skipTest("TODO")

    def testLogClose(self):
        cmd.log_close
        self.skipTest("TODO")

    def testLogOpen(self):
        cmd.log_open
        self.skipTest("TODO")

    def testQuit(self):
        cmd.quit
        self.skipTest("TODO")

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
        cmd.resume
        self.skipTest("TODO")

    def testSplash(self):
        cmd.splash
        self.skipTest("TODO")

    def testSync(self):
        cmd.sync
        self.skipTest("TODO")

