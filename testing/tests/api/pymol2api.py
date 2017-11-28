from __future__ import absolute_import

from pymol import testing, invocation

class TestPyMOL2(testing.PyMOLTestCase):

    @testing.requires_version('1.8.5')
    def testMultiInstance(self):
        import pymol
        import pymol2

        p1 = pymol2.PyMOL()
        p2 = pymol2.PyMOL()
        p3 = pymol # singleton

        p1.start()
        p2.start()

        p1.cmd.fragment('ala')
        p2.cmd.fragment('trp')
        p3.cmd.fragment('ile') # singleton

        self.assertEqual(p1.cmd.count_atoms(), 10)
        self.assertEqual(p2.cmd.count_atoms(), 24)
        self.assertEqual(p3.cmd.count_atoms(), 19) # singleton

        p1.stop()
        p2.stop()

    @testing.requires_version(
            '2.1' if invocation.options.incentive_product else
            '1.9')
    def testDel(self):
        import pymol2
        import weakref

        p1 = pymol2.PyMOL()
        p1.start()
        p1.cmd.fragment('ala')
        p1.stop()

        weak_p = weakref.ref(p1)
        weak_c = weakref.ref(p1.cmd)

        del p1

        self.assertEqual(None, weak_p())
        self.assertEqual(None, weak_c())
