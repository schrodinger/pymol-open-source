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

        # involve pymol._colortype (holds a bound method)
        self._test_colortype(p1.cmd)

        p1.stop()

        weak_p = weakref.ref(p1)
        weak_c = weakref.ref(p1.cmd)

        del p1

        self.assertEqual(None, weak_p())
        self.assertEqual(None, weak_c())

    def _test_colortype(self, _self):
        colorset = set()
        _self.set('sphere_color', 'blue', '(*)')
        _self.iterate('*', 'colorset.add(('
                'tuple(s.sphere_color),'
                'int(s.sphere_color)))',
                space={'colorset': colorset, 'tuple': tuple, 'int': int})
        self.assertEqual(list(colorset), [((0., 0., 1.), 2)])
