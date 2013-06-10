import random
from pymol import cmd, testing, stored

class TestViewing(testing.PyMOLTestCase):

    def assertViewIs(self, view):
        view = [round(v, 2) for v in view]
        ref = [round(v, 2) for v in cmd.get_view()]
        self.assertEqual(view, ref, 'view differs')

    def _load_default_scene(self):
        cmd.load('viewing-ref/ref.pse')
        cmd.set_view((
            1.0,    0.0,    0.0,
            0.0,    1.0,    0.0,
            0.0,    0.0,    1.0,
            0.0,    0.0,  -43.0,
            2.6,   -1.9,    0.0,
           36.4,   49.6,  -20.0))

    def _testZoomGeneric(self, func, views):
        self._load_default_scene()
        cmd.frame(2)

        for (state, view) in views.items():
            func(state=state)
            self.assertViewIs(view)

        func(state=-1)
        self.assertViewIs(views[2])

    def testZoom(self):
        views = {
            0: (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -30.417539596557617,
                1.666666865348816, -1.6666665077209473, 0.0, 23.981420516967773,
                36.853660583496094, -20.0),
            1: (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -14.178204536437988,
                0.0, 0.0, 0.0, 11.178204536437988, 17.178203582763672, -20.0),
            2: (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -14.178204536437988,
                4.999999046325684, 0.0, 0.0, 11.178204536437988, 17.178203582763672, -20.0),
        }
        self._testZoomGeneric(cmd.zoom, views)

        cmd.zoom('all', 2.0, 3, 1)
        self.assertViewIs((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
            0.0, 0.0, -25.093609, 0.0, -5.0, 0.0, 19.783993, 30.403225, -20.0))

    def testCenter(self):
        views = {
            0: (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
                0.0, 0.0, -43.0, 1.666666865348816, -1.6666665077209473, 0.0,
                36.400001525878906, 49.599998474121094, -20.0),
            1: (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
                0.0, 0.0, -43.0, 0.0, 0.0, 0.0,
                36.400001525878906, 49.599998474121094, -20.0),
            2: (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
                0.0, 0.0, -43.0, 4.999999046325684, 0.0, 0.0,
                36.400001525878906, 49.599998474121094, -20.0),
        }
        self._testZoomGeneric(cmd.center, views)

    def testOrient(self):
        views = {
            0: (0.6800623536109924, -0.7324044704437256, 0.03314966708421707, 0.7331299781799316,
                0.6789760589599609, -0.03888334706425667, 0.005970506463199854, 0.05074610933661461,
                0.9986937046051025, 0.0, 0.0, -30.417539596557617, 1.666666865348816,
                -1.6666665077209473, 0.0, 23.981420516967773, 36.853660583496094, -20.0),
            1: (-0.4805402159690857, -0.8505961298942566, 0.21346519887447357, 0.7515078186988831,
                -0.5248664617538452, -0.39968883991241455, 0.4520145058631897, -0.031645797193050385,
                0.8914490342140198, 0.0, -0.0, -14.178204536437988, 0.0, 0.0, 0.0, 11.178204536437988,
                17.178203582763672, -20.0),
            2: (-0.4805402159690857, -0.8505960702896118, 0.21346530318260193, 0.7515078186988831,
                -0.52486652135849, -0.3996887803077698, 0.4520145058631897, -0.03164568170905113,
                0.8914490938186646, 0.0, -0.0, -14.178204536437988, 4.999999046325684, 0.0, 0.0,
                11.178204536437988, 17.178203582763672, -20.0),
        }
        self._testZoomGeneric(cmd.orient, views)

    def testClip(self):
        for a in "xyz":
            cmd.turn(a, random.random() * 10 - 5)
            cmd.move(a, random.random() * 10 - 5)
        v = cmd.get_view()
        cmd.clip("near", -5)
        self.assertAlmostEqual(v[15] + 5, cmd.get_view()[15], delta=1e-3)
        cmd.clip("far", 10)
        self.assertAlmostEqual(v[16] - 10, cmd.get_view()[16], delta=1e-3)
        cmd.clip("move", -15)
        a = cmd.get_view()
        self.assertAlmostEqual(v[15], a[15] - 20, delta=1e-3)
        self.assertAlmostEqual(v[16], a[16] - 5, delta=1e-3)
        cmd.clip("slab", 20)
        v = cmd.get_view()
        self.assertAlmostEqual(v[16] - v[15], 20.0, delta=1e-3)
        cmd.pseudoatom()
        cmd.clip("atoms", 5, "all")
        v = cmd.get_view()
        self.assertAlmostEqual(v[16] - v[15], 10.0, delta=1e-3)

    def testOrigin(self):
        cmd.origin
        self.skipTest('TODO')

    def testMove(self):
        for a in "xyz":
            cmd.turn(a, random.random() * 10 - 5)
            cmd.move(a, random.random() * 10 - 5)
        v = list(cmd.get_view())
        d = (2,4,6)
        cmd.move("x", d[0])
        cmd.move("y", d[1])
        cmd.move("z", d[2])
        m = cmd.get_view()
        v[9] += d[0]
        v[10] += d[1]
        v[11] += d[2]
        v[15] -= d[2]
        v[16] -= d[2]
        self.assertArrayEqual(v, m, delta=1e-3)

    def testEnable(self):
        cmd.create('m1', 'none')
        cmd.create('m2', 'none')
        cmd.disable()
        self.assertEqual(cmd.get_names('public_objects', 1), [])
        cmd.enable('m1')
        self.assertEqual(cmd.get_names('public_objects', 1), ['m1'])
        cmd.enable()
        self.assertEqual(cmd.get_names('public_objects', 1), ['m1', 'm2'])

    def testDisable(self):
        # see testEnable
        pass

    def testToggle(self):
        cmd.toggle
        self.skipTest('TODO')

    def testShow(self):
        cmd.fragment('ala')
        self.assertEqual(cmd.count_atoms('rep sticks'), 0)
        cmd.show('sticks')
        self.assertEqual(cmd.count_atoms('rep sticks'), 10)
        cmd.hide('lines', 'not elem C')
        self.assertEqual(cmd.count_atoms('rep lines'), 3)

    def testHide(self):
        # see testShow
        pass

    def testView(self):
        cmd.turn('x', 30)
        a = cmd.get_view()
        cmd.view('A', 'store')
        cmd.turn('y', 30)
        self.assertNotEqual(a, cmd.get_view())
        cmd.view('A', 'recall', animate=0)
        self.assertEqual(a, cmd.get_view())

    def testScene(self):
        cmd.scene
        self.skipTest('TODO')

    def testStereo(self):
        cmd.stereo
        self.skipTest('TODO')

    def testTurn(self):
        cmd.turn('z', 90)
        v = cmd.get_view()
        self.assertArrayEqual(v[:9], [0.,1.,0.,-1.,0.,0.,0.,0.,1.], delta=1e-3)

    def testRock(self):
        cmd.rock
        self.skipTest('TODO')

    def testLabel(self):
        cmd.label
        self.skipTest('TODO')

    def testLabel2(self):
        cmd.label2
        self.skipTest('TODO')

    def testWindow(self):
        cmd.window
        self.skipTest('TODO')

    def testViewport(self):
        v = (100,50)
        cmd.viewport(*v)
        self.assertEqual(cmd.get_viewport(), v)

    def testCartoon(self):
        cmd.cartoon
        self.skipTest('TODO')

    def testCapture(self):
        cmd.capture
        self.skipTest('TODO')

    def testDraw(self):
        cmd.draw
        self.skipTest('TODO')

    def testRay(self):
        # tested in many other tests
        pass

    def testRefresh(self):
        cmd.refresh
        self.skipTest('TODO')

    def testReset(self):
        # view
        v = cmd.get_view()
        cmd.turn('x', 10)
        cmd.move('y', 10)
        self.assertNotEqual(v, cmd.get_view())
        cmd.reset()
        self.assertEqual(v, cmd.get_view())
        # object
        cmd.pseudoatom("m1")
        x = cmd.get_object_matrix("m1")
        cmd.translate([1,2,3], object="m1")
        self.assertNotEqual(x, cmd.get_object_matrix("m1"))
        cmd.reset("m1")
        self.assertEqual(x, cmd.get_object_matrix("m1"))

    def testDirty(self):
        cmd.dirty
        self.skipTest('TODO')

    def testRebuild(self):
        cmd.rebuild
        self.skipTest('TODO')

    def testRecolor(self):
        cmd.recolor
        self.skipTest('TODO')

    def testColor(self):
        cmd.fragment('ala')
        cmd.color(3)
        colors = set()
        cmd.iterate('all', 'colors.add(color)', space=locals())
        self.assertItemsEqual(colors, [3])

    def _testSpectrum_setup(self):
        cmd.pseudoatom('m1', pos=(-2,0,0))
        cmd.pseudoatom('m2', pos=(0,0,0))
        cmd.pseudoatom('m3', pos=(2,0,0))
        cmd.show_as('spheres')
        cmd.set('ambient', 1)
        cmd.set('specular', 0)
        cmd.set('reflect', 0)
        cmd.set('direct', 0)
        cmd.viewport(40,20)
        cmd.zoom()

    def testSpectrum(self):
        self._testSpectrum_setup()
        cmd.spectrum()
        img = self.get_imagearray()
        self.assertImageHasColor('blue', img, delta=1)
        self.assertImageHasColor('red', img)

    def testSpectrumany(self):
        self._testSpectrum_setup()
        cmd.spectrum('count', 'red blue')
        img = self.get_imagearray()
        self.assertImageHasColor('red', img)
        self.assertImageHasColor('blue', img)
        self.assertImageHasColor('0x7f007f', img, delta=30)
