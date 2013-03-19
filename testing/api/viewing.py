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
        cmd.clip
        self.skipTest('TODO')

    def testOrigin(self):
        cmd.origin
        self.skipTest('TODO')

    def testMove(self):
        cmd.move
        self.skipTest('TODO')

    def testEnable(self):
        cmd.enable
        self.skipTest('TODO')

    def testDisable(self):
        cmd.disable
        self.skipTest('TODO')

    def testToggle(self):
        cmd.toggle
        self.skipTest('TODO')

    def testShow(self):
        cmd.show
        self.skipTest('TODO')

    def testHide(self):
        cmd.hide
        self.skipTest('TODO')

    def testView(self):
        cmd.view
        self.skipTest('TODO')

    def testScene(self):
        cmd.scene
        self.skipTest('TODO')

    def testStereo(self):
        cmd.stereo
        self.skipTest('TODO')

    def testTurn(self):
        cmd.turn
        self.skipTest('TODO')

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
        cmd.viewport
        self.skipTest('TODO')

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
        cmd.ray
        self.skipTest('TODO')

    def testRefresh(self):
        cmd.refresh
        self.skipTest('TODO')

    def testReset(self):
        cmd.reset
        self.skipTest('TODO')

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
        cmd.color
        self.skipTest('TODO')

    def testSpectrum(self):
        cmd.spectrum
        self.skipTest('TODO')
