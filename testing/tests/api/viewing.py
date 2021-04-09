import sys
import random
from pymol import cmd, testing, stored

expr_list = [
    'segi', 'chain', 'resn', 'resi', 'name', 'alt', 'elem', 'text_type',
    'formal_charge', 'numeric_type', 'ID',
    'q', 'b', 'partial_charge', 'vdw',
]

class TestViewing(testing.PyMOLTestCase):

    def _assertIterate(self, expr, selection):
        stored.v = []
        cmd.iterate(selection, 'stored.v.append(%s)' % expr)
        self.assertEqual(len(stored.v), cmd.count_atoms(selection))
        return stored.v

    def assertEqualIterate(self, a, b, selection='*', msg=None):
        for va, vb in self._assertIterate('(%s,%s)' % (a, b), selection):
            for t in (float, int):
                if isinstance(va, t) or isinstance(vb, t):
                    va, vb = t(va), t(vb)
            self.assertEqual(va, vb, '%s=%s != %s=%s' % (a, repr(va), b, repr(vb)))

    def assertTrueIterate(self, expr, selection='*', msg=None):
        for v in self._assertIterate(expr, selection):
            self.assertTrue(v, expr)

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
        from chempy import cpv
        cmd.pseudoatom('m1')
        cmd.pseudoatom('m2')
        cmd.pseudoatom('m3', pos=[1,0,0])
        # by selection
        cmd.origin('m3')
        cmd.rotate('y', 90, 'm1')
        # by position
        cmd.origin(position=[-1,0,0])
        cmd.rotate('y', 90, 'm2')
        coords = []
        cmd.iterate_state(1, 'm1 m2', 'coords.append([x,y,z])', space=locals())
        d = cpv.distance(*coords)
        self.assertAlmostEqual(d, 2 * 2**0.5)

    @testing.requires("multi_undo")
    def testUndoOrigin(self):
        from chempy import cpv
        cmd.pseudoatom('m1')
        cmd.pseudoatom('m2')
        cmd.pseudoatom('m3', pos=[1,0,0])
        # by selection
        cmd.origin('m3')
        cmd.rotate('y', 90, 'm1')
        # by position
        cmd.origin(position=[-1,0,0])
        cmd.undo()
        cmd.rotate('y', 90, 'm2')
        coords = []
        cmd.iterate_state(1, 'm1 m2', 'coords.append([x,y,z])', space=locals())
        d = cpv.distance(*coords)
        self.assertAlmostEqual(d, 0)

        cmd.delete("*")
        cmd.pseudoatom('m1')
        cmd.pseudoatom('m2')
        cmd.pseudoatom('m3', pos=[1,0,0])
        # by selection
        cmd.origin('m3')
        cmd.rotate('y', 90, 'm1')
        # by position
        cmd.origin(position=[-1,0,0])
        cmd.undo()
        cmd.redo()
        cmd.rotate('y', 90, 'm2')
        coords = []
        cmd.iterate_state(1, 'm1 m2', 'coords.append([x,y,z])', space=locals())
        d = cpv.distance(*coords)
        self.assertAlmostEqual(d, 2 * 2**0.5)

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
        if testing.PYMOL_VERSION[1] > 1.84:
            cmd.set('auto_show_classified', 0)

        cmd.fragment('gly')
        cmd.show('sticks')
        cmd.show('spheres')
        cmd.toggle()
        self.assertEqual(cmd.count_atoms('rep lines'), 0)
        self.assertEqual(cmd.count_atoms('rep sticks'), 7)
        self.assertEqual(cmd.count_atoms('rep spheres'), 7)
        cmd.toggle()
        cmd.toggle('sticks', '(all)')
        self.assertEqual(cmd.count_atoms('rep lines'), 7)
        self.assertEqual(cmd.count_atoms('rep sticks'), 0)
        self.assertEqual(cmd.count_atoms('rep spheres'), 7)
        cmd.toggle('sticks', '(all)')
        cmd.toggle('spheres', '(all)')
        self.assertEqual(cmd.count_atoms('rep lines'), 7)
        self.assertEqual(cmd.count_atoms('rep sticks'), 7)
        self.assertEqual(cmd.count_atoms('rep spheres'), 0)

    def testShow(self):
        if testing.PYMOL_VERSION[1] > 1.84:
            cmd.set('auto_show_classified', 0)

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
        cmd.fragment('ala', 'm1')
        cmd.fragment('gly', 'm2')
        cmd.create('m2', 'm2', 1, 2)
        cmd.create('m2', 'm2', 1, 3)
        cmd.create('m2', 'm2', 1, 4)

        # store
        cmd.show_as('sticks', 'm1')
        cmd.show_as('spheres', 'm2')
        cmd.color('blue', 'm1')
        cmd.color('yellow', 'm2')
        view_001 = cmd.get_view()
        cmd.scene('new', 'store', 'hello world')

        # store
        cmd.frame(3)
        cmd.show_as('lines')
        cmd.color('red')
        cmd.turn('x', 45)
        view_002 = cmd.get_view()
        cmd.scene('new', 'store')

        # we actually don't know the auto naming counter
        # self.assertEqual(cmd.get_scene_list(), ['001', '002'])
        names = cmd.get_scene_list()

        # recall
        cmd.scene(names[0], 'recall', animate=0)
        self.assertArrayEqual(view_001, cmd.get_view(), delta=1e-3)
        self.assertEqual(cmd.count_atoms('m1'), cmd.count_atoms('color blue'))
        self.assertEqual(cmd.count_atoms('m1'), cmd.count_atoms('rep sticks'))
        self.assertEqual(cmd.count_atoms('m2'), cmd.count_atoms('color yellow'))
        self.assertEqual(cmd.count_atoms('m2'), cmd.count_atoms('rep spheres'))
        self.assertNotEqual(cmd.get_wizard(), None)
        self.assertEqual(cmd.get_state(), 1)

        # recall
        cmd.scene(names[1], 'recall', animate=0)
        self.assertArrayEqual(view_002, cmd.get_view(), delta=1e-3)
        self.assertEqual(0, cmd.count_atoms('color blue'))
        self.assertEqual(0, cmd.count_atoms('rep sticks'))
        self.assertEqual(cmd.count_atoms(), cmd.count_atoms('color red'))
        self.assertEqual(cmd.count_atoms(), cmd.count_atoms('rep lines'))
        self.assertEqual(cmd.get_wizard(), None)
        self.assertEqual(cmd.get_state(), 3)

        # with movie (not playing) it must not recall the state
        cmd.mset('1-4')
        cmd.frame(1)
        cmd.scene(names[1], 'recall', animate=0)
        self.assertEqual(cmd.get_state(), 1)

        # rename and delete
        cmd.scene('F2', 'store')
        cmd.scene(names[0], 'rename', new_key='F1')
        cmd.scene(names[1], 'delete')
        self.assertEqual(cmd.get_scene_list(), ['F1', 'F2'])

    def testSceneInsertBeforeAfter(self):
        # insert_before and insert_after
        cmd.scene('F1', 'store')
        cmd.scene('F2', 'store')
        cmd.scene('F1')
        cmd.scene('new', 'insert_after')
        cmd.scene('new', 'insert_before')
        names = cmd.get_scene_list()
        self.assertEqual(names[0], 'F1')
        self.assertEqual(names[3], 'F2')
        # we don't know the auto naming counter, but we know this:
        self.assertTrue(names[1] > names[2])

    def testSceneOrder(self):
        cmd.scene('Z0', 'store')
        cmd.scene('F1', 'store')
        cmd.scene('F10', 'store')
        cmd.scene('F2', 'store')
        cmd.scene_order('*', True)
        self.assertEqual(cmd.get_scene_list(), ['F1', 'F2', 'F10', 'Z0'])
        cmd.scene_order('F10 F1')
        self.assertEqual(cmd.get_scene_list(), ['F2', 'F10', 'F1', 'Z0'])
        cmd.scene_order('F10 F1', location='top')
        self.assertEqual(cmd.get_scene_list(), ['F10', 'F1', 'F2', 'Z0'])
        cmd.scene_order('F10 F1', location='bottom')
        self.assertEqual(cmd.get_scene_list(), ['F2', 'Z0', 'F10', 'F1'])

    @testing.requires_version('2.5')
    def testSceneOrderWithSpaces(self):
        cmd.scene('Z0 1 2', 'store')
        cmd.scene('F1 ABC', 'store')
        cmd.scene('F10 XX', 'store')
        cmd.scene_order('*', True)
        self.assertEqual(cmd.get_scene_list(), ['F1 ABC', 'F10 XX', 'Z0 1 2'])
        cmd.scene_order(['F10 XX', 'F1 ABC'])
        self.assertEqual(cmd.get_scene_list(), ['F10 XX', 'F1 ABC', 'Z0 1 2'])

    @testing.requires('gui')
    def testStereo(self):
        for (k,v) in cmd.stereo_dict.items():
            if v < 1:
                continue
            cmd.stereo(k)
            self.assertTrue(cmd.get_setting_int('stereo'))
            self.assertEqual(cmd.get_setting_int('stereo_mode'), v)
        cmd.stereo('off')
        self.assertFalse(cmd.get_setting_int('stereo'))
        cmd.stereo('on')
        self.assertTrue(cmd.get_setting_int('stereo'))
        shift = cmd.get_setting_float('stereo_shift')
        cmd.stereo('swap')
        self.assertAlmostEqual(cmd.get_setting_float('stereo_shift'), -shift)

    def testTurn(self):
        cmd.turn('z', 90)
        v = cmd.get_view()
        self.assertArrayEqual(v[:9], [0.,1.,0.,-1.,0.,0.,0.,0.,1.], delta=1e-3)

    def testRock(self):
        cmd.rock(1)
        self.assertTrue(cmd.get_setting_int('rock'))
        cmd.rock(0)
        self.assertFalse(cmd.get_setting_int('rock'))
        cmd.rock(-1)
        self.assertTrue(cmd.get_setting_int('rock'))
        cmd.rock(-1)
        self.assertFalse(cmd.get_setting_int('rock'))

    @testing.foreach(cmd.label, cmd.label2)
    def testLabel(self, func):
        cmd.pseudoatom('m1')
        for expr in expr_list + ['resn + " " + resi']:
            func('*', expr)
            self.assertEqualIterate('label', expr)
        func()
        self.assertTrueIterate('label == ""')

    @testing.requires_version('1.8.7')
    def testLabelUnicode(self):
        if sys.version_info[0] < 3:
            uchr = unichr
        else:
            uchr = chr

        stored.L_unicode = u''.join(uchr(i) for i in range(32, 0x1FF))
        stored.L_utf8 = stored.L_unicode.encode('utf-8')

        self.assertTrue(not isinstance(stored.L_unicode, bytes))
        self.assertTrue(isinstance(stored.L_utf8, bytes))

        cmd.pseudoatom('pseudo_unicode', label=stored.L_unicode)
        cmd.pseudoatom('pseudo_utf8', label=stored.L_utf8)
        cmd.pseudoatom('label_unicode')
        cmd.pseudoatom('label_utf8')
        cmd.label('label_unicode', 'stored.L_unicode')
        cmd.label('label_utf8', 'stored.L_utf8')

        stored.check_unicode = []
        cmd.iterate('all', 'stored.check_unicode.append(label == stored.L_unicode)')
        self.assertTrue(all(stored.check_unicode))

    def testWindow(self):
        cmd.window
        self.skipTest('TODO')

    def testViewport(self):
        v = (100,50)
        cmd.viewport(*v)
        self.assertEqual(cmd.get_viewport(), v)

    def testCartoon(self):
        cmd.load(self.datafile("1t46-frag.pdb"))
        cmd.cartoon("skip")
        cmd.cartoon("rect", "resi 590-593")
        cmd.cartoon("oval", "resi 595-598")
        cmd.cartoon("tube", "resi 600-602")
        cmd.cartoon("putty", "resi 606-620")
        cmd.cartoon("dumbbell", "resi 623-625")
        cmd.set("cartoon_highlight_color", "red")
        cmd.color("yellow")
        cmd.show_as("cartoon")
        cmd.viewport(200, 100)
        cmd.set_view(
            (-0.102381177, 0.9807688, -0.16615206, -0.781365454, 0.024080699,
             0.623609841, 0.615618646, 0.193670169, 0.763873398, -0.000125334,
             -0.000184208, -65.135032654, 24.068731308, 26.194671631,
             53.642494202, 40.266975403, 89.999656677, -20.))
        self.ambientOnly()
        self.assertImageEqual("viewing-ref/cartoon.png", delta=1)

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
        self.ambientOnly()
        cmd.viewport(100, 100)
        cmd.pseudoatom("m1", pos=(0, 0, 0), vdw=1)
        cmd.pseudoatom("m2", pos=(0, 0, -4), vdw=1)
        cmd.zoom()
        cmd.show('spheres')
        cmd.color("blue", "m1")
        cmd.color("red", "m2")
        self.assertImageHasNotColor('red')
        cmd.alter("m2", "vdw=2")
        self.assertImageHasNotColor('red')
        cmd.rebuild()
        self.assertImageHasColor('red')

    @testing.foreach("spheres", "nb_spheres", "sticks", "lines", "nonbonded",
                     "surface")
    def testRecolor(self, rep):
        self.ambientOnly()
        cmd.viewport(100, 100)
        cmd.fragment("gly")
        if rep in ('nb_spheres', 'nonbonded'):
            cmd.unbond("*", "*")
        if rep in ('lines', 'nonbonded'):
            cmd.set("line_width", 5)
        cmd.zoom()
        cmd.color("blue")
        cmd.show_as(rep)
        self.assertImageHasColor('blue')
        cmd.alter("all", "color = 0x40FF0000")
        self.assertImageHasColor('blue')
        cmd.recolor()
        self.assertImageHasColor('red')

    def testColor(self):
        cmd.fragment('ala')
        cmd.color(3)
        colors = set()
        cmd.iterate('all', 'colors.add(color)', space=locals())
        self.assertItemsEqual(colors, [3])

    @testing.requires_version('2.4')
    def test_color_deep(self):
        cmd.viewport(100, 70)
        self.ambientOnly()

        cmd.fragment('trp', 'm1')
        cmd.orient('m1')
        cmd.show_as('sticks')
        cmd.show('spheres', 'name C')
        cmd.set('stick_color', 'blue', 'm1')
        cmd.set('stick_color', 'red', 'name CH2+CZ3+CE3+CD2')
        cmd.set('sphere_color', 'yellow', 'name C')

        cmd.color('green')

        img = self.get_imagearray()
        self.assertImageHasColor('blue', img)
        self.assertImageHasColor('red', img)
        self.assertImageHasColor('yellow', img)
        self.assertImageHasNotColor('green', img)

        cmd.color_deep('green')

        img = self.get_imagearray()
        self.assertImageHasNotColor('blue', img)
        self.assertImageHasNotColor('red', img)
        self.assertImageHasNotColor('yellow', img)
        self.assertImageHasColor('green', img)

    def _testSpectrum_setup(self):
        cmd.pseudoatom('m1', pos=(-2,0,0))
        cmd.pseudoatom('m2', pos=(0,0,0))
        cmd.pseudoatom('m3', pos=(2,0,0))
        cmd.show_as('spheres')
        self.ambientOnly()
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

    def testSetColor(self):
        self._testSpectrum_setup()
        cmd.color('red')
        self.assertImageHasColor('red')

        # redefine existing name, [0-1] float range
        cmd.set_color('red', [0., 0., 1.])
        cmd.recolor()
        self.assertImageHasColor('blue')

        # add new name, [0-1] float range
        cmd.set_color('thomas', [0., 1., 0.])
        cmd.color('thomas')
        self.assertImageHasColor('green')

        # redefine existing (custom) name, [0-255] int range
        cmd.set_color('thomas', [0x00, 0x80, 0x80])
        cmd.recolor()
        self.assertImageHasColor('0x008080')

        # add new very long name, [0-255] int range
        longname = 'some' + 'very' * 10 + 'longcolorname'
        cmd.set_color(longname, [0xFF, 0xFF, 0x00])
        cmd.color(longname)
        self.assertImageHasColor('yellow')
