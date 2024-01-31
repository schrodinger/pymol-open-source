import pymol
from pymol import menu, cmd, testing

class TestMenu(testing.PyMOLTestCase):
    def _make_objects(self, n=2):
        for i, aa in enumerate(['gly', 'ala', 'cys', 'his'][:n], 1):
            cmd.fragment(aa, 'm%d' % i)

    def _make_multistate(self, n=3, name='multi1'):
        cmd.fragment('gly', name)
        for i in range(2, n + 1):
            cmd.create(name, name, 1, n)

    def _make_seles(self, n=1):
        self._make_objects(n)
        for i in range(1, n + 1):
            cmd.select('s%i' % i, 'm%d' % i)

    def _make_scenes(self, n=3):
        for i in range(1, n + 1):
            cmd.scene('%03d' % i, 'store')

    def _make_cgo(self, oname):
        from pymol import cgo
        cmd.load_cgo([cgo.STOP], oname)

    def _make_map(self, oname):
        molname = 'm_for_map'
        cmd.fragment('gly', molname)
        cmd.set('gaussian_b_floor', 20)
        cmd.map_new(oname, 'gaussian', 0.5, molname)
        cmd.delete(molname)

    def _make_mesh(self, oname):
        mapname = oname + '_map'
        self._make_map(mapname)
        cmd.isomesh(oname, mapname)

    def _to_dict(self, m):
        d = {}
        for (t, key, value) in m:
            if t == 1:
                d[key] = value
        return d

    def _get_first_evaluable(self, m):
        if not isinstance(m, list):
            return m
        for (t, key, value) in m:
            if t == 1:
                m_sub = self._get_first_evaluable(value)
                if m_sub:
                    return m_sub

    def _eval_first(self, m):
        from pymol import preset, util
        s = self._get_first_evaluable(m)
        if not s:
            raise UserWarning('empty menu')
        if isinstance(s, str):
            exec(s)
        else:
            s()

    def testAlignToObject(self):
        self._make_objects(3)
        m = menu.align_to_object(cmd, 'm2')
        labels = [x[1] for x in m]
        self.assertTrue('m1' in labels)
        self.assertTrue('m2' not in labels)
        self.assertTrue('m3' in labels)

    def testAlignToSele(self):
        self._make_seles(3)
        m = menu.align_to_sele(cmd, 's2')
        labels = [x[1] for x in m]
        self.assertTrue('s1' in labels)
        self.assertTrue('s2' not in labels)
        self.assertTrue('s3' in labels)

    def testAllAction(self):
        self._make_objects(3)
        sele = 'all'
        m = menu.all_action(cmd, sele)
        self.assertEqual(m[0], [2, 'Action:', ''])

    def testCameraMotion(self):
        nstate = 4
        nscene = 3
        self._make_multistate(nstate)
        self._make_scenes(nscene)
        m = menu.camera_motion(cmd)
        d = self._to_dict(m)
        dscene = self._to_dict(d['store with scene'])
        dstate = self._to_dict(d['store with state'])
        self.assertEqual(len(dscene), nscene)
        self.assertEqual(len(dstate), nstate + 1)
        self.assertTrue('current' in dstate)

    def testCgoHide(self):
        sele = 'cgo1'
        self._make_cgo(sele)
        m = menu.cgo_hide(cmd, sele)

    def testCgoShow(self):
        sele = 'cgo1'
        self._make_cgo(sele)
        m = menu.cgo_show(cmd, sele)

    def testColorramps(self):
        cmd.ramp_new('r1', 'none')
        cmd.ramp_new('r2', 'none')
        m = menu.colorramps(cmd, '')
        d = self._to_dict(m)
        self.assertEqual(sorted(d), ['r1', 'r2'])

    @testing.requires_version('1.9')
    def testCopyTo(self):
        nobj = 3
        self._make_objects(nobj)
        sele = 'm1'
        m = menu.copy_to(cmd, sele)
        d = self._to_dict(m)
        self.assertTrue('new' in d)
        self.assertEqual(len(d), nobj)

    def testGeneralColor(self):
        sele = 'cgo1'
        self._make_cgo(sele)
        m = menu.general_color(cmd, sele)
        self._eval_first(m)

    def testGroupAction(self):
        self._make_objects()
        sele = 'g1'
        cmd.group(sele, 'm*')
        m = menu.group_action(cmd, sele)
        d = self._to_dict(m)
        self.assertTrue('rename group' in d)

    @testing.requires_version('1.9')
    def testMainMenu(self):
        m = menu.main_menu(cmd, (0., 0., 0.), (123, 456))

    def testMapAction(self):
        sele = 'map1'
        self._make_map(sele)
        m = menu.map_action(cmd, sele)
        self._eval_first(m)
        m = menu.map_hide(cmd, sele)
        self._eval_first(m)
        m = menu.map_show(cmd, sele)
        self._eval_first(m)

    def testMeasurementColor(self):
        self._make_objects(1)
        sele = 'd1'
        cmd.distance(sele, 'first all', 'last all')
        m = menu.measurement_color(cmd, sele)
        self._eval_first(m)
        m = menu.measurement_hide(cmd, sele)
        self._eval_first(m)
        m = menu.measurement_show(cmd, sele)
        self._eval_first(m)

    def testMeshAction(self):
        sele = 'mesh1'
        self._make_mesh(sele)
        m = menu.mesh_action(cmd, sele)
        self._eval_first(m)
        m = menu.mesh_hide(cmd, sele)
        self._eval_first(m)
        m = menu.mesh_show(cmd, sele)
        self._eval_first(m)
        self._testMeshColor(sele)

    @testing.requires_version('1.8.4')
    def _testMeshColor(self, sele):
        m = menu.mesh_color(cmd, sele)
        self._eval_first(m)

    def testMolAction(self):
        self._make_objects()
        sele = 'm1'
        m = menu.mol_action(cmd, sele)
        self._eval_first(m)

    def testMouseConfig(self):
        m = menu.mouse_config(cmd)
        self._eval_first(m)

    def testObjMotion(self):
        self._make_objects(1)
        sele = 'm1'
        m = menu.obj_motion(cmd, sele)

    def testPickMenu(self):
        self._make_objects(1)
        m = menu.pick_menu(cmd, '/m1////', '(first m1)')

    def testPickSele(self):
        self._make_seles(1)
        sele = 's1'
        m = menu.pick_sele(cmd, sele, 'dummytitle')

    def testPresets(self):
        sele = 'm1'
        cmd.load(self.datafile('1rx1.pdb'), sele)
        m = menu.presets(cmd, sele)
        d = self._to_dict(m)
        self._eval_first(d['simple'])

    def testRampAction(self):
        sele = 'r1'
        cmd.ramp_new(sele, 'none')
        self._testRampColor(sele)
        m = menu.ramp_action(cmd, sele)
        self._eval_first(m)

    @testing.requires_version('1.8')
    def _testRampColor(self, sele):
        m = menu.ramp_color(cmd, sele)
        self._eval_first(m)

    def testSceneMenu(self):
        self._make_scenes()
        m = menu.scene_menu(cmd, '001')
        d = self._to_dict(m)
        self._eval_first(d['update'])

    def testSeleAction(self):
        self._make_seles()
        sele = 's1'
        m = menu.sele_action(cmd, sele)
        d = self._to_dict(m)
        self._eval_first(d['zoom'])

    def testSeqOption(self):
        self._make_objects(1)
        sele = '_seeker'
        cmd.select(sele, 'm1')
        m = menu.seq_option(cmd, sele, '/m1////strip')
        self._eval_first(m)

    def testSimpleAction(self):
        sele = 'cgo1'
        self._make_cgo(sele)
        m = menu.simple_action(cmd, sele)

    def testSliceAction(self):
        sele = 'slice1'
        mapname = 'map1'
        self._make_map(mapname)
        cmd.slice_new(sele, mapname)
        cmd.ramp_new('r2', mapname)
        m = menu.slice_color(cmd, sele)
        self._eval_first(m)
        m = menu.slice_hide(cmd, sele)
        self._eval_first(m)
        m = menu.slice_show(cmd, sele)
        self._eval_first(m)
        m = menu.slice_action(cmd, sele)
        self._eval_first(m)

    def testSurfaceAction(self):
        sele = 'surf1'
        mapname = 'map1'
        self._make_map(mapname)
        cmd.isosurface(sele, mapname)
        m = menu.surface_action(cmd, sele)
        self._eval_first(m)
        m = menu.surface_hide(cmd, sele)
        self._eval_first(m)
        m = menu.surface_show(cmd, sele)
        self._eval_first(m)

    def testVacuum(self):
        cmd.fab('ACD', 'm1')
        m = menu.vacuum(cmd, 'm1')
        self._eval_first(m)

    @testing.requires_version('1.7.2')
    def testVolColor(self):
        sele = 'vol1'
        mapname = 'map1'
        self._make_map(mapname)
        cmd.volume(sele, mapname)
        m = menu.vol_color(cmd, sele)
        # self._eval_first(m)
        m = menu.volume_hide(cmd, sele)
        self._eval_first(m)
        m = menu.volume_show(cmd, sele)
        self._eval_first(m)
