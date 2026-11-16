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
        # 15: ordered, one solvent row instead of twin rows, glass
        labels = [row[1] for row in m if row[0] == 1]
        self.assertEqual(labels, [
            'default', 'classified', 'simple', 'technical', 'ball and stick',
            'b factor putty', 'pretty', 'publication', 'glass', 'wire',
            'illustration', 'ligands',
            'ligand sites', 'protein interface', 'cartoon transparency',
            'show solvent'])
        from pymol import preset
        self.addCleanup(setattr, preset, 'solvent', False)
        cmd.do(d['show solvent'], echo=0)
        self.assertEqual(menu.presets(cmd, sele)[-1][1], 'hide solvent')
        cmd.do(d['glass'], echo=0)
        self.assertGreater(cmd.count_atoms('rep surface'), 0)
        cmd.do(self._to_dict(menu.looks(cmd))['comic'], echo=0)
        self.assertEqual(cmd.get_setting_int('ray_trace_mode'), 3)
        # cartoon transparency lasts until the next preset
        cmd.do(self._to_dict(d['cartoon transparency'])['50%'], echo=0)
        self.assertEqual(cmd.get_setting_float('cartoon_transparency', sele), 0.5)
        cmd.do(d['pretty'], echo=0)
        self.assertEqual(cmd.get_setting_float('cartoon_transparency', sele), 0)

    def testShades(self):
        # "shades" under spectrum and under by chain: one color, light to dark
        cmd.load(self.datafile('1oky-frag.pdb'), 'm1')
        for builder in (menu.spectrum, menu.by_chain):
            rows = self._to_dict(self._to_dict(builder(cmd, 'm1'))['shades'])
            self.assertEqual(len(rows), 8)
            cmd.color('white')
            cmd.do(rows['\\059marine'], echo=0)
            self.assertEqual(cmd.count_atoms('color white'), 0)

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
        self._eval_first(m)

    def testConsolidated(self):
        # 14c: one builder for objects and groups, one for selections
        self._make_seles(1)
        cmd.group('g1', 'm1')
        for name in ('sele_action2', 'test1', 'test2', 'mol_toggle',
                     'simple_show', 'simple_hide'):
            self.assertFalse(hasattr(menu, name), name)

        def labels(m):
            return [row[1] for row in m if row[0] == 1]

        def commands(m):
            return [row[2] for row in m if row[0] == 1 and isinstance(row[2], str)]

        mol = labels(menu.mol_action(cmd, 'm1'))
        self.assertEqual(mol[:2], ['rename object', 'duplicate object'])
        self.assertEqual(len(mol), 10)
        group = labels(menu.group_action(cmd, 'g1'))
        self.assertEqual(group[:2], ['rename group', 'duplicate group'])
        self.assertEqual(len(group), 13)
        more = self._to_dict(menu.mol_action(cmd, 'm1'))['more']
        self.assertIn('copy to object', labels(more))
        self.assertNotIn('copy to object', labels(self._to_dict(
            menu.group_action(cmd, 'g1'))['more']))
        full = labels(menu.sele_action(cmd, 's1'))
        picked = labels(menu.sele_action(cmd, 's1', 1))
        self.assertEqual(set(full) - set(picked), {'drag coordinates', 'clean'})
        # no menu has two rows with the same command
        for m in (menu.mol_action(cmd, 'm1'), more, menu.group_action(cmd, 'g1'),
                  menu.sele_action(cmd, 's1'), menu.sele_action(cmd, 's1', 1),
                  menu.exclude(cmd, 's1')):
            self.assertEqual(len(commands(m)), len(set(commands(m))))
        cmd.do(self._to_dict(menu.mol_action(cmd, 'm1'))['duplicate object'], echo=0)
        self.assertEqual(cmd.count_atoms('m1_cpy'), cmd.count_atoms('m1'))
        cmd.do(self._to_dict(menu.mol_action(cmd, 'm1'))['duplicate object'], echo=0)
        self.assertIn('m1_cpy01', cmd.get_names())
        # every direct command of these menus is valid Python
        for m in (menu.mol_action(cmd, 'm1'), more, menu.group_action(cmd, 'g1'),
                  menu.sele_action(cmd, 's1')):
            for command in commands(m):
                compile(command, 'menu', 'exec')

    def testExclude(self):
        cmd.fab('AG', 'm1')
        cmd.pseudoatom('m1', resn='ZN', elem='Zn', hetatm=1)
        for label, left in (('polymer', 'not polymer'), ('inorganic', 'not inorganic')):
            cmd.select('s1', 'm1')
            cmd.do(self._to_dict(menu.exclude(cmd, 's1'))[label], echo=0)
            self.assertEqual(cmd.count_atoms('s1'), cmd.count_atoms('m1 and ' + left))
            self.assertLess(cmd.count_atoms('s1'), cmd.count_atoms('m1'))
            self.assertGreater(cmd.count_atoms('s1'), 0)

    def testNoViewRows(self):
        # zoom, orient, center, origin: toolbar and double click, not menus
        self._make_seles(1)
        cmd.group('g1', 'm1')
        for m in [
            menu.all_action(cmd, 'all'),
            menu.mol_action(cmd, 'm1'),
            menu.group_action(cmd, 'g1'),
            menu.sele_action(cmd, 's1'),
            menu.main_menu(cmd, (0., 0., 0.), (0, 0)),
            menu.all_option(cmd, 'all'),
            menu.pick_sele(cmd, 's1', ''),
            menu.pick_option(cmd, 'm1', ''),
            menu.seq_option(cmd, 'm1', '/m1/'),
        ]:
            labels = {row[1].split(' ')[0] for row in m}
            self.assertFalse(labels & {'zoom', 'orient', 'center', 'origin'})

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
