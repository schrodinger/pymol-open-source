'''
Test loading MAE files
'''

from pymol import cmd, testing

mae_filenames = [
    'data/foo.mae',
    'data/foo.mae.gz',
    'data/foo.cms',
#    'data/foo.cms.gz',
]

mae_urls = [
    'http://thomas-holder.de/tmp/foo.mae',
    'http://thomas-holder.de/tmp/foo.mae.gz',
]

@testing.requires('incentive')
class TestLoadMAE(testing.PyMOLTestCase):

    @testing.foreach.zip(mae_filenames)
    def testLoad(self, filename):
        cmd.load(filename)

        v = cmd.get_object_list()
        self.assertEqual(v, ['foo'])

    @testing.foreach.zip(mae_urls)
    @testing.requires('network')
    def testLoadURL(self, filename):
        cmd.load(filename)

        v = cmd.get_object_list()
        self.assertEqual(v, ['foo'])

    @testing.requires_version('2.1.1')
    def testLoadCryst1(self):
        cmd.load(self.datafile('cryst1.mae'), 'm1')
        s = cmd.get_symmetry('m1')
        self.assertArrayEqual(s[:6], [50.84, 42.77, 28.95, 90., 90., 90.], delta=0.01)
        self.assertEqual(s[6], 'P 21 21 21')

    def _assertCountEqual(self, sele1, sele2, delta=0):
        self.assertAlmostEqual(
                cmd.count_atoms(sele1),
                cmd.count_atoms(sele2), delta=delta)

    @testing.requires_version('2.2')
    def testExportStyle(self):
        cmd.fab('ACDEF', 'm1')
        cmd.hide()
        cmd.show('cartoon', 'resi 1-3')
        cmd.show('lines', 'resn CYS')
        cmd.show('sticks', 'resn ASP+PHE')
        cmd.show('spheres', 'resn GLU')
        cmd.set('stick_ball', 1, 'resn PHE')
        cmd.set('stick_ball_ratio', 1.5, 'm1')
        testlabel = 'Hello "World"'
        cmd.label('name SG', repr(testlabel))

        with testing.mktemp('.mae') as filename:
            cmd.save(filename)
            cmd.delete('*')
            cmd.load(filename, 'm2')

        g_labels = []
        cmd.iterate('name SG', 'g_labels.append(label)', space=locals())
        cmd.alter('*', 'b = 1 if s.stick_ball else 0')

        self._assertCountEqual('rep cartoon & guide', 'resi 1-3 & guide')
        self._assertCountEqual('rep lines', 'resn CYS', delta=1)
        self._assertCountEqual('rep sticks', 'resn ASP+PHE')
        self._assertCountEqual('rep spheres', 'resn GLU')
        self._assertCountEqual('b > 0.5', 'resn PHE')
        self.assertTrue(cmd.get_setting_float('stick_ball_ratio', 'm2') > 1.1)
        self.assertEqual(g_labels[0], testlabel)

    @testing.foreach.product((0, 1), (0, 1))
    @testing.requires_version('2.3')
    def testLoadMae(self, multiplex, discrete):
        cmd.load(self.datafile('multimae.maegz'), 'm',
                multiplex=multiplex, discrete=discrete)

        nstate = 4
        natoms = 79
        nmodel = 1
        ndiscrete = discrete

        if multiplex:
            natoms *= nstate
            ndiscrete *= nstate
            nmodel = nstate
            nstate = 1
        elif discrete:
            self.assertEqual(cmd.count_atoms(state=1), natoms)
            natoms *= nstate

        self.assertEqual(cmd.count_states(), nstate)
        self.assertEqual(cmd.count_discrete('*'), ndiscrete)
        self.assertEqual(cmd.count_atoms(), natoms)
        self.assertEqual(len(cmd.get_object_list()), nmodel)

    @testing.requires_version('2.3')
    def testLoadMaeUserLabel(self):
        from pymol import stored
        cmd.load(self.datafile('userlabels2.mae'), 'm')
        cmd.iterate('rank 1', 'stored.label = label')
        self.assertEqual(stored.label, '1.31 TRP')

    @testing.requires_version('2.5')
    def testLabel(self):
        cmd.load(self.datafile('labels.mae'), 'm')
        labels = []
        cmd.iterate('all', 'labels.append(label)', space={'labels': labels})
        self.assertEqual(labels, [
            'N +1 (-1.2, 0.2, -0.2) 0.290 0.75 18.50 GLY 1 D N 1 1',
            'C  (0.2, 0.3, -0.5) -0.010 0.75 24.17 GLY 1 D CA 2 2',
            'C  (1.1, -0.4, 0.5) 0.620 0.75 11.07 GLY 1 D C 3 3',
            'O  (0.5, -1.0, 1.5) -0.570 0.75 31.09 GLY 1 D O 4 4',
            'H  (0.4, -0.1, -1.5) 0.090 0.75 24.17 GLY 1 D HA 9 16',
            'H  (0.5, 1.4, -0.5) 0.090 0.75 24.17 GLY 1 D HA3 8 15',
            'H  (-1.5, 1.0, 0.4) 0.160 0.75 18.50 GLY 1 D H1 7 14',
            'H  (-1.7, 0.2, -1.1) 0.160 0.75 18.50 GLY 1 D H2 5 12',
            'H  (-1.4, -0.7, 0.3) 0.160 0.75 18.50 GLY 1 D H3 6 13',
            'N  (2.3, -0.4, 0.4) -0.380 1.00 48.64 CYS 2 D N 10 5',
            'C  (3.1, -1.1, 1.4) -0.160 1.00 41.94 CYS 2 D CA 11 6',
            'C  (4.6, -0.9, 1.1) 0.750 1.00 0.21 CYS 2 D C 12 7',
            'O  (5.0, -0.3, 0.1) -0.800 1.00 3.55 CYS 2 D O 13 8',
            'C  (2.7, -2.5, 1.5) -0.200 1.00 38.71 CYS 2 D CB 14 9',
            'S  (3.5, -3.4, 2.9) -0.310 1.00 18.13 CYS 2 D SG 15 10',
            'O -1 (5.6, -1.6, 2.0) -0.800 1.00 0.21 CYS 2 D OXT 16 11',
            'H  (2.8, 0.1, -0.3) 0.270 1.00 48.64 CYS 2 D H 21 21',
            'H  (3.0, -0.6, 2.4) 0.140 1.00 41.94 CYS 2 D HA 17 17',
            'H  (2.9, -3.0, 0.6) 0.140 1.00 38.71 CYS 2 D HB2 18 18',
            'H  (1.6, -2.6, 1.7) 0.140 1.00 38.71 CYS 2 D HB3 19 19',
            'H  (3.8, -4.6, 2.5) 0.210 1.00 18.13 CYS 2 D HG 20 20',
        ])
