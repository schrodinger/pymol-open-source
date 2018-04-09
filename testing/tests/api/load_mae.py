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
