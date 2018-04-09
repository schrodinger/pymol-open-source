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
