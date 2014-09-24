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
