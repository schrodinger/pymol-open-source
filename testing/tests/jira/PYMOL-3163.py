'''
Insertion codes CMS import
'''

from pymol import cmd, testing

@testing.requires_version('2.3')
class TestCMSinscodes(testing.PyMOLTestCase):
    def test(self):
        cmd.load(self.datafile('inscodes.cms'))
        self.assertEqual(29, cmd.count_atoms())
        self.assertEqual(7, cmd.count_atoms('resi 52A'))
