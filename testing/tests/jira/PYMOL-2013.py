'''
new feature: cartoon_gap_cutoff

Not very sophisticated test, just check if images with and without
setting differ.
'''

from pymol import cmd, CmdException, testing, stored

class TestPYMOL2013(testing.PyMOLTestCase):

    def test(self):
        self.ambientOnly()

        cmd.viewport(150, 150)
        cmd.load(self.datafile('1oky-frag.pdb'))
        cmd.cartoon('dash')
        cmd.show_as('cartoon')
        cmd.orient()

        # no gaps
        img_nogaps = self.get_imagearray()

        # gap of length 0
        cmd.unbond('96/C', '97/N')
        img_gap = self.get_imagearray()
        self.assertFalse((img_gap == img_nogaps).all())

        # close gap with setting, should match first image (special case)
        cmd.set('cartoon_gap_cutoff', 1) # exact cutoff
        img_gap = self.get_imagearray()
        self.assertTrue((img_gap == img_nogaps).all())

        # gap of length 1
        cmd.remove('92/')
        img_gap = self.get_imagearray()
        cmd.set('cartoon_gap_cutoff', 2) # exact cutoff
        img_nogaps = self.get_imagearray()
        self.assertFalse((img_gap == img_nogaps).all())

        # gap of length 3
        cmd.remove('102-104/')
        cmd.cartoon('auto') # default
        img_gap = self.get_imagearray()
        cmd.set('cartoon_gap_cutoff', 10) # arbitrary larger cutoff
        img_nogaps = self.get_imagearray()
        self.assertFalse((img_gap == img_nogaps).all())
