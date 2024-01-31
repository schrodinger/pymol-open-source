'''
FreeMOL tests
'''

from pymol import testing

@testing.requires('freemol')
class TestFreemol(testing.PyMOLTestCase):

    def testValidate(self):
        import freemol.apbs
        import freemol.mengine
        import freemol.mpeg_encode

        # psize.py not used with PyMOL 2.0
        freemol.apbs._psize_py = '<ignore>'

        self.assertTrue(freemol.apbs.validate())
        self.assertTrue(freemol.mengine.validate())
        self.assertTrue(freemol.mpeg_encode.validate())
