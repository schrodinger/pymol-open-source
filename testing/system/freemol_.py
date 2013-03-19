'''
FreeMOL tests
'''

from pymol import testing

class TestFreemol(testing.PyMOLTestCase):

    def testValidate(self):
        import freemol.apbs
        import freemol.mengine
        import freemol.mpeg_encode
        self.assertTrue(freemol.apbs.validate())
        self.assertTrue(freemol.mengine.validate())
        self.assertTrue(freemol.mpeg_encode.validate())
