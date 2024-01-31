'''
ray+png in batch mode did ray trace twice
'''

from pymol import cmd, testing, stored

class Test365(testing.PyMOLTestCase):

    def test(self):
        cmd.fragment('ala')
        cmd.ray(100, 100)
        img = self.get_imagearray() # cmd.png without optional arguments

        # bug was: ray tracing twice, second time would not be 100x100
        self.assertEqual(img.shape[:2], (100,100))
