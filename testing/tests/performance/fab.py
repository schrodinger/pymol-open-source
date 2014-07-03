from pymol import cmd, testing, stored

class StressFab(testing.PyMOLTestCase):

    def test(self):
        with self.timing(max=0.1):
            cmd.fab('AAAA')
