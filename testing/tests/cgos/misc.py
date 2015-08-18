from pymol import cmd, testing, cgo
import unittest

class Test(testing.PyMOLTestCase):

    def test_pse(self):
        self.ambientOnly()

        viewport = (50, 50)
        cmd.viewport(*viewport)

        objname = 'simpletri'
        obj = []
        obj.extend([
            cgo.BEGIN, cgo.TRIANGLES,
            cgo.COLOR, 1.0, 0., 0.,
            cgo.VERTEX, 0.0, 0.0, 1.0,
            cgo.VERTEX, 1.0, 0.0, 0.0,
            cgo.VERTEX, 0.0, 1.0, 0.0,
            cgo.END
            ])
        cmd.load_cgo(obj, objname, 1)

        s = cmd.get_session()
        cmd.set_session(s)

        self.assertEqual([objname], cmd.get_names())
        self.assertImageHasColor('red')
