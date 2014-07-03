from pymol import cmd, CmdException, testing, stored

class TestPYMOL1293(testing.PyMOLTestCase):

    def test(self):
        cmd.pseudoatom('m1')
        cmd.mset('1x1')
        cmd.create('m2', 'm1')
        cmd.ray()
        v = cmd.get_object_matrix('m2')
        self.assertArrayEqual(v, [
            1., 0., 0., 0.,
            0., 1., 0., 0.,
            0., 0., 1., 0.,
            0., 0., 0., 1.,
            ], 0.001, 'object matrix not identity')
