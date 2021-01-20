'''
Bonds to symmetry mates
'''

from pymol import cmd, testing


@testing.requires_version('2.5')
class TestBondSymOp(testing.PyMOLTestCase):
    def test_commands(self):
        sym = [5, 2, 3, 60, 90, 90, 'P M 1 1']
        cmd.set('orthoscopic')
        cmd.set('stick_ball')
        cmd.set('stick_ball_color', 'blue')
        cmd.pseudoatom('m1', name='A1', pos=(0.5, 0, 0))
        cmd.pseudoatom('m1', name='A2', pos=(2, 2, 0.5))
        cmd.pseudoatom('m1', name='A3', pos=(2, 2, 2.0))
        cmd.set_symmetry('m1', *sym)
        cmd.bond('name A1', 'name A1', symop="2_555")
        cmd.bond('name A2', 'name A3')
        cmd.valence(1, 'name A2', '*', symop="1_564")
        cmd.color('red')
        cmd.color('0x00cc00', 'name A2')
        cmd.color('yellow', 'name A3')
        cmd.show_as('sticks')
        cmd.viewport(200, 200)
        cmd.set_view((0., 1., 0., 0., 0., 1., 1., 0., 0., -0.223610640,
                      -0.169926167, -10.333316803, 1., 1., 1., 6., 13., 20.))
        self.ambientOnly()
        self.assertImageEqual("symop-ref/bond_symop.png")

        m = cmd.get_model()
        self.assertEqual(m.bond[0].symmetry_2, "2_555")
        self.assertEqual(m.bond[1].symmetry_2, "1_564")

        # load model
        cmd.delete('*')
        cmd.load_model(m, 'm2', zoom=0)
        cmd.set_symmetry('m2', *sym)
        cmd.color('red')
        cmd.color('0x00cc00', 'name A2')
        cmd.color('yellow', 'name A3')
        cmd.show_as('sticks')
        self.assertImageEqual("symop-ref/bond_symop.png")

        # load session
        for binary in (0, 1):
            cmd.set_session(cmd.get_session(binary=binary))
            self.assertImageEqual("symop-ref/bond_symop.png")

    @testing.requires("incentive")
    def test_load_mae(self):
        cmd.load(self.datafile("pbc.mae"))
        cmd.viewport(200, 150)
        cmd.set_view(
            (0.993128419, -0.107736677, -0.045689747, 0.110594988, 0.736444175,
             0.667396009, -0.038255211, -0.667863667, 0.743299425,
             -0.000001433, -0.000000438, -8.649258614, 0.962046325,
             0.647872984, 0.307970583, 5.649246216, 11.649243355, -20.))
        self.ambientOnly()
        self.assertImageEqual("symop-ref/pbc.png")
