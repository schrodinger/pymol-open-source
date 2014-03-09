'''
Crash with sticks, valence and order=0
'''

from pymol import cmd, testing, stored

class Test1916(testing.PyMOLTestCase):

    def test(self):
        from chempy import Atom, Bond, models

        m = models.Indexed()

        for i in range(2):
            a = Atom()
            a.coord = [float(i), 0., 0.]
            m.add_atom(a)

        b = Bond()
        b.index = [0, 1]
        b.order = 0
        m.add_bond(b1)

        cmd.load_model(m, 'foo')

        cmd.set('valence')
        cmd.show('sticks')
