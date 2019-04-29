'''
MAE export chain=" "
'''

from pymol import cmd, testing

@testing.requires('incentive')
class Test3231(testing.PyMOLTestCase):

    def _test(self, literal):
        cmd.fragment('ala')
        cmd.alter('all', 'chain=" "') # space
        cmd.alter('index 1', 'chain=""') # empty
        cmd.alter('index 2', "chain='\"'") # quote
        cmd.alter('elem O', "name='O\"'") # quote
        with testing.mktemp('.mae') as filename:
            cmd.save(filename)
            cmd.delete('*')

            if literal:
                cmd.set('pdb_literal_names')
                strip = lambda s: s
            else:
                strip = str.strip

            cmd.load(filename)
            self.assertEqual(10, cmd.count_atoms())

            atomnames = []
            cmd.iterate('all', 'atomnames.append(name)', space=locals())
            self.assertTrue(strip(' O" ') in atomnames)
            self.assertTrue(strip(' N  ') in atomnames)
            self.assertTrue(strip('1HB ') in atomnames)

    @testing.requires_version('2.3.2')
    def test(self):
        self._test(False)

    @testing.requires_version('2.4')
    def test_literal(self):
        self._test(True)
