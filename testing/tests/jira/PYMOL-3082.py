'''
chempy.Atom defaults
'''

import os
import sys
from pymol import cmd, CmdException, testing, stored

@testing.requires_version('1.8.2')
class TestChempyAtom(testing.PyMOLTestCase):
    def _assert_iterate_equal(self, sele, key, value):
        stored.values = []
        cmd.iterate(sele, 'stored.values.append(' + key + ')')
        for v in stored.values:
            self.assertEqual(v, value)

    def test(self):
        cmd.fragment('gly', 'm1')

        cmd.alter('*', 'chain = "ABC"')
        cmd.alter('*', 'segi = "DEF"')

        m = cmd.get_model()

        self.assertEqual(m.atom[0].resn, "GLY")
        self.assertEqual(m.atom[0].chain, "ABC")
        self.assertEqual(m.atom[0].segi, "DEF")

        cmd.alter('*', 'resn = ""')
        cmd.alter('*', 'chain = ""')
        cmd.alter('*', 'segi = ""')
        cmd.alter('*', 'resi = "0"')

        m = cmd.get_model()

        self.assertEqual(m.atom[0].resn, "")
        self.assertEqual(m.atom[0].chain, "")
        self.assertEqual(m.atom[0].segi, "")
        self.assertEqual(m.atom[0].resi, "0")
        self.assertEqual(m.atom[0].resi_number, 0)

        cmd.alter('*', 'resn = "UNK"')
        cmd.alter('*', 'resi = "1A"')

        m = cmd.get_model()

        self.assertEqual(m.atom[0].resn, "UNK")
        self.assertEqual(m.atom[0].resi, "1A")
        self.assertEqual(m.atom[0].resi_number, 1)

        if m.atom[0].has('ins_code'):  # not exported before 2.3
            self.assertEqual(m.atom[0].ins_code, "A")

        cmd.delete('*')
        cmd.load_model(m, 'm1')

        self._assert_iterate_equal('*', 'resi', '1A')
        self._assert_iterate_equal('*', 'resn', 'UNK')

        for a in m.atom:
            a.ins_code = "B"

        cmd.delete('*')
        cmd.load_model(m, 'm1')

        self._assert_iterate_equal('*', 'resi', '1B')

        cmd.alter('*', 'resi = "2"')

        m = cmd.get_model()

        self.assertEqual(m.atom[0].resi, "2")

        if m.atom[0].has('ins_code'):  # not exported before 2.3
            self.assertEqual(m.atom[0].ins_code, "")

        for a in m.atom:
            del a.resi_number  # pre 2.3
            a.resi = "3C"

            # >= 2.3 sets resi_number for resi
            if a.resi_number == 3:
                self.assertEqual(a.ins_code, 'C')

        cmd.delete('*')
        cmd.load_model(m, 'm1')

        self._assert_iterate_equal('*', 'resi', '3C')
