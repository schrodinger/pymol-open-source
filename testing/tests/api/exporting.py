'''
unit tests for pymol.exporting
'''

import os
import tempfile
import Image

import pymol.exporting
from pymol import cmd, testing, stored

class TestExporting(testing.PyMOLTestCase):

    def testCache(self):
        for action in pymol.exporting.cache_action_dict:
            cmd.cache(action)

    def testCopyImage(self):
        cmd.copy_image
        self.skipTest("TODO")

    def testExportCoords(self):
        cmd.fragment('ala')
        c = cmd.export_coords('ala', 1)
        self.assertTrue(bool(c))
        # TODO: pymol.experimenting.load_coords is crude and buggy, for
        # example the state is decremented twice.

    def testGetFastastr(self):
        seq, name = 'ACD', 'm1'
        cmd.fab(seq, name)
        s = cmd.get_fastastr()
        self.assertEqual(s.split(), ['>' + name, seq])

    def testGetPdbstr(self):
        cmd.pseudoatom()
        lines = cmd.get_pdbstr().splitlines()
        self.assertTrue(lines[0].startswith('HETATM'))
        self.assertTrue(lines[1].startswith('END'))

    def testGetSession(self):
        cmd.fragment('ala')
        x = cmd.count_atoms()
        s = cmd.get_session()
        cmd.reinitialize()
        cmd.set_session(s)
        self.assertEqual(['ala'], cmd.get_names())
        self.assertEqual(x, cmd.count_atoms())

    def testMultisave(self):
        cmd.multisave
        self.skipTest("TODO")

    @testing.foreach(
        ('0.5in', '0.25in', 200, (100, 50)),   # inch
        ('2.54cm', '1.27cm', 100, (100, 50)),  # centimeter
        (100, 0, -1, (100, 75)),               # px, only width
        (0, 75, -1, (100, 75)),                # px, only height
    )
    def testPng(self, w, h, dpi, size):
        with testing.mktemp('.png') as filename:
            cmd.png(filename, w, h, dpi, ray=1)
            self.assertEqual(size, Image.open(filename).size)

    def testSave(self):
        cmd.fragment('ala')
        with testing.mktemp('.pdb') as filename:
            cmd.save(filename)
            filesize = os.path.getsize(filename)
            self.assertTrue(filesize > 0)

