'''
unit tests for pymol.exporting
'''

import os
import sys
import tempfile
import Image
import unittest

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
        if sys.version_info.major > 2:
            self.skipTest("not Python3 compatible")
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
        names = ['ala', 'gly', 'his', 'arg']

        for name in names:
            cmd.fragment(name)

        with testing.mktemp('.pdb') as filename:
            cmd.multisave(filename, names[0])                       # new
            cmd.multisave(filename, names[1])                       # new (overwrite)
            cmd.multisave(filename, ' '.join(names[2:]), append=1)  # append
            cmd.delete('*')
            cmd.load(filename)

        self.assertEqual(cmd.get_object_list(), names[1:])

    @testing.foreach(
        ('0.5in', '0.25in', 200, (100, 50)),   # inch
        ('2.54cm', '1.27cm', 100, (100, 50)),  # centimeter
        (100, 0, -1, (100, 75)),               # px, only width
        (0, 75, -1, (100, 75)),                # px, only height
    )
    @testing.requires('no_edu') # ray
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

    def testSaveRef(self):
        cmd.fragment('ala', 'm1')
        cmd.copy('m2', 'm1')
        cmd.copy('m3', 'm1')

        cmd.rotate('y', 90, 'm2')
        cmd.align('m3', 'm2')

        # with ref=m3
        with testing.mktemp('.pdb') as filename:
            cmd.save(filename, 'm2', ref='m3')
            cmd.load(filename, 'm4')
            self.assertAlmostEqual(cmd.rms_cur('m4', 'm1'), 0.00, delta=1e-2)
            self.assertAlmostEqual(cmd.rms_cur('m4', 'm2'), 1.87, delta=1e-2)

        # without ref
        with testing.mktemp('.pdb') as filename:
            cmd.save(filename, 'm2')
            cmd.load(filename, 'm5')
            self.assertAlmostEqual(cmd.rms_cur('m5', 'm2'), 0.00, delta=1e-2)
            self.assertAlmostEqual(cmd.rms_cur('m5', 'm1'), 1.87, delta=1e-2)

    # cmp_atom : compares all fields in Atom (see chempy/__init__.py)
    #            except the id (which is unique to the instance)
    def cmp_atom(self, selfobj,other):
        cmp = lambda a, b: a != b # py3k
        return \
                cmp(type(selfobj), type(other)) or \
                cmp(selfobj.segi, other.segi) or \
                cmp(selfobj.chain, other.chain) or \
                cmp(selfobj.resi_number, other.resi_number) or \
                cmp(selfobj.resi, other.resi) or \
                cmp(selfobj.resn, other.resn) or \
                cmp(selfobj.symbol, other.symbol) or \
                cmp(selfobj.name, other.name)

    # cmp_bond : compares all fields in Bond (see chempy/__init__.py)
    def cmp_bond(self, selfobj,other):
        cmp = lambda a, b: a != b # py3k
        return \
                cmp(selfobj.order, other.order) or \
                cmp(selfobj.stereo, other.stereo)

    def assertModelsAreSame(self, m1, m2):
        self.assertTrue(len(m1.atom) == len(m2.atom))
        idx = 0
        for m1atomidx in m1.atom:
            self.assertTrue(self.cmp_atom(m1atomidx, m2.atom[idx]) == 0)
            idx = idx + 1
        idx = 0
        for m1bondidx in m1.bond:
            self.assertTrue(self.cmp_bond(m1bondidx, m2.bond[idx]) == 0)
            idx = idx + 1

    @testing.requires_version('1.7.6')
    @unittest.skipIf(sys.version_info[0] > 2, 'pse_binary_dump not py3k ready')
    def testPSEBulkImport(self):
        cmd.load(self.datafile('1rx1_1766_bulk.pse.gz'))
        m1 = cmd.get_model()
        cmd.load(self.datafile('1rx1_176.pse.gz'))
        m2 = cmd.get_model()
        self.assertModelsAreSame(m1, m2)

    @testing.foreach.product((0, 1.7, 1.76, 1.8, 1.82, 1.9), (0, 1))
    @testing.requires_version('1.7.6.5')
    def testPSEBulkExportImport(self, pse_export_version, pse_binary_dump):
        with testing.mktemp('.pse') as filename:
            cmd.load(self.datafile("1oky-frag.pdb"))
            m1 = cmd.get_model()
            cmd.set("pse_export_version", pse_export_version)
            cmd.set("pse_binary_dump", pse_binary_dump)
            cmd.save(filename)
            cmd.reinitialize()
            cmd.load(filename)
            m2 = cmd.get_model()
            self.assertModelsAreSame(m1, m2)
