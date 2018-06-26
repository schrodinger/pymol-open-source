'''
unit tests for pymol.wizard.nucmutagenesis
'''

import os
import sys
import tempfile
import unittest

from pymol import cmd, testing

class TestNucMutagenesis(testing.PyMOLTestCase):

    def test_CanGetSourceSequence(self):
        cmd.load(self.datafile("1rna.cif"))
        seq = cmd.get_fastastr('/1rna/A/A').splitlines()[1]
        self.assertEqual(seq, "UUAUAUAUAUAUAA")

    def test_CanGetFragment(self):
        cmd.fragment("gtp")
        numAtoms = cmd.count_atoms("gtp")
        self.assertEqual(numAtoms, 48)

    def test_CanInit(self):
        cmd.wizard("nucmutagenesis")
        self.assertTrue(cmd.get_wizard() != None)

    def test_CorrectChiDihedral(self):
        cmd.load(self.datafile("1rna.cif"))
        src_dihedral = cmd.get_dihedral("/1rna/A/A/14 & name O4'",
                                        "/1rna/A/A/14 & name C1'",
                                        "/1rna/A/A/14 & name N9",
                                        "/1rna/A/A/14 & name C4'")
        cmd.wizard("nucmutagenesis")
        cmd.select("/1rna/A/A/14")
        cmd.get_wizard().mode = "(d)G"
        cmd.get_wizard().do_select("sele")
        cmd.get_wizard().apply()
        des_dihedral = cmd.get_dihedral("/1rna/A/A/14 & name O4'",
                                        "/1rna/A/A/14 & name C1'",
                                        "/1rna/A/A/14 & name N9",
                                        "/1rna/A/A/14 & name C4'")
        self.assertAlmostEqual(src_dihedral, des_dihedral, delta = 2.0)

    def test_CanGetNewRNASequence(self):
        cmd.load(self.datafile("1rna.cif"))
        cmd.wizard("nucmutagenesis")
        cmd.select("/1rna/A/A/14")
        cmd.get_wizard().mode = "(d)G"
        cmd.get_wizard().do_select("sele")
        cmd.get_wizard().apply()
        seq = cmd.get_fastastr("/1rna/A/A").splitlines()[1]
        self.assertEqual(seq, "UUAUAUAUAUAUAG")

    def test_CanGetNewDNASequence(self):
        cmd.load(self.datafile("1bna.cif"))
        cmd.wizard("nucmutagenesis")
        cmd.select("/1bna/A/A/1")
        cmd.get_wizard().mode = "(d)A"
        cmd.get_wizard().do_select("sele")
        cmd.get_wizard().apply()
        seq = cmd.get_fastastr("/1bna/A/A").splitlines()[1]
        self.assertEqual(seq, "AGCGAATTCGCG")

    def test_CanMutateNonCanonicalNucleo(self):
        cmd.load(self.datafile("1k5e.cif"))
        cmd.wizard("nucmutagenesis")
        cmd.select("/1k5e/A/A/6")
        cmd.get_wizard().mode = "(d)A"
        cmd.get_wizard().do_select("sele")
        cmd.get_wizard().apply()
        seq = cmd.get_fastastr("/1k5e/A/A").splitlines()[1]
        self.assertEqual(seq, "CGGACAAGAAG")