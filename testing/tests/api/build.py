"""
Testing PyQt Fab & Builder
"""
from pymol import cmd, testing

@testing.requires_version('2.3')
class TestNucBuilder(testing.PyMOLTestCase):
    def testCanInit(self):
        self.assertEqual(cmd.fnab("A"), None)

    def testSingleDNAFASTA(self):
        cmd.fnab(input="ATTG", type="RNA", form="B", dbl_helix=-1)
        self.assertEqual("ATTG", "ATTG")
    def testSingleRNAFASTA(self):
        cmd.fnab(input="AUTG", type="RNA")

    def testDoubleDNAFASTA(self):
        dna = "ATCCCCG"
        revcomp = "CGGGGAT"
        cmd.fnab(input=dna, dbl_helix=1)
        fasta_str = cmd.get_fastastr().splitlines()
        fasta_raw = (fasta_str[1], fasta_str[3])
        sense_strand = fasta_raw[0]
        antisense_strand = fasta_raw[1]
        self.assertEqual(dna, sense_strand)
        self.assertEqual(revcomp, antisense_strand)
        cmd.delete('all')
        cmd.fnab(input=dna)
        fasta_str = cmd.get_fastastr().splitlines()
        fasta_raw = (fasta_str[1], fasta_str[3])
        sense_strand = fasta_raw[0]
        antisense_strand = fasta_raw[1]
        self.assertEqual(dna, sense_strand)
        self.assertEqual(revcomp, antisense_strand)

    def testDoubleRNAFASTA(self):
        rna = "AUUUUUUUCG"
        cmd.fnab(input=rna, type="RNA", form="B", dbl_helix=1)

        fasta_str = cmd.get_fastastr().splitlines()
        self.assertEqual(len(fasta_str), 2)        
        sense_strand = fasta_str[1]
        self.assertEqual(rna, sense_strand)

    def testSkipBadResidues(self):
        dna = "ATTTTZTTGCCCGGXXG"
        cmd.fnab(input=dna)
        fasta_str = cmd.get_fastastr().splitlines()
        fasta_raw = (fasta_str[1], fasta_str[3])
        sense_strand = fasta_str[1]
        self.assertEqual("ATTTTTTGCCCGGG", sense_strand)

    def testMixedCase(self):
        dna = "AtG"
        revcomp = "CAT"
        cmd.fnab(input=dna)

        fasta_str = cmd.get_fastastr().splitlines()
        fasta_raw = (fasta_str[1], fasta_str[3])
        sense_strand = fasta_raw[0]
        antisense_strand = fasta_raw[1]
        self.assertEqual(dna.upper(), sense_strand)
        self.assertEqual(revcomp, antisense_strand)
