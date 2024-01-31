"""
Testing PyQt Fab & Builder
"""
from pymol import cmd, testing

@testing.requires('incentive')
@testing.requires_version('2.3')
class TestNucBuilder(testing.PyMOLTestCase):
    def testCanInit(self):
        self.assertEqual(cmd.fnab("A"), None)

    def testSingleDNAFASTA(self):
        dna = "ATGC"
        cmd.fnab(input=dna, mode="DNA", form="B", dbl_helix=-1)
        fasta_str = cmd.get_fastastr().splitlines()
        self.assertEqual(dna, fasta_str[1])

    def testSingleRNAFASTA(self):
        rna = "AUGC"
        cmd.fnab(input=rna, mode="RNA")
        fasta_str = cmd.get_fastastr().splitlines()
        self.assertEqual(rna, fasta_str[1])

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
        cmd.fnab(input=rna, mode="RNA", form="B", dbl_helix=1)

        fasta_str = cmd.get_fastastr().splitlines()
        self.assertEqual(len(fasta_str), 2)        
        sense_strand = fasta_str[1]
        self.assertEqual(rna, sense_strand)

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
