import os
import sys
from pymol import cmd, testing, stored

@testing.requires_version('2.3')
class TestSeqalign(testing.PyMOLTestCase):

    def testLoadAlnMatchingIds(self):
        cmd.fab('ACDEFGHIKLMNPQRS', 'seq1')
        cmd.fab('ACDIKLMNP', 'seq2')
        cmd.fab('GHIKPQRS', 'seq3')
        cmd.load(self.datafile('alignment.aln'), 'aln')
        self.assertEqual(cmd.count_atoms('guide & aln & seq1'), 11)
        self.assertEqual(cmd.count_atoms('guide & aln & seq2'), 7)
        self.assertEqual(cmd.count_atoms('guide & aln & seq3'), 6)

    def testLoadAlnMappingStr(self):
        cmd.fab('ACDEFGHIKLMNPQRS', 'm1')
        cmd.fab('ACDIKLMNP', 'm2')

        from pymol.seqalign import load_aln_multi

        load_aln_multi(self.datafile('alignment.aln'), 'aln', mapping=
            'seq1 m1 '
            'seq2 m2 '
        )

        self.assertEqual(cmd.count_atoms('guide & aln & m1'), 7)
        self.assertEqual(cmd.count_atoms('guide & aln & m2'), 7)

    def testLoadAlnMappingDict(self):
        cmd.fab('ACDIKLMNP', 'm2')
        cmd.fab('GHIKPQRS', 'm3')

        from pymol.seqalign import load_aln_multi

        load_aln_multi(self.datafile('alignment.aln'), 'aln', mapping={
            'seq2': 'm2',
            'seq3': 'm3',
        })

        self.assertEqual(cmd.count_atoms('guide & aln & m2'), 2)
        self.assertEqual(cmd.count_atoms('guide & aln & m3'), 2)
