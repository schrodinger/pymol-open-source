
from pymol import cmd, testing, stored

class TestSelecting(testing.PyMOLTestCase):

    def testDeselect(self):
        cmd.load(self.datafile("1oky.pdb.gz"), "m1")
        cmd.select("all")
        self.assertEquals("sele" in cmd.get_names("all",enabled_only=1), True)
        cmd.deselect()
        self.assertEquals("sele" in cmd.get_names("all",enabled_only=1), False)

    def testIndicate(self):
        cmd.load(self.datafile("1oky.pdb.gz"), "m1")
        cmd.indicate("polymer")
        self.assertEquals(cmd.count_atoms("polymer"),cmd.count_atoms("indicate"))
        self.assertEquals(cmd.count_atoms("indicate & solvent & organic & inorganic"), 0)

    def testPop(self):
        cmd.fragment("ala")
        cmd.select("src", "ala")
        # iterate across the 10 atoms in ala
        cnt = 0
        while cmd.pop("tmp","src"):
            cnt+=1
        self.assertEquals(cnt,10)

    def testSelect(self):
        cmd.select
        self.skipTest("TODO: This will take forever.")

    def testSelectList(self):
        cmd.select_list
        self.skipTest("TODO")

    def testMacros(self):
        '''
        Test selection macros: /model/segi/chain/resn`resi/name`alt

        See also tests/jira/PYMOL-2720.py
        '''
        cmd.load(self.datafile("1oky.pdb.gz"), "m1")
        cmd.alter('*', 'segi="X"')
        cmd.alter('resi 200-238', 'chain="B"')
        cmd.alter('resi 239-360', 'chain="C"')
        cmd.alter('resi 150-220', 'segi="Y"')
        cmd.alter('resi 220-', 'segi="Z"')

        # all atoms
        self.assertEqual(cmd.count_atoms('/'), cmd.count_atoms())
        self.assertEqual(cmd.count_atoms('/////'), cmd.count_atoms())

        # model
        self.assertEqual(cmd.count_atoms('/m1'), cmd.count_atoms('all'))
        self.assertEqual(cmd.count_atoms('m1////'), cmd.count_atoms('all'))

        # segi
        self.assertEqual(cmd.count_atoms('//Y'), cmd.count_atoms('segi Y'))
        self.assertEqual(cmd.count_atoms('//Y+Z'), cmd.count_atoms('segi Y+Z'))

        # chains
        self.assertEqual(cmd.count_atoms('///B'), cmd.count_atoms('chain B'))
        self.assertEqual(cmd.count_atoms('///A+B'), cmd.count_atoms('chain A+B'))

        # resn/resi
        self.assertEqual(cmd.count_atoms('////100'), cmd.count_atoms('resi 100'))
        self.assertEqual(cmd.count_atoms('////`100'), cmd.count_atoms('resi 100'))
        self.assertEqual(cmd.count_atoms('////100-110'), cmd.count_atoms('resi 100-110'))
        self.assertEqual(cmd.count_atoms('////ARG`100'), cmd.count_atoms('resi 100'))
        self.assertEqual(cmd.count_atoms('////ALA'), cmd.count_atoms('resn ALA'), )
        self.assertEqual(cmd.count_atoms('ALA/'), cmd.count_atoms('resn ALA'), )
        self.assertEqual(cmd.count_atoms('////ALA`'), cmd.count_atoms('resn ALA'))
        self.assertEqual(cmd.count_atoms('////ALA+GLU'), cmd.count_atoms('resn ALA+GLU'), )

        # name/alt
        self.assertEqual(cmd.count_atoms('/////CG'), cmd.count_atoms('name CG'))
        self.assertEqual(cmd.count_atoms('*/CG'), cmd.count_atoms('name CG'))
        self.assertEqual(cmd.count_atoms('/////CG`B'), cmd.count_atoms('name CG & alt B'))
        self.assertEqual(cmd.count_atoms('/////`B'), cmd.count_atoms('alt B'))

        # combies
        self.assertEqual(cmd.count_atoms('100/CA'), cmd.count_atoms('resi 100 & name CA'))
        self.assertEqual(cmd.count_atoms('A//CA'), cmd.count_atoms('chain A & name CA'))
        self.assertEqual(cmd.count_atoms('A//`B'), cmd.count_atoms('chain A & alt B'))
        self.assertEqual(cmd.count_atoms('/m1/Y/B/*/CA'), cmd.count_atoms('segi Y & chain B & name CA'))

        # or-logic
        ref = cmd.count_atoms('(resi 100 & name CA) (resi 110 & name CB)')
        self.assertEqual(cmd.count_atoms('100/CA|110/CB'), ref)
        self.assertEqual(cmd.count_atoms('100/CA or 110/CB'), ref)
        self.assertEqual(cmd.count_atoms('100/CA 110/CB'), ref)

    @testing.requires_version('1.8.7')
    def test_not_enabled(self):
        cmd.fragment('ala')
        cmd.fragment('arg')
        cmd.fragment('asp')

        cmd.fragment('gly')
        cmd.fragment('glu')

        cmd.delete('not a*')  # use "not " prefix

        self.assertEqual(cmd.get_names(), ['ala', 'arg', 'asp'])

        cmd.disable('arg')
        cmd.delete('!enabled') # use "!" prefix

        self.assertEqual(cmd.get_names(), ['ala', 'asp'])

        cmd.disable('asp')
        cmd.delete('enabled')

        self.assertEqual(cmd.get_names(), ['asp'])

    @testing.requires_version('2.1')
    def test_protein_nucleic(self):
        npro = 340
        nnuc = 112
        cmd.load(self.datafile("1oky-frag.pdb"), "1p")
        cmd.load(self.datafile('1ehz-5.pdb'), "1n")
        self.assertEqual(cmd.count_atoms('polymer'), npro + nnuc)
        self.assertEqual(cmd.count_atoms('polymer.protein'), npro)
        self.assertEqual(cmd.count_atoms('polymer.nucleic'), nnuc)

    @testing.requires_version('2.2')
    def test_wildcards(self):
        cmd.pseudoatom('foo')
        cmd.pseudoatom('foobar')
        cmd.pseudoatom('foo_bar')
        cmd.pseudoatom('foo_bar_bla')
        cmd.pseudoatom('bar')
        cmd.pseudoatom('_bar')
        self.assertEqual(1, cmd.count_atoms('foo'))
        self.assertEqual(1, cmd.count_atoms('*foo'))
        self.assertEqual(4, cmd.count_atoms('*foo*'))
        self.assertEqual(1, cmd.count_atoms('bar'))
        self.assertEqual(4, cmd.count_atoms('*bar'))
        self.assertEqual(5, cmd.count_atoms('*bar*'))
        self.assertEqual(2, cmd.count_atoms('foo*bar'))
        self.assertEqual(3, cmd.count_atoms('*foo*bar*'))
        self.assertEqual(1, cmd.count_atoms('_*'))
        self.assertEqual(6, cmd.count_atoms('*'))
        self.assertEqual(6, cmd.count_atoms('*'))
        cmd.alter('*', 'name = "ABCD"')
        cmd.alter('foo', 'name = "ABED"')
        cmd.alter('bar', 'name = "ABEC"')
        self.assertEqual(6, cmd.count_atoms('name A*'))
        self.assertEqual(2, cmd.count_atoms('name ABE*'))
        self.assertEqual(1, cmd.count_atoms('name ABED*'))
        self.assertEqual(1, cmd.count_atoms('name AB*C'))
        self.assertEqual(5, cmd.count_atoms('name AB*C*'))
        self.assertEqual(5, cmd.count_atoms('name AB*D'))
        self.assertEqual(5, cmd.count_atoms('name AB*D*'))

    @testing.requires_version('2.2')
    def test_wildcard_sets_ranges(self):
        for i in range(10): cmd.pseudoatom('m%d' % i, chain=chr(65 + i))
        cmd.alter('m5', 'chain = "AB"')
        cmd.alter('m6', 'chain = "BA"')
        cmd.alter('m7', 'chain = "CC"')
        cmd.alter('m8', 'chain = "ZA"')
        cmd.alter('m9', 'chain = "ABA"')
        # A patterns
        self.assertEqual(1, cmd.count_atoms('chain A'))
        self.assertEqual(3, cmd.count_atoms('chain A*'))
        self.assertEqual(4, cmd.count_atoms('chain *A'))
        self.assertEqual(5, cmd.count_atoms('chain *A*'))
        self.assertEqual(1, cmd.count_atoms('chain A*A'))
        # B patterns
        self.assertEqual(2, cmd.count_atoms('chain B*'))
        self.assertEqual(2, cmd.count_atoms('chain *B'))
        self.assertEqual(4, cmd.count_atoms('chain *B*'))
        # X patterns (no matches)
        self.assertEqual(0, cmd.count_atoms('chain X*'))
        self.assertEqual(0, cmd.count_atoms('chain *X'))
        self.assertEqual(0, cmd.count_atoms('chain *X*'))
        # list with wildcards
        self.assertEqual(5, cmd.count_atoms('chain B*+A*'))
        self.assertEqual(3, cmd.count_atoms('chain B*+A*A'))
        self.assertEqual(3, cmd.count_atoms('chain B*+A*A+*X'))

        # lexicographical alpha ranges, A:C, will match AB (A <= AB <= C) but not CC (C < CC)
        # no wildcards in alpha ranges possible
        self.assertEqual(6, cmd.count_atoms('chain A:C'))
        self.assertEqual(6, cmd.count_atoms('chain A:CA'))
        self.assertEqual(7, cmd.count_atoms('chain A:CX')) # include CC
        self.assertEqual(6, cmd.count_atoms('chain A:C+Z'))
        self.assertEqual(7, cmd.count_atoms('chain A:C+Z*'))

    def test_sets_ranges(self):
        cmd.fab('ACDEFGHIKL')
        cmd.alter('resi 9', 'resi="9A"') # insertion code
        self.assertEqual(3, cmd.count_atoms('guide & resi 2-4'))
        self.assertEqual(3, cmd.count_atoms('guide & resi 2:4'))
        self.assertEqual(2, cmd.count_atoms('guide & resi 2+4'))
        self.assertEqual(4, cmd.count_atoms('guide & resi 2-4+6'))
        self.assertEqual(6, cmd.count_atoms('guide & resi 2-4+6-8'))
        self.assertEqual(0, cmd.count_atoms('guide & resi 9'))
        self.assertEqual(1, cmd.count_atoms('guide & resi 9A'))
        self.assertEqual(1, cmd.count_atoms('guide & resi 10'))
        self.assertEqual(0, cmd.count_atoms('guide & resi 10A'))
        self.assertEqual(2, cmd.count_atoms('guide & resi 9-10'))
        self.assertEqual(2, cmd.count_atoms('guide & resi 9A-10A'))
        self.assertEqual(10 + 9, cmd.count_atoms('name CA+CB'))
        self.assertEqual(10 + 9, cmd.count_atoms('name CA+CB+XYZ'))
        self.assertEqual(10, cmd.count_atoms('name C'))
        self.assertEqual(50, cmd.count_atoms('name C*'))
        self.assertEqual(10, cmd.count_atoms('name H'))
        self.assertEqual(24, cmd.count_atoms('name H*'))
        self.assertEqual(10, cmd.count_atoms('name *H'))
        self.assertEqual(74, cmd.count_atoms('name *H*'))
        self.assertEqual(20, cmd.count_atoms('name C+N'))
        self.assertEqual(23, cmd.count_atoms('name C+N*'))
        self.assertEqual(60, cmd.count_atoms('name C*+N'))
        self.assertEqual(63, cmd.count_atoms('name C*+N*'))
