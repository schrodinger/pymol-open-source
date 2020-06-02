
from pymol import cmd, testing, stored

class TestSelecting(testing.PyMOLTestCase):

    @testing.requires_version('2.4')
    def test_error_propagation(self):
        try:
            cmd.center("(all) and no_such_name")
        except Exception as e:
            self.assertTrue("Invalid selection name" in e.message)
            self.assertTrue("no_such_name" in e.message)
            return
        self.fail("did not raise")

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
        cmd.fragment("gly", "m1")
        NC = 2

        # auto_number_selections=0
        cmd.select("elem C")
        self.assertEquals(NC, cmd.count_atoms("sele"))
        self.assertEquals(0, cmd.get_setting_int("sel_counter"))

        # auto_number_selections=1
        cmd.set('auto_number_selections', 1)
        cmd.set('sel_counter', 3)
        cmd.select("elem C")
        self.assertEquals(NC, cmd.count_atoms("sel04"))
        self.assertEquals(4, cmd.get_setting_int("sel_counter"))

        cmd.set('auto_number_selections', 1)
        cmd.select(None, "elem C")
        self.assertEquals(NC, cmd.count_atoms("sel05"))
        self.assertEquals(5, cmd.get_setting_int("sel_counter"))

        # name=None always numbers the selection
        cmd.set('auto_number_selections', 0)
        cmd.select(None, "elem C")
        self.assertEquals(NC, cmd.count_atoms("sel06"))
        self.assertEquals(6, cmd.get_setting_int("sel_counter"))

        # default
        cmd.select("foo", "elem C")
        self.assertEquals(NC, cmd.count_atoms("foo"))
        self.assertEquals(["foo"], cmd.get_names("selections", enabled_only=1))

        # merge with non-existing, enable=0
        cmd.delete('foo')
        cmd.select("foo", "elem C", 0, merge=1)
        self.assertEquals(NC, cmd.count_atoms("foo"))
        self.assertEquals([], cmd.get_names("selections", enabled_only=1))

        # merge, enable=1
        cmd.select("foo", "elem N", 1)
        self.assertEquals(1, cmd.count_atoms("foo"))
        self.assertEquals(["foo"], cmd.get_names("selections", enabled_only=1))

        cmd.select("foo", "elem C", -1, merge=1)
        self.assertEquals(NC + 1, cmd.count_atoms("foo"))
        self.assertEquals(["foo"], cmd.get_names("selections", enabled_only=1))

        cmd.select("foo", "elem O", -1, merge=2)
        self.assertEquals(NC + 2, cmd.count_atoms("foo"))
        self.assertEquals(["foo"], cmd.get_names("selections", enabled_only=1))

        # merge, enable=0
        cmd.select("foo", "elem N", 0)
        self.assertEquals(1, cmd.count_atoms("foo"))
        self.assertEquals([], cmd.get_names("selections", enabled_only=1))

        cmd.select("foo", "elem C", -1, merge=1)
        self.assertEquals(NC + 1, cmd.count_atoms("foo"))
        self.assertEquals([], cmd.get_names("selections", enabled_only=1))

        cmd.select("foo", "elem O", -1, merge=2)
        self.assertEquals(1, cmd.count_atoms("foo"))
        self.assertEquals([], cmd.get_names("selections", enabled_only=1))

        # state
        cmd.create('m1', 'm1 & elem C', 1, 2)
        self.assertEquals(7, cmd.select('present'))
        self.assertEquals(7, cmd.select('present', state=1))
        self.assertEquals(2, cmd.select('present', state=2))
        self.assertEquals(0, cmd.select('present', state=3))

        # domain
        cmd.select("foo", "elem C")
        cmd.select("bar", "name CA+N+O", domain="foo")
        self.assertEquals(1, cmd.count_atoms("bar"))

    @testing.requires("multi_undo")
    def testUndoSelect(self):
        cmd.fragment("gly", "m1")
        NC = 2

        # auto_number_selections=0
        cmd.select("elem C")
        self.assertEquals(NC, cmd.count_atoms("sele"))
        self.assertEquals(["sele"], cmd.get_names("selections", enabled_only=1))
        cmd.undo2()
        self.assertEquals([], cmd.get_names("selections", enabled_only=1))
        cmd.redo2()
        self.assertEquals(NC, cmd.count_atoms("sele"))
        self.assertEquals(["sele"], cmd.get_names("selections", enabled_only=1))

        self.assertEquals(0, cmd.get_setting_int("sel_counter"))

        # auto_number_selections=1
        cmd.set('auto_number_selections', 1)
        cmd.set('sel_counter', 3)
        cmd.select("elem C")
        self.assertEquals(NC, cmd.count_atoms("sel04"))
        cmd.undo2()
        self.assertTrue("sel04" not in cmd.get_names("selections", enabled_only=1))
        cmd.redo2()
        self.assertEquals(NC, cmd.count_atoms("sel04"))
        self.assertEquals(4, cmd.get_setting_int("sel_counter"))

        cmd.set('auto_number_selections', 1)
        cmd.select(None, "elem C")
        self.assertEquals(NC, cmd.count_atoms("sel05"))
        cmd.undo2()
        self.assertTrue("sel05" not in cmd.get_names("selections", enabled_only=1))
        cmd.redo2()
        self.assertEquals(NC, cmd.count_atoms("sel05"))
        self.assertEquals(5, cmd.get_setting_int("sel_counter"))

        # name=None always numbers the selection
        cmd.set('auto_number_selections', 0)
        cmd.select(None, "elem C")
        self.assertEquals(NC, cmd.count_atoms("sel06"))
        cmd.undo2()
        self.assertTrue("sel06" not in cmd.get_names("selections", enabled_only=1))
        cmd.redo2()
        self.assertEquals(NC, cmd.count_atoms("sel06"))
        self.assertEquals(6, cmd.get_setting_int("sel_counter"))

        # default
        cmd.select("foo", "elem C")
        self.assertEquals(NC, cmd.count_atoms("foo"))
        self.assertEquals(["foo"], cmd.get_names("selections", enabled_only=1))
        cmd.undo2()
        self.assertTrue("foo" not in cmd.get_names("selections", enabled_only=1))
        cmd.redo2()
        self.assertEquals(NC, cmd.count_atoms("foo"))
        self.assertEquals(["foo"], cmd.get_names("selections", enabled_only=1))

        # merge with non-existing, enable=0
        cmd.delete("foo")
        cmd.select("foo", "elem C", 0, merge=1)
        self.assertEquals(NC, cmd.count_atoms("foo"))
        self.assertEquals([], cmd.get_names("selections", enabled_only=1))
        self.assertTrue("foo" in cmd.get_names("selections", enabled_only=0))
        cmd.undo2()
        self.assertTrue("foo" not in cmd.get_names("selections", enabled_only=0))
        cmd.redo2()
        self.assertEquals(NC, cmd.count_atoms("foo"))
        self.assertEquals([], cmd.get_names("selections", enabled_only=1))
        self.assertTrue("foo" in cmd.get_names("selections", enabled_only=0))

        # merge, enable=1
        cmd.select("foo", "elem N", 1)
        self.assertEquals(1, cmd.count_atoms("foo"))
        self.assertEquals(["foo"], cmd.get_names("selections", enabled_only=1))
        cmd.undo2()
        self.assertEquals([], cmd.get_names("selections", enabled_only=1))
        cmd.redo2()
        self.assertEquals(1, cmd.count_atoms("foo"))
        self.assertEquals(["foo"], cmd.get_names("selections", enabled_only=1))

        cmd.select("foo", "elem C", -1, merge=1)
        self.assertEquals(NC + 1, cmd.count_atoms("foo"))
        self.assertEquals(["foo"], cmd.get_names("selections", enabled_only=1))
        cmd.undo2()
        self.assertEquals(NC - 1, cmd.count_atoms("foo"))
        self.assertEquals(["foo"], cmd.get_names("selections", enabled_only=1))
        cmd.redo2()
        self.assertEquals(NC + 1, cmd.count_atoms("foo"))
        self.assertEquals(["foo"], cmd.get_names("selections", enabled_only=1))

        cmd.select("foo", "elem O", -1, merge=2)
        self.assertEquals(NC + 2, cmd.count_atoms("foo"))
        self.assertEquals(["foo"], cmd.get_names("selections", enabled_only=1))

        # merge, enable=0
        cmd.select("foo", "elem N", 1)
        cmd.select("foo", "elem N", 0)
        self.assertEquals(1, cmd.count_atoms("foo"))
        self.assertEquals([], cmd.get_names("selections", enabled_only=1))
        cmd.undo2()
        self.assertEquals(1, cmd.count_atoms("foo"))
        self.assertEquals(["foo"], cmd.get_names("selections", enabled_only=1))
        cmd.redo2()
        self.assertEquals(1, cmd.count_atoms("foo"))
        self.assertEquals([], cmd.get_names("selections", enabled_only=1))

        cmd.select("foo", "elem C", -1, merge=1)
        self.assertEquals(NC + 1, cmd.count_atoms("foo"))
        self.assertEquals([], cmd.get_names("selections", enabled_only=1))
        cmd.undo2()
        self.assertEquals(NC - 1, cmd.count_atoms("foo"))
        self.assertEquals([], cmd.get_names("selections", enabled_only=1))
        cmd.redo2()
        self.assertEquals(NC + 1, cmd.count_atoms("foo"))
        self.assertEquals([], cmd.get_names("selections", enabled_only=1))

        cmd.select("foo", "elem O", -1, merge=2)
        self.assertEquals(NC - 1, cmd.count_atoms("foo"))
        self.assertEquals([], cmd.get_names("selections", enabled_only=1))
        cmd.undo2()
        self.assertEquals(NC + 1, cmd.count_atoms("foo"))
        self.assertEquals([], cmd.get_names("selections", enabled_only=1))
        cmd.redo2()
        self.assertEquals(NC - 1, cmd.count_atoms("foo"))
        self.assertEquals([], cmd.get_names("selections", enabled_only=1))

        # state
        cmd.delete("sele")
        cmd.create('m1', 'm1 & elem C', 1, 2)
        cmd.select('present')
        self.assertEquals(7, cmd.count_atoms('present'))
        self.assertEquals(["sele"], cmd.get_names("selections", enabled_only=1))
        cmd.undo2()
        self.assertEquals(7, cmd.count_atoms('present'))
        self.assertEquals([], cmd.get_names("selections", enabled_only=1))
        cmd.redo2()
        self.assertEquals(7, cmd.count_atoms('present'))
        self.assertEquals(["sele"], cmd.get_names("selections", enabled_only=1))


        cmd.delete("sele")
        cmd.select('present', state=1)
        self.assertEquals(7, cmd.count_atoms('present', state=1))
        self.assertEquals(["sele"], cmd.get_names("selections", enabled_only=1))
        cmd.undo2()
        self.assertEquals(7, cmd.count_atoms('present', state=1))
        self.assertEquals([], cmd.get_names("selections", enabled_only=1))
        cmd.redo2()
        self.assertEquals(7, cmd.count_atoms('present', state=1))
        self.assertEquals(["sele"], cmd.get_names("selections", enabled_only=1))


        cmd.delete("sele")
        cmd.select('present', state=2)
        self.assertEquals(2, cmd.count_atoms('present', state=2))
        self.assertEquals(["sele"], cmd.get_names("selections", enabled_only=1))
        cmd.undo2()
        self.assertEquals(2, cmd.count_atoms('present', state=2))
        self.assertEquals([], cmd.get_names("selections", enabled_only=1))
        cmd.redo2()
        self.assertEquals(2, cmd.count_atoms('present', state=2))
        self.assertEquals(["sele"], cmd.get_names("selections", enabled_only=1))

        cmd.delete("sele")
        cmd.select('present', state=3)
        self.assertEquals(0, cmd.count_atoms('present', state=3))
        self.assertEquals(["sele"], cmd.get_names("selections", enabled_only=1))
        cmd.undo2()
        self.assertEquals(0, cmd.count_atoms('present', state=3))
        self.assertEquals([], cmd.get_names("selections", enabled_only=1))
        cmd.redo2()
        self.assertEquals(0, cmd.count_atoms('present', state=3))
        self.assertEquals(["sele"], cmd.get_names("selections", enabled_only=1))

        # domain
        cmd.delete("foo")
        cmd.select("foo", "elem C")
        cmd.select("bar", "name CA+N+O", domain="foo") #foo is disabled; must test enabled state
        self.assertEquals(1, cmd.count_atoms("bar"))
        self.assertTrue("foo" not in cmd.get_names("selections", enabled_only=1))
        self.assertTrue("bar" in cmd.get_names("selections", enabled_only=1))
        cmd.undo2()
        self.assertTrue("foo" in cmd.get_names("selections", enabled_only=1))
        self.assertTrue("bar" not in cmd.get_names("selections", enabled_only=1))
        cmd.redo2()
        self.assertEquals(1, cmd.count_atoms("bar"))
        self.assertTrue("foo" not in cmd.get_names("selections", enabled_only=1))
        self.assertTrue("bar" in cmd.get_names("selections", enabled_only=1))

    def testSelectList(self):
        cmd.fragment('ala', 'm1')
        cmd.alter('m1', 'rank += 10')
        cmd.alter('m1', 'ID += 20')
        cmd.select_list('s1', 'm1', [1, 4], mode='index')
        self.assertEqual(cmd.count_atoms('s1'), 2)
        self.assertEqual(cmd.count_atoms('s1 & elem N'), 1)
        self.assertEqual(cmd.count_atoms('s1 & elem O'), 1)
        cmd.select_list('s1', 'm1', [15, 16, 17, 18, 19], mode='rank')
        self.assertEqual(cmd.count_atoms('s1'), 5)
        self.assertEqual(cmd.count_atoms('s1 & elem H'), 5)
        cmd.select_list('s1', 'm1', [21, 22, 23], mode='id')
        self.assertEqual(cmd.count_atoms('s1'), 3)
        self.assertEqual(cmd.count_atoms('s1 & elem C'), 3)

    @testing.requires_version('2.4')
    def testSelectList_state(self):
        cmd.fragment('ala', 'm1')
        cmd.create('m1', 'm1 & elem C', 1, 2)
        cmd.select_list('s1', 'm1', [100] * 100 + [2, 3, 4], state=1, mode='rank')
        self.assertEqual(cmd.count_atoms('s1'), 3)
        cmd.select_list('s1', 'm1', [100] * 100 + [2, 3, 4], state=2, mode='rank')
        self.assertEqual(cmd.count_atoms('s1'), 2)

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

    @testing.foreach('origin', 'center')
    def test_dummy_selectors(self, op):
        self.assertEqual(cmd.count_atoms('all'), 0)
        self.assertEqual(cmd.count_atoms(op), 0)
        self.assertEqual(cmd.count_atoms(op + ' around 2'), 0)
        cmd.fragment('arg')
        self.assertEqual(cmd.count_atoms('all'), 24)
        self.assertEqual(cmd.count_atoms(op), 0)
        self.assertEqual(cmd.count_atoms(op + ' around 2'), 8)
        self.assertEqual(cmd.count_atoms(op + ' expand 2'), 8)
        self.assertEqual(cmd.count_atoms(op + ' within 2 of arg'), 0)
        self.assertEqual(cmd.count_atoms('arg within 2 of ' + op), 8)
        self.assertEqual(cmd.count_atoms('arg near_to 2 of ' + op), 8)
        self.assertEqual(cmd.count_atoms('arg beyond 3 of ' + op), 13)
        self.assertEqual(cmd.count_atoms(op + ' gap 3'), 5)

    def test_offcenter_origin_selector(self):
        cmd.fragment('arg')
        cmd.origin('elem O')
        self.assertEqual(cmd.count_atoms('origin around 2'), 2)
        self.assertEqual(cmd.count_atoms('origin expand 2'), 2)
        self.assertEqual(cmd.count_atoms('arg within 2 of origin'), 2)
        self.assertEqual(cmd.count_atoms('arg near_to 2 of origin'), 2)
        self.assertEqual(cmd.count_atoms('arg beyond 3 of origin'), 19)
        self.assertEqual(cmd.count_atoms('origin gap 3'), 14)
        # TODO Also test 'center' selector here. I'd expect that it gives the
        # same result as the test_dummy_selectors test, but numbers are
        # different!
