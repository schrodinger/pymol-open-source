
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
