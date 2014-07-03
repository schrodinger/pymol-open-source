
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

