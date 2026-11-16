"""
Group commands of modules/pymol/grouping.py: stepping, enable and disable,
copy, sort, chains, and their rows in the group action menu.
"""

from pymol import CmdException, cmd, grouping, menu, testing


class TestGrouping(testing.PyMOLTestCase):
    def setUp(self):
        super().setUp()
        for name in ("m2", "M1", "m3"):
            cmd.fragment("ala", name)
        cmd.group("g", "m2 M1 m3")

    def shown(self):
        enabled = cmd.get_names("objects", enabled_only=1)
        return [n for n in grouping.members("g") if n in enabled]

    def test_step(self):
        cmd.group_disable("g")
        self.assertEqual(self.shown(), [])
        self.assertEqual(cmd.group_next(), "m2")  # the only group is found
        self.assertEqual(cmd.group_next("g"), "M1")
        self.assertEqual(self.shown(), ["M1"])
        cmd.group_next("g")
        self.assertEqual(cmd.group_next("g"), "m2")  # wraps around
        self.assertEqual(cmd.group_previous("g"), "m3")
        cmd.group_disable("g")
        self.assertEqual(cmd.group_previous("g"), "m3")  # starts at the end
        cmd.disable("g")
        cmd.group_next("g")  # shows the group as well
        self.assertIn("g", cmd.get_names("objects", enabled_only=1))
        # two groups with visible members: the group has to be named
        cmd.fragment("gly", "other")
        cmd.group("g2", "other")
        self.assertRaises(CmdException, cmd.group_next)
        cmd.disable("other")
        self.assertEqual(cmd.group_next(), "M1")  # g is the one in use
        self.assertRaises(CmdException, cmd.group_next, "m2")

    def test_enable(self):
        self.assertEqual(cmd.group_disable("g"), ["m2", "M1", "m3"])
        self.assertEqual(cmd.group_enable("g", 2), ["m2", "m3"])
        self.assertEqual(self.shown(), ["m2", "m3"])
        self.assertEqual(cmd.group_enable("g", 0), ["m2", "M1", "m3"])

    def test_copy(self):
        cmd.pseudoatom("p")
        cmd.group("inner", "p")
        cmd.group("g", "inner")
        cmd.disable("M1")
        made = cmd.copy_group("g")
        self.assertEqual(made, ["m2_cpy", "M1_cpy", "m3_cpy", "inner_cpy", "p_cpy"])
        self.assertEqual(
            grouping.members("g_cpy"), ["m2_cpy", "M1_cpy", "m3_cpy", "inner_cpy"]
        )
        self.assertEqual(grouping.members("inner_cpy"), ["p_cpy"])
        self.assertEqual(cmd.count_atoms("m2_cpy"), cmd.count_atoms("m2"))
        enabled = cmd.get_names("objects", enabled_only=1)
        self.assertNotIn("M1_cpy", enabled)
        self.assertIn("m2_cpy", enabled)
        # a second copy would overwrite: nothing is made
        before = cmd.get_names()
        self.assertRaises(CmdException, cmd.copy_group, "g")
        self.assertEqual(cmd.get_names(), before)
        cmd.copy_group("g", "backup", "_bak")
        self.assertEqual(len(grouping.members("backup")), 4)

    def test_sort(self):
        self.assertFalse(grouping.is_sorted("g"))
        self.assertEqual(cmd.sort_group("g"), ["M1", "m2", "m3"])
        self.assertEqual(grouping.members("g"), ["M1", "m2", "m3"])
        self.assertTrue(grouping.is_sorted("g"))
        cmd.sort_group("g", 1)
        self.assertEqual(grouping.members("g"), ["m3", "m2", "M1"])

    def test_extract_chain(self):
        cmd.alter("m2 and name CA+CB", "chain='A'")
        cmd.alter("m2 and not chain A", "chain='B'")
        cmd.alter("m3", "chain='A'")
        self.assertEqual(cmd.extract_chain("m2"), ["m2_A", "m2_B"])
        self.assertEqual(cmd.count_atoms("m2_A"), 2)
        self.assertEqual(cmd.count_atoms("m2"), 10)  # source unchanged
        cmd.delete("m2_*")
        made = cmd.extract_chain("g", "A", "chains_A")
        self.assertEqual(made, ["m2_A", "m3_A"])
        self.assertEqual(grouping.members("chains_A"), made)
        self.assertRaises(CmdException, cmd.extract_chain, "g", "A")  # exists
        self.assertRaises(CmdException, cmd.extract_chain, "g", "Z")
        self.assertEqual(cmd.extract_chain("M1"), ["M1_"])  # blank chain

    def test_menu(self):
        rows = {row[1]: row[2] for row in menu.group_action(cmd, "g") if row[0] == 1}
        cmd.do(rows["sort members"], echo=0)
        self.assertEqual(grouping.members("g"), ["M1", "m2", "m3"])
        rows = {row[1]: row[2] for row in menu.group_action(cmd, "g") if row[0] == 1}
        cmd.do(rows["sort members"], echo=0)  # sorted already: the row reverses
        self.assertEqual(grouping.members("g"), ["m3", "m2", "M1"])
        cmd.do(rows["disable members"], echo=0)
        self.assertEqual(self.shown(), [])
        cmd.do(rows["enable members"], echo=0)
        self.assertEqual(len(self.shown()), 3)
        cmd.do(rows["duplicate group"], echo=0)
        self.assertEqual(len(grouping.members("g_cpy")), 3)

    def test_registered(self):
        from pymol.shortcut_dict import shortcut_dict_ref

        for name in (
            "group_next",
            "group_previous",
            "group_enable",
            "group_disable",
            "copy_group",
            "sort_group",
            "extract_chain",
        ):
            self.assertIn(name, cmd.keyword)
            self.assertIn("DESCRIPTION", getattr(cmd, name).__doc__)
        self.assertEqual(shortcut_dict_ref["CTRL-N"][0], "group_next")
        self.assertEqual(shortcut_dict_ref["CTRL-P"][0], "group_previous")
