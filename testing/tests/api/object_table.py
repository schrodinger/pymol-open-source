"""
Object table model: rows, visibility, counts, sorting, filtering, refresh
after delete. Headless (QtCore only).
"""

from pymol import cmd, testing


class TestObjectTable(testing.PyMOLTestCase):
    def test_model(self):
        from pymol.Qt import QtCore
        from pmg_qt.object_table import COLUMNS, ObjectTableModel, SortRole

        cmd.fragment("gly", "g1")
        cmd.fragment("ala", "a2")
        cmd.disable("a2")

        model = ObjectTableModel(cmd)
        self.assertTrue(model.refresh())
        self.assertEqual(model.rowCount(), 2)
        self.assertEqual(model.columnCount(), len(COLUMNS))
        rows = {r[0]: r for r in model.rows}
        self.assertEqual(rows["a2"][1:], ("molecule", False, cmd.count_atoms("a2"), 1))
        self.assertTrue(rows["g1"][2])

        proxy = QtCore.QSortFilterProxyModel()
        proxy.setSourceModel(model)
        proxy.setSortRole(SortRole)
        proxy.sort(3, QtCore.Qt.SortOrder.DescendingOrder)
        self.assertEqual(proxy.data(proxy.index(0, 0)), "a2")  # ala > gly atoms
        proxy.setFilterKeyColumn(0)
        proxy.setFilterFixedString("g")
        self.assertEqual(proxy.rowCount(), 1)

        cmd.delete("g1")
        self.assertTrue(model.refresh(full=False))
        self.assertEqual(model.rowCount(), 1)
        self.assertFalse(model.refresh())
