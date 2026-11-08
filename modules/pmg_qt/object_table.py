"""
Object table dock: one row per object with name, type, visibility, atom
count and state count. Sort by any column, filter by name, Zoom/Orient/Show/
Hide the selected rows, copy rows as tab-separated text.
"""

from pymol.Qt import QtCore, QtWidgets

Qt = QtCore.Qt

COLUMNS = ("Name", "Type", "Shown", "Atoms", "States")
SortRole = Qt.ItemDataRole.UserRole + 1


def snapshot(cmd, cache=None):
    """Rows for every object. `cache` maps name -> (atoms, states) to skip
    recounting atoms of objects that were already listed."""
    cache = cache or {}
    enabled = set(cmd.get_names("objects", enabled_only=1))
    rows = []
    for name in cmd.get_names("objects"):
        kind = cmd.get_type(name).replace("object:", "")
        if name in cache:
            atoms, states = cache[name]
        else:
            try:
                atoms = cmd.count_atoms("%" + name) if kind == "molecule" else None
                states = cmd.count_states("%" + name)
            except Exception:
                atoms, states = None, None
        rows.append((name, kind, name in enabled, atoms, states))
    return rows


class ObjectTableModel(QtCore.QAbstractTableModel):
    def __init__(self, cmd, parent=None):
        super().__init__(parent)
        self.cmd = cmd
        self.rows = []

    def refresh(self, full=True):
        """Reload. full=False keeps cached atom/state counts for known names
        (cheap, used by the poll). Returns True if the table changed."""
        cache = {} if full else {r[0]: (r[3], r[4]) for r in self.rows}
        rows = snapshot(self.cmd, cache)
        if rows == self.rows:
            return False
        self.beginResetModel()
        self.rows = rows
        self.endResetModel()
        return True

    def rowCount(self, parent=QtCore.QModelIndex()):
        return 0 if parent.isValid() else len(self.rows)

    def columnCount(self, parent=QtCore.QModelIndex()):
        return len(COLUMNS)

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):
        if (
            role == Qt.ItemDataRole.DisplayRole
            and orientation == Qt.Orientation.Horizontal
        ):
            return COLUMNS[section]
        return None

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        row = self.rows[index.row()]
        col = index.column()
        value = row[col]
        if role == Qt.ItemDataRole.DisplayRole:
            if col == 2:
                return None
            return "" if value is None else value
        if role == Qt.ItemDataRole.CheckStateRole and col == 2:
            return Qt.CheckState.Checked if value else Qt.CheckState.Unchecked
        if role == Qt.ItemDataRole.TextAlignmentRole and col in (3, 4):
            return int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        if role == Qt.ItemDataRole.UserRole:
            return row[0]
        if role == SortRole:
            if col in (3, 4):
                return -1 if value is None else value
            return int(value) if col == 2 else value
        return None


class ObjectTableDock(QtWidgets.QDockWidget):
    def __init__(self, window):
        super().__init__("Objects", window)
        self.setObjectName("dock_objects")
        self.cmd = window.cmd

        self.model = ObjectTableModel(self.cmd, self)
        self.proxy = QtCore.QSortFilterProxyModel(self)
        self.proxy.setSourceModel(self.model)
        self.proxy.setSortRole(SortRole)
        self.proxy.setFilterCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
        self.proxy.setFilterKeyColumn(0)

        self.view = QtWidgets.QTableView()
        self.view.setModel(self.proxy)
        self.view.setSortingEnabled(True)
        self.view.sortByColumn(0, Qt.SortOrder.AscendingOrder)
        self.view.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows
        )
        self.view.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self.view.verticalHeader().hide()
        self.view.horizontalHeader().setSectionResizeMode(
            0, QtWidgets.QHeaderView.ResizeMode.Stretch
        )
        self.view.setContextMenuPolicy(Qt.ContextMenuPolicy.ActionsContextMenu)
        self.view.doubleClicked.connect(lambda _: self.zoom_selected())

        filter_edit = QtWidgets.QLineEdit()
        filter_edit.setPlaceholderText("Filter by name")
        filter_edit.setClearButtonEnabled(True)
        filter_edit.textChanged.connect(self.proxy.setFilterFixedString)
        refresh = QtWidgets.QToolButton()
        refresh.setText("Refresh")
        refresh.setToolTip("Recount atoms and states")
        refresh.clicked.connect(self.refresh)

        top = QtWidgets.QHBoxLayout()
        top.addWidget(filter_edit)
        top.addWidget(refresh)
        layout = QtWidgets.QVBoxLayout()
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addLayout(top)
        layout.addWidget(self.view)
        body = QtWidgets.QWidget()
        body.setLayout(layout)
        self.setWidget(body)

        for text, shortcut, slot in (
            ("Zoom", "Space", self.zoom_selected),
            ("Orient", "Ctrl+Shift+O", self.orient_selected),
            ("Show", "Ctrl+Shift+S", self.show_selected),
            ("Hide", "Ctrl+Shift+H", self.hide_selected),
            ("Copy rows", "Ctrl+C", self.copy_selected),
        ):
            action = QtWidgets.QAction(text, self.view)
            action.setShortcut(shortcut)
            action.setShortcutContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
            action.triggered.connect(slot)
            self.view.addAction(action)

        # ponytail: 1 s poll of names/visibility (no atom data) while visible;
        # upgrade path is a change notification from the C++ executive.
        self._timer = QtCore.QTimer(self)
        self._timer.setInterval(1000)
        self._timer.timeout.connect(lambda: self._refresh(full=False))
        self.visibilityChanged.connect(
            lambda visible: self._timer.start() if visible else self._timer.stop()
        )

    def refresh(self):
        self._refresh(full=True)

    def _refresh(self, full):
        """Refresh the model and keep the row selection (by name) across the
        model reset."""
        names = set(self.selected_names())
        if self.model.refresh(full=full) and names:
            selection = self.view.selectionModel()
            selection.clear()
            for row in range(self.proxy.rowCount()):
                index = self.proxy.index(row, 0)
                if self.proxy.data(index, Qt.ItemDataRole.UserRole) in names:
                    selection.select(
                        index,
                        QtCore.QItemSelectionModel.SelectionFlag.Select
                        | QtCore.QItemSelectionModel.SelectionFlag.Rows,
                    )

    def selected_names(self):
        return [
            self.proxy.data(self.proxy.index(i.row(), 0), Qt.ItemDataRole.UserRole)
            for i in self.view.selectionModel().selectedRows()
        ]

    def _selection(self):
        return " or ".join("%" + n for n in self.selected_names())

    def zoom_selected(self):
        if self._selection():
            self.cmd.zoom(self._selection(), animate=1.0)

    def orient_selected(self):
        if self._selection():
            self.cmd.orient(self._selection(), animate=1.0)

    def show_selected(self):
        for name in self.selected_names():
            self.cmd.enable(name)
        self._refresh(full=False)

    def hide_selected(self):
        for name in self.selected_names():
            self.cmd.disable(name)
        self._refresh(full=False)

    def copy_selected(self):
        names = set(self.selected_names())
        lines = ["\t".join(COLUMNS)] + [
            "\t".join("" if v is None else str(v) for v in row)
            for row in self.model.rows
            if row[0] in names
        ]
        QtWidgets.QApplication.clipboard().setText("\n".join(lines))
