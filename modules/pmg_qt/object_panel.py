"""
Objects panel: the Qt replacement for the object list and the wizard buttons
of the viewport's internal GUI. One row per object or selection (group members
nested) with a visibility dot and the A/S/H/L/C menus from pymol.menu; a click
on the row shows or hides the object. Under the list: the states panel (step,
play, speed) while there is more than one state or frame.
Everything is Qt text, so it follows the application font, the theme and
the toolbar's UI scale, and the dock is resized by dragging its edge.
"""

import re

from pymol import CmdException, _cmd
from pymol import menu as menus
from pymol.Qt import QtCore, QtGui, QtWidgets

from .icons import menu_icon

Qt = QtCore.Qt
Checked, Unchecked = Qt.CheckState.Checked, Qt.CheckState.Unchecked
UserRole = Qt.ItemDataRole.UserRole

BUTTONS = (
    ("A", "Action"),
    ("S", "Show"),
    ("H", "Hide"),
    ("L", "Label"),
    ("C", "Color"),
)

# pymol.menu function per row type and button; same table as the C++ panel
# (layer3/Executive.cpp, CExecutive::click). None: no such button.
_MOL = ("mol_show", "mol_hide", "mol_labels", "mol_color")
_SIMPLE = ("simple_action", "cgo_show", "cgo_hide", None, "general_color")
MENUS = {
    "all": ("all_action",) + _MOL,
    "selection": ("sele_action",) + _MOL,
    "object:molecule": ("mol_action",) + _MOL,
    "object:group": ("group_action",) + _MOL,
    "object:map": ("map_action", "map_show", "map_hide", None, "general_color"),
    "object:mesh": ("mesh_action", "mesh_show", "mesh_hide", None, "mesh_color"),
    "object:surface": (
        "surface_action",
        "surface_show",
        "surface_hide",
        None,
        "mesh_color",
    ),
    "object:slice": (
        "slice_action",
        "slice_show",
        "slice_hide",
        None,
        "slice_color",
    ),
    "object:measurement": (
        "simple_action",
        "measurement_show",
        "measurement_hide",
        None,
        "measurement_color",
    ),
    "object:cgo": _SIMPLE,
    "object:alignment": _SIMPLE,
    "object:volume": (
        "simple_action",
        "volume_show",
        "volume_hide",
        None,
        "vol_color",
    ),
    "object:ramp": ("ramp_action", None, None, None, "ramp_color"),
    "object:": ("simple_action", None, None, None, None),  # callback objects
}


def snapshot(cmd):
    """Rows (name, kind, enabled, parent group) in panel order, "all" first."""
    enabled = set(cmd.get_names("public", enabled_only=1))
    parent = {}
    for group in cmd.get_names("public_group_objects"):
        with cmd.lockcm:
            members = _cmd.get_group_members(cmd._COb, group)
        for member in members.split():
            parent[member] = group
    rows = [("all", "all", bool(enabled), "")]
    for name in cmd.get_names("public"):
        rows.append((name, cmd.get_type(name), name in enabled, parent.get(name, "")))
    return rows


def menu_entries(cmd, name, kind, button):
    """The pymol.menu list behind one A/S/H/L/C button, or None."""
    function = MENUS.get(kind, MENUS["object:"])[button]
    if not function:
        return None
    if kind == "object:surface" and button == 4:
        return menus.mesh_color(cmd, name, "surface")
    return getattr(menus, function)(cmd, name)


def describe(cmd, name, kind):
    """(type, atoms, states) for tooltips and copied rows; None if unknown."""
    target = name if kind == "all" else "%" + name
    try:
        atoms = cmd.count_atoms(target) if kind != "object:group" else None
        states = cmd.count_states(target)
    except CmdException:  # maps, CGOs and friends have no atoms
        atoms = states = None
    return kind.replace("object:", "") or "object", atoms, states


SELECTING = (
    "Atoms",
    "Residues",
    "Chains",
    "Segments",
    "Objects",
    "Molecules",
    "C-alphas",
)


_COLOR_CODE = re.compile(r"\\(\d{3}|---)")


def plain_label(text):
    """Menu labels carry \\RGB color markup (digits 0-9 per channel). Returns
    the text without it and the first color, or None."""
    match = re.search(r"\\(\d)(\d)(\d)", text)
    color = (
        QtGui.QColor(*(int(c) * 255 // 9 for c in match.groups())) if match else None
    )
    return _COLOR_CODE.sub("", text).strip().replace("&", "&&"), color


def _swatch(color):
    pixmap = QtGui.QPixmap(12, 12)
    pixmap.fill(color)
    return QtGui.QIcon(pixmap)


def fill_menu(menu, entries, cmd):
    """Fill a QMenu from a pymol.menu list of [kind, label, command]:
    kind 0 separator, 1 item, 2 title; command is a command string, a
    submenu list, or a callable returning one (filled when first opened)."""

    def fill_once(sub, command):
        if sub.isEmpty():
            fill_menu(sub, command() if callable(command) else command, cmd)

    for entry in entries:
        kind = entry[0]
        if kind == 0:
            menu.addSeparator()
            continue
        text, color = plain_label(entry[1])
        command = entry[2] if len(entry) > 2 else ""
        if kind == 2:  # a caption: not selectable, drawn muted
            title = QtWidgets.QWidgetAction(menu)
            title.setDefaultWidget(QtWidgets.QLabel(text, objectName="menu_title"))
            title.setEnabled(False)
            menu.addAction(title)
            continue
        if isinstance(command, str):
            action = menu.addAction(text)
            action.triggered.connect(lambda _=False, c=command: cmd.do(c, echo=0))
        else:
            sub = menu.addMenu(text)
            sub.aboutToShow.connect(lambda s=sub, c=command: fill_once(s, c))
            action = sub.menuAction()
        icon = menu_icon(text, color or menu.palette().text().color())
        if icon or color:
            action.setIcon(icon or _swatch(color))


def context_menu(parent, cmd, name, kind):
    """The right-click menu of a row: show as and color, then the rename,
    duplicate and delete rows of its action menu, then its other menus."""
    popup = QtWidgets.QMenu(parent)
    buttons = MENUS.get(kind, MENUS["object:"])

    def entries(button):
        return menu_entries(cmd, name, kind, button) or []

    def later(title, rows):  # a submenu, filled when it is first opened
        sub = popup.addMenu(title)
        sub.aboutToShow.connect(lambda: sub.isEmpty() and fill_menu(sub, rows(), cmd))

    show_as = [row[2] for row in entries(1) if row[1] == "as"]
    if show_as:
        later("Show as", lambda: show_as[0])
    if buttons[4]:
        later("Color", lambda: entries(4))
    popup.addSeparator()
    first_words = ("rename", "duplicate", "delete")
    flat = [
        row
        for row in entries(0)
        if row[0] == 1 and plain_label(row[1])[0].split(" ")[0] in first_words
    ]
    fill_menu(popup, flat, cmd)
    popup.addSeparator()
    for button, (_, title) in enumerate(BUTTONS[:4]):
        if buttons[button]:
            later(title, lambda b=button: entries(b))
    return popup


def popup_at_cursor(parent, entries, cmd):
    """Show a pymol.menu list as a context menu at the mouse pointer, once
    the caller has returned to the event loop. The core calls this through
    pymol.menu.popup in the middle of handling a click."""

    def show():
        menu = QtWidgets.QMenu(parent)
        menu.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        fill_menu(menu, entries, cmd)
        menu.popup(QtGui.QCursor.pos())

    QtCore.QTimer.singleShot(0, show)
    return True


# chip colors A S H L C: three blues, violet, green
ACCENTS = ("#1e4fb8", "#2a73d9", "#3f93e6", "#8e5bd6", "#4caf50")


class ChipDelegate(QtWidgets.QStyledItemDelegate):
    """Paints the A/S/H/L/C letters as rounded chips. The chips of a hidden
    object are dimmed until the mouse is over its row."""

    def paint(self, painter, option, index):
        letter = index.data()
        if not letter:
            return
        first = index.siblingAtColumn(0)
        # the role comes back as a plain int, not as the enum
        hidden = (
            Qt.CheckState(first.data(Qt.ItemDataRole.CheckStateRole)) != Checked
            and first.data(UserRole)[1] != "all"
        )
        painter.save()
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        if hidden and not option.state & QtWidgets.QStyle.StateFlag.State_MouseOver:
            painter.setOpacity(0.4)
        side = min(option.rect.width(), option.rect.height()) - 4
        chip = QtCore.QRect(0, 0, side, side)
        chip.moveCenter(option.rect.center())
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QtGui.QColor(ACCENTS[index.column() - 1]))
        painter.drawRoundedRect(chip, 4, 4)
        font = QtGui.QFont(option.font)
        font.setBold(True)
        painter.setFont(font)
        painter.setPen(Qt.GlobalColor.white)
        painter.drawText(chip, Qt.AlignmentFlag.AlignCenter, letter)
        painter.restore()

    def sizeHint(self, option, index):
        side = option.fontMetrics.height() + 10
        return QtCore.QSize(side, side)


class WizardPanel(QtWidgets.QGroupBox):
    """Title, buttons and menus of the active wizard (wizard.get_panel()).
    Hidden while no wizard runs."""

    changed = QtCore.Signal()

    def __init__(self, cmd):
        super().__init__("Wizard")
        self.cmd = cmd
        self.panel = None
        self.buttons = QtWidgets.QVBoxLayout(self)
        self.hide()

    def refresh(self):
        wizard = self.cmd.get_wizard()
        panel = wizard.get_panel() if wizard else None
        if panel == self.panel:
            return
        self.panel = panel
        while self.buttons.count():
            self.buttons.takeAt(0).widget().deleteLater()
        for kind, label, code in panel or ():
            text = plain_label(label)[0]
            if kind == 1:  # title
                self.setTitle(text)
                continue
            button = QtWidgets.QPushButton(text)
            if kind == 3:  # popup: code names a menu of the wizard
                popup = QtWidgets.QMenu(button)
                popup.aboutToShow.connect(lambda m=popup, c=code: self._fill(m, c))
                popup.triggered.connect(lambda _: self.changed.emit())
                button.setMenu(popup)
            else:
                button.clicked.connect(lambda _=False, c=code: self._run(c))
            self.buttons.addWidget(button)
        self.setVisible(bool(panel))

    def _fill(self, popup, code):
        popup.clear()
        fill_menu(popup, self.cmd.get_wizard().get_menu(code) or (), self.cmd)

    def _run(self, code):
        self.cmd.do(code, echo=0)
        self.changed.emit()


class StatesPanel(QtWidgets.QWidget):
    """First, previous, play/stop, next, last, a counter to jump to a number
    and the playback speed (movie_fps). Without a movie a frame is a state,
    so the frame commands drive both. Hidden with a single state."""

    BUTTONS = (  # icon (icons.LABELS), tooltip, command
        ("first", "First", "rewind"),
        ("previous", "Previous", "backward"),
        ("play", "Play", "mtoggle"),
        ("forward", "Next", "forward"),
        ("last", "Last", "ending"),
    )
    SPEEDS = (1, 5, 10, 15, 30, 60)

    def __init__(self, cmd):
        super().__init__()
        self.cmd = cmd
        self.playing = False
        self.title = QtWidgets.QLabel("States")
        self.speed = QtWidgets.QToolButton()
        self.speed.setToolTip("Playback speed (movie_fps)")
        self.speed.setPopupMode(QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup)
        self.speed.setMenu(QtWidgets.QMenu(self.speed))
        for fps in self.SPEEDS:
            self.speed.menu().addAction(
                "%d fps" % fps, lambda f=fps: self._run("set", "movie_fps", f)
            )
        self.counter = QtWidgets.QSpinBox()
        self.counter.setMinimum(1)
        self.counter.setToolTip("Type a number to jump there")
        self.counter.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.counter.setKeyboardTracking(False)  # jump on Enter, not per digit
        self.counter.valueChanged.connect(lambda n: self._run("frame", n))

        top = QtWidgets.QHBoxLayout()
        top.addWidget(self.title, 1)
        top.addWidget(self.speed)
        row = QtWidgets.QHBoxLayout()
        row.setSpacing(2)
        row.addWidget(self.counter, 1)
        self.buttons = {}
        for icon, tip, command in self.BUTTONS:
            button = self.buttons[icon] = QtWidgets.QToolButton()
            button.setToolTip(tip)
            button.setAutoRaise(True)
            button.setAutoRepeat(command in ("backward", "forward"))
            button.clicked.connect(lambda _=False, c=command: self._run(c))
            row.addWidget(button)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 4, 0, 4)
        layout.addLayout(top)
        layout.addLayout(row)
        self._set_icons()

        # the window polls twice a second: too coarse for a running counter
        self.timer = QtCore.QTimer(self, interval=100)
        self.timer.timeout.connect(self.refresh)
        self.hide()

    def _run(self, command, *args):
        getattr(self.cmd, command)(*args)
        self.refresh()

    def _set_icons(self):
        color = self.palette().buttonText().color()
        for label, button in self.buttons.items():
            if label == "play" and self.playing:
                label = "stop"
            button.setIcon(menu_icon(label, color))

    def changeEvent(self, event):
        if event.type() == QtCore.QEvent.Type.PaletteChange:
            self._set_icons()  # they are tinted: follow Display > Dark Theme
        super().changeEvent(event)

    # ponytail: polled like the object list; a state change made by command
    # shows up to 0.5 s late. Upgrade path: a notification from the C++ scene.
    def refresh(self):
        cmd = self.cmd
        count = cmd.count_frames()
        self.setVisible(count > 1)
        playing = count > 1 and bool(cmd.get_movie_playing())
        if playing != self.playing:
            self.playing = playing
            self.buttons["play"].setToolTip("Stop" if playing else "Play")
            self._set_icons()
        if playing and self.isVisible():
            self.timer.isActive() or self.timer.start()
        else:
            self.timer.stop()
        if count < 2:
            return
        self.title.setText("Frames" if cmd.get_movie_length(images=0) else "States")
        self.speed.setText("%g fps" % cmd.get_setting_float("movie_fps"))
        # set only what changed: setting a spin box drops text being typed
        self.counter.blockSignals(True)
        if self.counter.maximum() != count:
            self.counter.setMaximum(count)
            self.counter.setSuffix("/%d" % count)
        if self.counter.value() != cmd.get_frame():
            self.counter.setValue(cmd.get_frame())
        self.counter.blockSignals(False)


class ObjectPanel(QtWidgets.QDockWidget):
    def __init__(self, cmd, parent=None):
        super().__init__("Objects", parent)
        self.setObjectName("dock_objects")
        self.cmd = cmd
        self.rows = []
        self.items = {}

        self.tree = tree = QtWidgets.QTreeWidget()
        tree.setHeaderLabels(["Name"] + [letter for letter, _ in BUTTONS])
        tree.setUniformRowHeights(True)
        tree.setHeaderHidden(True)
        # visibility is a dot, lit when the object is shown
        tree.setStyleSheet(
            "QTreeView::item { padding: 4px 0; }"
            "QTreeView::indicator { width: 10px; height: 10px; border: none;"
            " border-radius: 5px; margin: 0 4px; background-color: palette(mid); }"
            "QTreeView::indicator:checked { background-color: #3ddc5a; }"
            # the dimmed name of a hidden object lights up under the mouse
            "QTreeView::item:hover { color: palette(text); }"
        )
        tree.viewport().setAttribute(Qt.WidgetAttribute.WA_Hover)
        tree.setExpandsOnDoubleClick(False)
        tree.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection
        )
        header = tree.header()
        header.setStretchLastSection(False)
        header.setMinimumSectionSize(tree.fontMetrics().height() * 3 // 2)
        header.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.Stretch)
        for column in range(1, len(BUTTONS) + 1):
            header.setSectionResizeMode(
                column, QtWidgets.QHeaderView.ResizeMode.ResizeToContents
            )
            tree.setItemDelegateForColumn(column, ChipDelegate(tree))
        tree.itemClicked.connect(self._clicked)
        tree.itemDoubleClicked.connect(self._double_clicked)
        tree.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        tree.customContextMenuRequested.connect(self._context_menu)
        tree.viewport().installEventFilter(self)

        copy = QtGui.QShortcut(QtGui.QKeySequence.StandardKey.Copy, tree)
        copy.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        copy.activated.connect(self.copy_selected)
        space = QtGui.QShortcut(Qt.Key.Key_Space, tree)
        space.setContext(Qt.ShortcutContext.WidgetShortcut)
        space.activated.connect(
            lambda: tree.currentItem() and self._toggle(tree.currentItem())
        )

        self.filter = QtWidgets.QLineEdit()
        self.filter.setPlaceholderText("Filter by name")
        self.filter.setClearButtonEnabled(True)
        self.filter.textChanged.connect(self._apply_filter)

        self.states = StatesPanel(cmd)
        # toggle bar: the window fills it with its buttons
        self.footer = QtWidgets.QHBoxLayout()
        self.footer.setSpacing(2)
        self.wizard = WizardPanel(cmd)
        self.wizard.changed.connect(self.refresh_soon)

        body = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(body)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(self.filter)
        layout.addWidget(tree, 1)
        layout.addWidget(self.wizard)
        layout.addWidget(self.states)
        layout.addLayout(self.footer)
        self.setWidget(body)

    # ponytail: polled from the main window's 0.5 s timer (names, types and
    # visibility only, no atom data); upgrade path is a change notification
    # from the C++ executive.
    def poll(self):
        self.wizard.refresh()
        if self.wizard.panel:
            self.show()  # a wizard is unusable without its buttons
        if self.isVisible():
            self.refresh()
            self.states.refresh()

    def refresh_soon(self):
        QtCore.QTimer.singleShot(0, self.poll)

    def refresh(self):
        rows = snapshot(self.cmd)
        if rows == self.rows:
            return
        rebuilt = [(r[0], r[1], r[3]) for r in rows] != [
            (r[0], r[1], r[3]) for r in self.rows
        ]
        self.rows = rows
        if rebuilt:
            self._rebuild(rows)
        dim = self.palette().brush(
            QtGui.QPalette.ColorGroup.Disabled, QtGui.QPalette.ColorRole.Text
        )
        for name, _, enabled, _ in rows:
            item = self.items[name]
            item.setCheckState(0, Checked if enabled else Unchecked)
            item.setForeground(0, QtGui.QBrush() if enabled else dim)
        self.tree.viewport().update()  # chips follow the visibility
        if rebuilt:
            self._apply_filter()

    def _rebuild(self, rows):
        """Recreate the items; keep collapsed groups and the row selection."""
        collapsed = {
            n for n, i in self.items.items() if i.childCount() and not i.isExpanded()
        }
        selected = {n for n, i in self.items.items() if i.isSelected()}
        self.tree.clear()
        self.items = {}
        for name, kind, _, _ in rows:
            item = QtWidgets.QTreeWidgetItem(
                ["(%s)" % name if kind == "selection" else name]
            )
            item.setData(0, UserRole, (name, kind))
            # the whole row toggles in _clicked, not only the dot
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsUserCheckable)
            for column, function in enumerate(MENUS.get(kind, MENUS["object:"]), 1):
                if function:
                    letter, tip = BUTTONS[column - 1]
                    item.setText(column, letter)
                    item.setToolTip(column, tip)
                    item.setTextAlignment(column, Qt.AlignmentFlag.AlignCenter)
            self.items[name] = item
        root = self.tree.invisibleRootItem()
        for name, _, _, parent in rows:
            self.items.get(parent, root).addChild(self.items[name])
        for name, item in self.items.items():
            item.setExpanded(name not in collapsed)
            item.setSelected(name in selected)

    def _apply_filter(self):
        text = self.filter.text().lower()

        def show(item):
            # a group stays visible while any member matches
            children = [show(item.child(i)) for i in range(item.childCount())]
            visible = text in item.data(0, UserRole)[0].lower() or any(children)
            item.setHidden(not visible)
            return visible

        for i in range(self.tree.topLevelItemCount()):
            show(self.tree.topLevelItem(i))

    def _toggle(self, item):
        """Show or hide the object of this row in the viewport."""
        name = item.data(0, UserRole)[0]
        if item.checkState(0) == Checked:
            self.cmd.disable(name)
        else:
            self.cmd.enable(name)
        self.refresh()

    def _clicked(self, item, column):
        if column:
            rect = self.tree.visualRect(self.tree.indexFromItem(item, column))
            self._popup(
                item, column - 1, self.tree.viewport().mapToGlobal(rect.bottomLeft())
            )
        # Ctrl or Shift click only selects rows, for copying
        elif not QtWidgets.QApplication.keyboardModifiers():
            self.tree.clearSelection()
            self._toggle(item)

    def _double_clicked(self, item, column):
        if column == 0:  # its first click toggled the row: show it, then zoom
            name = item.data(0, UserRole)[0]
            self.cmd.enable(name)
            self.cmd.zoom(name, animate=-1)
            self.refresh()

    def _popup(self, item, button, position):
        entries = menu_entries(self.cmd, *item.data(0, UserRole), button)
        if entries:
            popup = QtWidgets.QMenu(self.tree)
            fill_menu(popup, entries, self.cmd)
            popup.exec(position)
            self.refresh_soon()

    def _context_menu(self, position):
        """Right click or the Menu key on a row."""
        item = self.tree.itemAt(position) or self.tree.currentItem()
        if item is None:
            return
        popup = context_menu(self.tree, self.cmd, *item.data(0, UserRole))
        popup.exec(self.tree.viewport().mapToGlobal(position))
        self.refresh_soon()

    def eventFilter(self, watched, event):
        # counts are looked up when asked for, never while polling
        if event.type() == QtCore.QEvent.Type.ToolTip:
            item = self.tree.itemAt(event.pos())
            if item and self.tree.columnAt(event.pos().x()) == 0:
                kind, atoms, states = describe(self.cmd, *item.data(0, UserRole))
                parts = [kind]
                if atoms is not None:
                    parts.append("%d atoms" % atoms)
                if states is not None:
                    parts.append("%d states" % states)
                QtWidgets.QToolTip.showText(event.globalPos(), ", ".join(parts))
                return True
        return super().eventFilter(watched, event)

    def copy_selected(self):
        lines = ["Name\tType\tShown\tAtoms\tStates"]
        for name, kind, enabled, _ in self.rows:
            if self.items[name].isSelected():
                row = (name, *describe(self.cmd, name, kind))
                row = row[:2] + (enabled,) + row[2:]
                lines.append("\t".join("" if v is None else str(v) for v in row))
        QtWidgets.QApplication.clipboard().setText("\n".join(lines))
