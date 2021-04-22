import json
import os
from textwrap import fill

from pymol import setting
from pymol import save_shortcut
from pymol.Qt import QtGui, QtWidgets
from pymol.Qt import QtCore, QtCoreModels
from pymol.shortcut_manager import ShortcutManager, ShortcutIndex
from pymol.keyboard import get_default_keys
Qt = QtCore.Qt
QSI = QtGui.QStandardItem  # For brevity


def get_shortcut_key_map():
    shortcut_key_map = {}
    for key, value in vars(Qt).items():
        if isinstance(value, Qt.Key):
            shortcut_key_map[value] = key.partition('_')[2]
    return shortcut_key_map


_SHORTCUT_KEY_MAP = get_shortcut_key_map()

_SHORTCUT_MODIFIER_MAP = {
    Qt.ControlModifier: _SHORTCUT_KEY_MAP[Qt.Key_Control],
    Qt.AltModifier: _SHORTCUT_KEY_MAP[Qt.Key_Alt],
    Qt.ShiftModifier: _SHORTCUT_KEY_MAP[Qt.Key_Shift],
    Qt.MetaModifier: _SHORTCUT_KEY_MAP[Qt.Key_Meta],
}

_REPLACE_KEYS = {
    'PageUp': 'pgup',
    'PageDown': 'pgdn',
    'Home': 'home',
    'Insert': 'insert',
    'Up': 'up',
    'Down': 'down',
    'Left': 'left',
    'Right': 'right',
    'End':'end'
}

class PyMOLShortcutMenu(QtWidgets.QWidget):
    '''
    Keyboard shortcut dialog for PyMOL. This displays all assigned shortcuts
    and allows them to be changed or new shortcuts to be created.
    '''

    def __init__(self, parent, saved_shortcuts, cmd):
        QtWidgets.QWidget.__init__(self, parent, Qt.Window)
        self.resize(700, 700)
        self.cmd = cmd
        self.shortcut_manager = ShortcutManager(saved_shortcuts, cmd)

        self.build_panel_elements(parent)

    def build_panel_elements(self, parent):
        '''
        Responsible for creating all panel elements in order and adding them to the layout.
        '''
        self.create_new_form = parent.load_form("create_shortcut", None)
        self.help_form = parent.load_form("help_shortcut", None)
        self.confirm_change = parent.load_form("change_confirm", None)

        self.model = QtGui.QStandardItemModel(self)
        self.proxy_model = QtCoreModels.QSortFilterProxyModel(self)
        self.proxy_model.setSourceModel(self.model)
        self.proxy_model.setFilterCaseSensitivity(Qt.CaseInsensitive)
        self.proxy_model.setFilterKeyColumn(-1)

        self.setWindowTitle('Keyboard Shortcut Menu')
        layout = QtWidgets.QVBoxLayout(self)
        self.setLayout(layout)

        # Create layout for filter bar and refresh button
        top_layout = QtWidgets.QGridLayout()
        layout.addLayout(top_layout)

        # Filter
        self.filter_le = QtWidgets.QLineEdit(self)
        top_layout.addWidget(self.filter_le)
        self.filter_le.setPlaceholderText("Filter")
        self.filter_le.textChanged.connect(self.proxy_model.setFilterRegExp)

        self.refresh_button = QtWidgets.QPushButton(self)
        self.refresh_button.resize(26, 26)
        top_layout.addWidget(self.refresh_button, 0, 1)
        # themed icons only available by default on X11
        if self.refresh_button.icon().isNull():
            self.refresh_button.setIcon(QtGui.QIcon(
                os.path.expandvars('$PYMOL_DATA/pmg_qt/icons/refresh.svg')))
        self.refresh_button.setToolTip(
            "Refresh the table to reflect any external changes")
        self.refresh_button.clicked.connect(self.refresh_populate)

        # Table
        self.table = QtWidgets.QTableView(self)
        self.table.setModel(self.proxy_model)
        layout.addWidget(self.table)
        self.intial_populate()
        self.formatTable()

        # Add layout for buttons
        button_layout = QtWidgets.QGridLayout()
        layout.addLayout(button_layout)

        # Buttons
        self.create_new_button = QtWidgets.QPushButton(self)
        button_layout.addWidget(self.create_new_button, 0, 0)
        self.create_new_button.setText("Create New")
        self.create_new_button.setToolTip(
            "Add a key binding that does not currently appear on the table")
        self.create_new_button.clicked.connect(
            lambda: self.create_new_form._dialog.show())

        self.delete_selected_button = QtWidgets.QPushButton(self)
        button_layout.addWidget(self.delete_selected_button, 0, 1)
        self.delete_selected_button.setText("Delete Selected")
        self.delete_selected_button.setToolTip(
            "Unbind selected key bindings and remove any that have been created")
        self.delete_selected_button.clicked.connect(self.delete_selected)
        self.delete_selected_button.setEnabled(False)

        self.reset_selected_button = QtWidgets.QPushButton(self)
        button_layout.addWidget(self.reset_selected_button, 0, 2)
        self.reset_selected_button.setText("Reset Selected")
        self.reset_selected_button.setToolTip(
            "Restore selected key bindings to their default values")
        self.reset_selected_button.clicked.connect(self.reset_selected)
        self.reset_selected_button.setEnabled(False)

        self.reset_all_button = QtWidgets.QPushButton(self)
        button_layout.addWidget(self.reset_all_button, 0, 3)
        self.reset_all_button.setText("Reset All")
        self.reset_all_button.setToolTip(
            "Restore all key bindings to their default values and remove any that have been created")
        self.reset_all_button.clicked.connect(self.reset_all_default)

        self.save_button = QtWidgets.QPushButton(self)
        button_layout.addWidget(self.save_button, 0, 4)
        self.save_button.setText("Save")
        self.save_button.setToolTip(
            "Save the current key bindings to be loaded automatically when opening PyMOL")
        self.save_button.clicked.connect(self.shortcut_manager.save_shortcuts)

        # Ensuring that confirmed key and binding remain in scope
        self.confirm_new_key = ''
        self.confirm_new_binding = ''

        # Connect create new and confirm menus
        self.create_new_shortcut_menu_connect()
        self.confirm_menu_connect()

        self.model.itemChanged.connect(self.itemChanged)

    def populateData(self):
        '''
        Fill the model with data from shortcut_dict.
        '''
        self.model.clear()
        self.model.setHorizontalHeaderLabels(
            ['Key', 'Command (click to edit)', 'Description'])

        for key, shortcut_list in self.shortcut_manager.cmd.shortcut_dict.items():
            key_item = QSI(key)
            command_item = QSI()
            descript_item = QSI()
            key_item.setFlags(Qt.ItemIsEnabled)
            descript_item.setFlags(Qt.ItemIsEditable)

            if shortcut_list[ShortcutIndex.USER_DEF]:
                if shortcut_list[ShortcutIndex.USER_DEF] != "Deleted":
                    command_text = shortcut_list[ShortcutIndex.USER_DEF]
                    descript_text = "user defined"
                else:
                    command_text = "Deleted"
                    descript_text = "Deleted"
            else:
                command_text = shortcut_list[ShortcutIndex.COMMAND]
                descript_text = shortcut_list[ShortcutIndex.DESCRIPT]

            command_item.setText(command_text)
            descript_item.setText(descript_text)

            self.model.appendRow([key_item, command_item, descript_item])
        self.formatTable()
        self.table.selectionModel().selectionChanged.connect(self.selection_changed)

    def selection_changed(self):
        item_selected = bool(self.table.selectionModel().selectedIndexes())
        self.delete_selected_button.setEnabled(item_selected)
        self.reset_selected_button.setEnabled(item_selected)

    def intial_populate(self):
        '''
        Runs when the menu is first opened. Separate from populateData so that
        the saved dictionary and current state of key_mappings can be checked.
        '''
        self.shortcut_manager.check_saved_dict()
        self.shortcut_manager.check_key_mappings()
        self.populateData()

    def refresh_populate(self):
        '''
        Called through the refresh button. 
        '''
        self.shortcut_manager.check_key_mappings()
        self.populateData()

    def reset_all_default(self):
        '''
        Iterates over all values to restore their default commands.
        This will restore them to the values from keyboard.py
        '''
        self.shortcut_manager.reset_all_default()
        self.populateData()

    def delete_selected(self):
        '''
        Removes selected keybindings and updates table to say "Deleted".
        Keys that don't have default values will be removed completely.
        '''
        selection_model = self.table.selectionModel()
        list_indexes = selection_model.selectedIndexes()
        delete_keys = []

        for ind, table_index in enumerate(list_indexes):
            table_colm_key = self.table.model().index(table_index.row(), 0)
            table_colm_command = self.table.model().index(table_index.row(), 1)
            table_colm_descipt = self.table.model().index(table_index.row(), 2)

            delete_key = table_colm_key.data()
            self.cmd.set_key(delete_key, '')
            if delete_key in self.shortcut_manager.default_bindings:
                self.table.model().setData(table_colm_command, 'Deleted')
                self.table.model().setData(table_colm_descipt, 'Deleted')
                self.shortcut_manager.cmd.shortcut_dict[delete_key][ShortcutIndex.USER_DEF] = 'Deleted'
            else:
                self.model.removeRow(table_index.row())
                print(delete_key, " has been deleted and will be removed from the table")
                delete_keys.append(delete_key)

        for key in delete_keys:
            del self.shortcut_manager.cmd.shortcut_dict[key]

    def reset_selected(self):
        '''
        Restores default key bindings for items selected in the tables selection model.
        '''
        selection_model = self.table.selectionModel()
        list_indexes = selection_model.selectedIndexes()

        for ind, table_index in enumerate(list_indexes):
            table_colm_key = self.table.model().index(table_index.row(), 0)
            table_colm_command = self.table.model().index(table_index.row(), 1)
            table_colm_descipt = self.table.model().index(table_index.row(), 2)

            reset_key = table_colm_key.data()
            if reset_key not in self.shortcut_manager.default_bindings:
                print("This key does not have a default value.")
            else:
                reset_binding = self.shortcut_manager.default_bindings[reset_key]
                reset_command = self.shortcut_manager.cmd.shortcut_dict[
                    reset_key][ShortcutIndex.COMMAND]
                reset_description = self.shortcut_manager.cmd.shortcut_dict[
                    reset_key][ShortcutIndex.DESCRIPT]

                self.table.model().setData(table_colm_command, reset_command)
                self.table.model().setData(table_colm_descipt, reset_description)

                self.shortcut_manager.cmd.shortcut_dict[reset_key][ShortcutIndex.USER_DEF] = ''

                self.cmd.set_key(reset_key, reset_binding)

    def create_new_shortcut_menu_connect(self):
        self.create_new_form.createButton.clicked.connect(
            self.create_new_shortcut_caller)
        self.create_new_form.helpButton.clicked.connect(
            self.help_menu_shortcut)
        self.create_new_form.keyEdit.installEventFilter(self)
        self.create_new_form.helpButton.setDefault(False)
        self.create_new_form.helpButton.setAutoDefault(False)

    def eventFilter(self, source, event):
        '''
        Event filter for creating new shortcuts. Processes the key event before passing it on.
        '''
        if (event.type() == QtCore.QEvent.KeyPress and source is self.create_new_form.keyEdit):
            raw_string = self.keyevent_to_string(event)
            processed_string = self.process_keyevent_string(raw_string)

            if processed_string in self.shortcut_manager.reserved_keys:
                return 0

            if processed_string:
                self.create_new_form.keyEdit.setText(processed_string)

        return super().eventFilter(source, event)

    def keyevent_to_string(self, event):
        '''
        Generates string from captured key event for process_keyevent_string.
        '''
        keyevent_list = []
        for mod, event_text in _SHORTCUT_MODIFIER_MAP.items():
            if event.modifiers() & mod:
                keyevent_list.append(event_text)

        key = _SHORTCUT_KEY_MAP.get(event.key(), event.text())

        if key not in keyevent_list:
            keyevent_list.append(key)

        return ' '.join(keyevent_list)

    def process_keyevent_string(self, raw_string):
        '''
        Returns string of keyevent used to populate create key menu.
        Provides filtering to the keyevents, returning a list of the accepted strings.
        raw_string: string from keyevent_to_string 
        '''
        split_string = raw_string.split()
        process_list = []
        prefix_key = split_string[0]
        if len(split_string) >= 2:
            prefix_key = split_string[0]
            suffix_key = split_string[1]
            if (prefix_key == 'Control' or prefix_key == 'Meta') and split_string[1]:
                if suffix_key == 'Shift' and len(split_string) > 2:
                    process_list.append('CTSH')
                    suffix_key = split_string[2]
                else:
                    process_list.append('CTRL')
            elif prefix_key == 'Alt':
                process_list.append('ALT')
            elif prefix_key == 'Shift':
                process_list.append('SHFT')
            if suffix_key in _REPLACE_KEYS:
                suffix_key = _REPLACE_KEYS[suffix_key]
            process_list.append(suffix_key)
        elif prefix_key in _REPLACE_KEYS:
            process_list.append(_REPLACE_KEYS[prefix_key])
        return('-'.join(process_list))

    def help_menu_shortcut(self):
        self.help_form._dialog.show()

    def confirm_menu_connect(self):
        self.confirm_change.confirmButton.clicked.connect(
            lambda: self.shortcut_manager.create_new_shortcut(self.confirm_new_key, self.confirm_new_binding))
        self.confirm_change.confirmButton.clicked.connect(lambda: self.populateData())
        self.confirm_change.confirmButton.clicked.connect(
            lambda: self.confirm_change._dialog.hide())
        self.confirm_change.cancelButton.clicked.connect(
            lambda: self.confirm_change._dialog.hide())

    def create_new_shortcut_caller(self):
        '''
        Creates a new shortcut after checking existing and reserved keys.
        '''
        new_key = self.create_new_form.keyEdit.text()
        new_binding = self.create_new_form.commandEdit.text()

        if new_key == '':
            pass
        elif new_key in self.shortcut_manager.cmd.shortcut_dict:
            hide_confirm_menu = self.confirm_change.doNotShowCheckBox.isChecked()
            if not hide_confirm_menu:
                self.confirm_new_key = new_key
                self.confirm_new_binding = new_binding
                self.confirm_change._dialog.show()
            else:
                self.shortcut_manager.create_new_shortcut(new_key, new_binding)
                self.populateData()
        else:
            self.shortcut_manager.create_new_shortcut(new_key, new_binding)
            self.populateData()
            self.table.scrollToBottom()

    def formatTable(self):
        '''
        Set up the table to look appropriately
        '''
        hh = self.table.horizontalHeader()
        hh.setStretchLastSection(True)
        self.table.verticalHeader().setVisible(False)
        self.table.setFocus()
        self.table.hide()
        self.table.resizeColumnsToContents()
        self.table.show()

    def itemChanged(self, item):
        """
        Called every time an item in the table is changed, only command items
        are changeable.
        @param item: The item which has changed
        @type  item: QStandardItem
        """
        try:
            if item.column() == 1 and item.text() != "Deleted":
                changed_key = self.model.index(item.row(), 0).data()
                changed_index = self.table.model().index(item.row(), 0)

                self.cmd.set_key(changed_key, item.text())
                self.shortcut_manager.cmd.shortcut_dict[changed_key][2] = item.text()

                filter_active = bool(self.filter_le.text())
                if not filter_active:
                    self.table.model().setData(self.table.model().index(item.row(), 2), 'user defined')
                else:
                    pass
        except Exception as e:
            print(e)
            print("Failed to change key binding")
