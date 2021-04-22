import os

from pymol.Qt import QtGui, QtCore, QtWidgets
from pymol.Qt.utils import UpdateLock, PopupOnException
import pymol
Qt = QtCore.Qt

from pymol.setting import name_dict

# Functions to convert atom property values to strings for display in the GUI
strfunctions = {
    'color': lambda c: hex(c) if (c >= 0x40000000) else str(c),
    'reps': bin,
    'flags': bin,
}

class UneditableDelegate(QtWidgets.QStyledItemDelegate):
    def createEditor(self, parent, option, index):
        return None

class FunctionSuspender:
    def __init__(self, func):
        self.func = func
    def __enter__(self):
        self.func.suspended = True
    def __exit__(self, exc_type, exc_value, traceback):
        self.func.suspended = False

def suspendable(func):
    def wrapper(*args, **kwargs):
        if not func.suspended:
            return func(*args, **kwargs)
    func.suspended = False
    wrapper.suspend = FunctionSuspender(func)
    return wrapper

def get_object_names(_self):
    # was _self.get_names('public_objects') in PyMOL 2.1
    # but that throws exceptions for groups/isosurfaces/etc.
    names = _self.get_object_list()
    return names

class PropsDialog(QtWidgets.QWidget):
    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self)
        self.cmd = parent.cmd
        self.form = parent.load_form('props')
        self.setup_tree_widget()
        self.setup_behavior()
        self.item_changed_skip = False

    def make_entry(self, parent, label):
        item = QtWidgets.QTreeWidgetItem(parent)
        item.setText(0, str(label))
        item.setFlags(QtCore.Qt.ItemIsEditable |
                      QtCore.Qt.ItemIsEnabled |
                      QtCore.Qt.ItemIsSelectable )
        return item

    def make_cat(self, parent, label):
        item = QtWidgets.QTreeWidgetItem(parent)
        item.setText(0, str(label))
        item.setFirstColumnSpanned(True)
        item.setExpanded(True)
        item.setChildIndicatorPolicy(
            QtWidgets.QTreeWidgetItem.ShowIndicator)
        return item

    def setup_tree_widget(self):
        self.form.treeWidget.setItemDelegateForColumn(0, UneditableDelegate())

        self.item_object = self.make_cat(self.form.treeWidget, "Object-Level")
        self.item_object_ttt = self.make_entry(self.item_object, "TTT Matrix")
        self.item_object_settings = self.make_cat(self.item_object, "Settings")

        self.item_ostate = self.make_cat(self.form.treeWidget, "Object-State-Level")
        self.item_ostate_title = self.make_entry(self.item_ostate, "Title")
        self.item_ostate_matrix = self.make_entry(self.item_ostate, "State Matrix")
        self.item_ostate_settings = self.make_cat(self.item_ostate, "Settings")
        self.item_ostate_properties = Ellipsis  # Incentive PyMOL only

        self.item_atom = self.make_cat(self.form.treeWidget, "Atom-Level")
        self.item_atom_identifiers = self.make_cat(self.item_atom, "Identifiers")
        self.item_atom_builtins = self.make_cat(self.item_atom, "Properties (built-in)")
        self.item_atom_settings = self.make_cat(self.item_atom, "Settings")
        self.item_atom_properties = Ellipsis  # Incentive PyMOL only

        self.item_astate = self.make_cat(self.form.treeWidget, "Atom-State-Level")
        self.item_astate_builtins = self.make_cat(self.item_astate, "Properties (built-in)")
        self.item_astate_settings = self.make_cat(self.item_astate, "Settings")

        self.keys_atom_identifiers = ['model', 'index', 'segi', 'chain', 'resi',
                                'resn',
                                'oneletter',  # read-only
                                'name', 'alt', 'ID', 'rank']
        self.keys_atom_builtins = ['elem', 'q', 'b', 'type', 'formal_charge',
                            'partial_charge', 'numeric_type', 'text_type',
                            # avoid stereo auto-assignment errors
                            # 'stereo',
                            'vdw', 'ss', 'color', 'reps',
                            'flags',
                            'label', 'cartoon',
                            'protons', 'geom', 'valence', 'elec_radius']
        self.keys_astate_builtins = ['state', 'x', 'y', 'z']

        self.items = {}
        for key in self.keys_atom_identifiers:
            self.items[key] = self.make_entry(self.item_atom_identifiers, key)
        for key in self.keys_atom_builtins:
            self.items[key] = self.make_entry(self.item_atom_builtins, key)
        for key in self.keys_astate_builtins:
            self.items[key] = self.make_entry(self.item_astate_builtins, key)

        self.items['model'].setDisabled(True)
        self.items['index'].setDisabled(True)
        self.items['state'].setDisabled(True)
        self.items['oneletter'].setDisabled(True)

    def setup_behavior(self):
        # init input fields
        self.form.input_model.addItems(get_object_names(self.cmd))
        self.form.input_state.setValue(self.cmd.get_state())

        # select pk1 atom if available
        self.update_from_pk1()

        # Install an eventfilter
        self.form.treeWidget.installEventFilter(self)

        # hook up events
        self.form.input_model.currentIndexChanged.connect(self.update_treewidget_model)
        self.form.input_state.valueChanged.connect(self.update_treewidget_state)
        self.form.input_index.valueChanged.connect(self.update_treewidget)
        self.form.button_refresh.clicked.connect(self.update_model_list)

        # themed icons only available by default on X11
        if self.form.button_refresh.icon().isNull():
            self.form.button_refresh.setIcon(QtGui.QIcon(
                os.path.expandvars('$PYMOL_DATA/pmg_qt/icons/refresh.svg')))

        # update and show
        self.update_treewidget_model()
        self.form.treeWidget.setColumnWidth(0, 200)

        self.form.treeWidget.itemChanged.connect(self.item_changed)

    def get_dialog(self):
        return self.form._dialog

    def item_changed(self, item, column):
        """
        Edits current item.
        """

        if self.item_changed_skip:
            return

        model = self.form.input_model.currentText()
        state = self.form.input_state.value()
        key = item.text(0)
        new_value = item.text(1)
        parent = item.parent()
        result = False

        if not result:
            if item is self.item_object_ttt:
                try:
                    if new_value:
                        result = self.cmd.set_object_ttt(model, new_value)
                except (ValueError, IndexError):
                    result = False
            elif item is self.item_ostate_title:
                result = self.cmd.set_title(model, state, new_value)
            elif item is self.item_ostate_matrix:
                self.cmd.matrix_reset(model, state)
                try:
                    new_value = self.cmd.safe_eval(new_value)
                    result = self.cmd.transform_object(model, new_value, state)
                except: # CmdTransformObject-DEBUG: bad matrix
                    result = False
            elif parent is self.item_object_settings:
                with PopupOnException():
                    self.cmd.set(key, new_value, model, quiet=0)
            elif parent is self.item_ostate_settings:
                with PopupOnException():
                    self.cmd.set(key, new_value, model, state, quiet=0)
            elif parent is self.item_ostate_properties:
                self.cmd.set_property(key, new_value, model, state, quiet=0)
            else:
                is_state = False

                if parent is self.item_atom_properties:
                    key = 'p.' + key
                elif parent is self.item_atom_settings:
                    key = 's.' + key
                elif parent is self.item_astate_settings:
                    key = 's.' + key
                    is_state = True
                elif key in self.keys_astate_builtins:
                    is_state = True

                def get_new_value(old_value):
                    if isinstance(old_value, (tuple, list, bool)):
                        return self.cmd.safe_eval(new_value)
                    try:
                        # cast to old type (required for e.g. 'resv = "3"')
                        if isinstance(old_value, int):
                            return int(new_value, 0)
                        return type(old_value)(new_value)
                    except ValueError:
                        # return str and let PyMOL handle it (e.g. change
                        # type of user property)
                        return new_value.encode('utf-8')

                alter_args = ('pk1', key + '= get_new_value(' + key + ')', 0,
                        {'get_new_value': get_new_value})

                with PopupOnException():
                    if is_state:
                        result = self.cmd.alter_state(state, *alter_args)
                    else:
                        result = self.cmd.alter(*alter_args)

        if not result:
            self.update_treewidget_model()

        self.item_changed_skip = False

    def unset_item(self, item):
        '''
        Calls unset on item if applicable. Matrices will be reset.
        '''
        model = self.form.input_model.currentText()
        state = self.form.input_state.value()
        key = item.text(0)
        old_value = item.text(1)
        old_value_type = type(old_value)
        parent = item.parent()

        unset_result = True

        if old_value == "":
            unset_result = False
            return unset_result

        if item is self.item_object_ttt:
            self.cmd.matrix_reset(model, mode=1) # mode=1 for TTT matrix reset
        elif item is self.item_ostate_title:
            self.cmd.set_title(model, state, '')
        elif item is self.item_ostate_matrix:
            self.cmd.matrix_reset(model, state, mode=2) # mode=2 for state matrix reset
        elif parent is self.item_object_settings:
            self.cmd.unset(key, model, quiet=0)
        elif parent is self.item_ostate_settings:
            self.cmd.unset(key, model, state, quiet=0)
        elif parent is self.item_ostate_properties:
            self.cmd.set_property(key, None, model, state, quiet=0)
        elif parent is self.item_atom_properties:
            key = 'p.' + key
            alter_args = ('pk1', key + '= None', 0)
            unset_result = self.cmd.alter(*alter_args)
        elif parent is self.item_atom_settings:
            key = 's.' + key
            alter_args = ('pk1', key + '= None', 0)
            unset_result = self.cmd.alter(*alter_args)
        else:
            unset_result = False

        if unset_result:
            self.update_treewidget_model()

        return unset_result

    def unset_caller(self):
        item_list_selected = self.form.treeWidget.selectedItems()
        self.unset_item(item_list_selected[0])

    def eventFilter(self, source, event):
        '''
        Event filter for creating new shortcuts. Processes the key event before passing it on.
        '''
        if (event.type() == QtCore.QEvent.KeyPress and source is self.form.treeWidget):
            if (event.key() == QtCore.Qt.Key_Delete):
                self.unset_caller()
                return 0
        return super().eventFilter(source, event)

    def update_object_settings(self, parent, model, state):
        parent.takeChildren()
        for sitem in (self.cmd.get_object_settings(model, state) or []):
            key = name_dict.get(sitem[0], sitem[0])
            item = self.make_entry(parent, key)
            item.setText(1, str(sitem[2]))

    def update_atom_settings(self, wrapper, parent):
        parent.takeChildren()
        for key in wrapper:
            item = self.make_entry(parent, name_dict.get(key, key))
            value = wrapper[key]
            item.setText(1, str(value))

    def update_atom_properties(self, wrapper, parent):
        parent.takeChildren()
        for key in wrapper:
            item = self.make_entry(parent, key)
            value = wrapper[key]
            item.setText(1, str(value))

    def update_atom_fields(self, ns):
        for key in self.keys_atom_identifiers + self.keys_atom_builtins:
            try:
                value = ns[key]
            except Exception as e:
                value = 'ERROR: ' + str(e)
            self.items[key].setText(1, strfunctions.get(key, str)(value))
        self.update_atom_settings(ns['s'], self.item_atom_settings)

    def update_astate_fields(self, ns):
        for key in self.keys_astate_builtins:
            value = ns[key]
            self.items[key].setText(1, str(value))
        self.update_atom_settings(ns['s'], self.item_astate_settings)

    def update_from_pk1(self):
        pk1_atom = []
        if self.cmd.iterate('?pk1', 'pk1_atom[:] = [model, index]', space=locals()) > 0:
            self.form.input_model.setCurrentIndex(self.form.input_model.findText(pk1_atom[0]))
            self.form.input_index.setValue(pk1_atom[1])

    def update_pk1(self):
        model = self.form.input_model.currentText()
        index = self.form.input_index.value()

        if model and index:
            try:
                self.cmd.edit((model, index))
                return True
            except pymol.CmdException:
                pass

        return False

    def update_treewidget(self, *args):
        if not self.update_pk1():
            return

        state = self.form.input_state.value()

        self.item_changed_skip = True
        space = {'update_atom_fields': self.update_atom_fields,'update_astate_fields': self.update_astate_fields,'locals': locals}
        count = self.cmd.iterate(
            'pk1', 'update_atom_fields(locals())', space=space)
        self.item_atom.setDisabled(not count)
        if count:
            count = self.cmd.iterate_state(
                state, 'pk1',
                'update_astate_fields(locals())', space=space)
        self.item_astate.setDisabled(not count)
        self.item_changed_skip = False

    def update_treewidget_state(self, *args):
        model = self.form.input_model.currentText()
        state = self.form.input_state.value()
        if not (model and state):
            return

        self.item_changed_skip = True
        self.item_ostate_title.setText(
            1, str(self.cmd.get_title(model, state) or ''))
        self.item_ostate_matrix.setText(
            1, str(
                self.cmd.get_object_matrix(
                    model, state, 0) or ''))

        self.update_object_settings(self.item_ostate_settings, model, state)

        self.update_treewidget()
        self.item_changed_skip = False

    @suspendable
    def update_treewidget_model(self, *args):
        self.item_changed_skip = True
        model = self.form.input_model.currentText()

        if not model:
            return

        self.item_object_ttt.setText(1, str(self.cmd.get_object_ttt(model) or ''))
        self.update_object_settings(self.item_object_settings, model, 0)

        natoms = self.cmd.count_atoms('?' + model)
        nstates = self.cmd.count_states('?' + model)

        self.form.input_state.setMinimum(1)
        self.form.input_state.setMaximum(nstates)
        self.form.input_index.setMinimum(1)
        self.form.input_index.setMaximum(natoms)

        self.item_atom.setHidden(natoms == 0)
        self.item_astate.setHidden(natoms == 0)
        self.item_ostate.setHidden(nstates == 0)

        self.update_treewidget_state()
        self.item_changed_skip = False

    def update_model_list(self, *args):
        if 'pk1' not in self.cmd.get_names('selections'):
            self.update_pk1()

        with self.update_treewidget_model.suspend:
            self.form.input_model.clear()
            self.form.input_model.addItems(get_object_names(self.cmd))
            self.update_from_pk1()

        self.update_treewidget_model()
