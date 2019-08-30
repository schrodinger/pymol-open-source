import os

from pymol.Qt import QtGui, QtCore, QtWidgets
from pymol.Qt.utils import UpdateLock, PopupOnException
import pymol

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

def props_dialog(parent):  #noqa
    from pymol.setting import name_dict

    cmd = parent.cmd
    form = parent.load_form('props', 'floating')
    parent.addDockWidget(QtCore.Qt.RightDockWidgetArea, form._dialog)

    def make_entry(parent, label):
        item = QtWidgets.QTreeWidgetItem(parent)
        item.setText(0, str(label))
        item.setFlags(QtCore.Qt.ItemIsEditable |
                      QtCore.Qt.ItemIsEnabled |
                      QtCore.Qt.ItemIsSelectable )
        return item

    def make_cat(parent, label):
        item = QtWidgets.QTreeWidgetItem(parent)
        item.setText(0, str(label))
        item.setFirstColumnSpanned(True)
        item.setExpanded(True)
        item.setChildIndicatorPolicy(
            QtWidgets.QTreeWidgetItem.ShowIndicator)
        return item

    # make first column uneditable
    form.treeWidget.setItemDelegateForColumn(0, UneditableDelegate())

    item_object = make_cat(form.treeWidget, "Object-Level")
    item_object_ttt = make_entry(item_object, "TTT Matrix")
    item_object_settings = make_cat(item_object, "Settings")

    item_ostate = make_cat(form.treeWidget, "Object-State-Level")
    item_ostate_title = make_entry(item_ostate, "Title")
    item_ostate_matrix = make_entry(item_ostate, "State Matrix")
    item_ostate_settings = make_cat(item_ostate, "Settings")
    item_ostate_properties = Ellipsis  # Incentive PyMOL only

    item_atom = make_cat(form.treeWidget, "Atom-Level")
    item_atom_identifiers = make_cat(item_atom, "Identifiers")
    item_atom_builtins = make_cat(item_atom, "Properties (built-in)")
    item_atom_settings = make_cat(item_atom, "Settings")
    item_atom_properties = Ellipsis  # Incentive PyMOL only

    item_astate = make_cat(form.treeWidget, "Atom-State-Level")
    item_astate_builtins = make_cat(item_astate, "Properties (built-in)")
    item_astate_settings = make_cat(item_astate, "Settings")

    keys_atom_identifiers = ['model', 'index', 'segi', 'chain', 'resi',
                             'resn', 'name', 'alt', 'ID', 'rank']
    keys_atom_builtins = ['elem', 'q', 'b', 'type', 'formal_charge',
                          'partial_charge', 'numeric_type', 'text_type',
                          # avoid stereo auto-assignment errors
                          # 'stereo',
                          'vdw', 'ss', 'color', 'reps',
                          'protons', 'geom', 'valence', 'elec_radius']
    keys_astate_builtins = ['state', 'x', 'y', 'z']

    items = {}
    for key in keys_atom_identifiers:
        items[key] = make_entry(item_atom_identifiers, key)
    for key in keys_atom_builtins:
        items[key] = make_entry(item_atom_builtins, key)
    for key in keys_astate_builtins:
        items[key] = make_entry(item_astate_builtins, key)

    items['model'].setDisabled(True)
    items['index'].setDisabled(True)
    items['state'].setDisabled(True)

    def item_changed(item, column):
        """
        Edits current item.
        """

        if item_changed.skip:
            return

        model = form.input_model.currentText()
        state = form.input_state.value()
        key = item.text(0)
        new_value = item.text(1)
        parent = item.parent()

        result = False
        if item is item_object_ttt:
            try:
                if new_value:
                    result = cmd.set_object_ttt(model, new_value)
            except (ValueError, IndexError):
                result = False
        elif item is item_ostate_title:
            result = cmd.set_title(model, state, new_value)
        elif item is item_ostate_matrix:
            cmd.matrix_reset(model, state)
            try:
                new_value = cmd.safe_eval(new_value)
                result = cmd.transform_object(model, new_value, state)
            except: # CmdTransformObject-DEBUG: bad matrix
                result = False
        elif parent is item_object_settings:
            with PopupOnException():
                cmd.set(key, new_value, model, quiet=0)
        elif parent is item_ostate_settings:
            with PopupOnException():
                cmd.set(key, new_value, model, state, quiet=0)
        elif parent is item_ostate_properties:
            cmd.set_property(key, new_value, model, state, quiet=0)
        else:
            is_state = False

            if parent is item_atom_properties:
                key = 'p.' + key
            elif parent is item_atom_settings:
                key = 's.' + key
            elif parent is item_astate_settings:
                key = 's.' + key
                is_state = True
            elif key in keys_astate_builtins:
                is_state = True

            def get_new_value(old_value):
                if isinstance(old_value, (tuple, list)):
                    return cmd.safe_eval(new_value)

                try:
                    # cast to old type (required for e.g. 'resv = "3"')
                    return type(old_value)(new_value)
                except ValueError:
                    # return str and let PyMOL handle it (e.g. change
                    # type of user property)
                    return new_value.encode('ascii')

            alter_args = ('pk1', key + '= get_new_value(' + key + ')', 0,
                    {'get_new_value': get_new_value})

            if is_state:
                result = cmd.alter_state(state, *alter_args)
            else:
                result = cmd.alter(*alter_args)

        if not result:
            update_treewidget_model()

    item_changed.skip = False

    def update_object_settings(parent, model, state):
        parent.takeChildren()
        for sitem in (cmd.get_object_settings(model, state) or []):
            key = name_dict.get(sitem[0], sitem[0])
            item = make_entry(parent, key)
            item.setText(1, str(sitem[2]))

    def update_atom_settings(wrapper, parent):
        parent.takeChildren()
        for key in wrapper:
            item = make_entry(parent, name_dict.get(key, key))
            value = wrapper[key]
            item.setText(1, str(value))

    def update_atom_properties(wrapper, parent):
        parent.takeChildren()
        for key in wrapper:
            item = make_entry(parent, key)
            value = wrapper[key]
            item.setText(1, str(value))

    def update_atom_fields(ns):
        for key in keys_atom_identifiers + keys_atom_builtins:
            try:
                value = ns[key]
            except Exception as e:
                value = 'ERROR: ' + str(e)
            items[key].setText(1, str(value))
        update_atom_settings(ns['s'], item_atom_settings)

    def update_astate_fields(ns):
        for key in keys_astate_builtins:
            value = ns[key]
            items[key].setText(1, str(value))
        update_atom_settings(ns['s'], item_astate_settings)

    space = {
        'update_atom_fields': update_atom_fields,
        'update_astate_fields': update_astate_fields,
        'locals': locals,
    }

    def update_from_pk1():
        pk1_atom = []
        if cmd.iterate('?pk1', 'pk1_atom[:] = [model, index]', space=locals()) > 0:
            form.input_model.setCurrentIndex(form.input_model.findText(pk1_atom[0]))
            form.input_index.setValue(pk1_atom[1])

    def update_pk1():
        model = form.input_model.currentText()
        index = form.input_index.value()

        if model and index:
            try:
                cmd.edit((model, index))
                return True
            except pymol.CmdException:
                pass

        return False

    def update_treewidget(*args):
        if not update_pk1():
            return

        state = form.input_state.value()

        item_changed.skip = True
        count = cmd.iterate(
            'pk1', 'update_atom_fields(locals())', space=space)
        item_atom.setDisabled(not count)
        if count:
            count = cmd.iterate_state(
                state, 'pk1',
                'update_astate_fields(locals())', space=space)
        item_astate.setDisabled(not count)
        item_changed.skip = False

    def update_treewidget_state(*args):
        model = form.input_model.currentText()
        state = form.input_state.value()
        if not (model and state):
            return

        item_changed.skip = True
        item_ostate_title.setText(
            1, str(cmd.get_title(model, state) or ''))
        item_ostate_matrix.setText(
            1, str(
                cmd.get_object_matrix(
                    model, state, 0) or ''))

        update_object_settings(item_ostate_settings, model, state)

        update_treewidget()
        item_changed.skip = False

    @suspendable
    def update_treewidget_model(*args):
        item_changed.skip = True
        model = form.input_model.currentText()

        if not model:
            return

        item_object_ttt.setText(1, str(cmd.get_object_ttt(model) or ''))
        update_object_settings(item_object_settings, model, 0)

        natoms = cmd.count_atoms('?' + model)
        nstates = cmd.count_states('?' + model)

        form.input_state.setMinimum(1)
        form.input_state.setMaximum(nstates)
        form.input_index.setMinimum(1)
        form.input_index.setMaximum(natoms)

        item_atom.setHidden(natoms == 0)
        item_astate.setHidden(natoms == 0)
        item_ostate.setHidden(nstates == 0)

        update_treewidget_state()
        item_changed.skip = False

    def update_model_list(*args):
        if 'pk1' not in cmd.get_names('selections'):
            update_pk1()

        with update_treewidget_model.suspend:
            form.input_model.clear()
            form.input_model.addItems(get_object_names(cmd))
            update_from_pk1()

        update_treewidget_model()

    # init input fields
    form.input_model.addItems(get_object_names(cmd))
    form.input_state.setValue(cmd.get_state())

    # select pk1 atom if available
    update_from_pk1()

    # hook up events
    form.input_model.currentIndexChanged.connect(update_treewidget_model)
    form.input_state.valueChanged.connect(update_treewidget_state)
    form.input_index.valueChanged.connect(update_treewidget)
    form.button_refresh.clicked.connect(update_model_list)

    # themed icons only available by default on X11
    if form.button_refresh.icon().isNull():
        form.button_refresh.setIcon(QtGui.QIcon(
            os.path.expandvars('$PYMOL_DATA/pmg_qt/icons/refresh.svg')))

    # update and show
    update_treewidget_model()
    form.treeWidget.setColumnWidth(0, 200)

    form.treeWidget.itemChanged.connect(item_changed)

    return form._dialog
