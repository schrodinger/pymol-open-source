'''
PyMOL APBS GUI Plugin

(c) Schrodinger, Inc.
'''

import sys
import os
import importlib

from pymol.Qt import QtCore, QtWidgets
from pymol.Qt.utils import loadUi, AsyncFunc, MainThreadCaller

getOpenFileNames = QtWidgets.QFileDialog.getOpenFileNames

from . import electrostatics
from .qtwidgets import ResizableMessageBox as QMessageBox


class SilentAbort(Exception):
    pass


class StdOutCapture:
    '''
Redirect stdout and/or stderr to a temporary file until 'release' is called.
    '''

    def __init__(self, stdout=True, stderr=True):
        import tempfile
        self.temp = tempfile.TemporaryFile('w+')
        self.save_stdout = None
        self.save_stderr = None
        if stdout:
            self.save_stdout = sys.stdout
            sys.stdout = self.temp
        if stderr:
            self.save_stderr = sys.stderr
            sys.stderr = self.temp

    def release(self):
        if self.save_stdout is not None:
            sys.stdout = self.save_stdout
        if self.save_stderr is not None:
            sys.stderr = self.save_stderr
        self.temp.seek(0)
        content = self.temp.read()
        self.temp.close()
        return content


def run_impl(form, _self):
    '''
Execute the pipeline (prep, apbs, surface vis)
    '''
    selection = form.input_sele.currentText()
    prep_name = ''
    map_name = ''
    ramp_name = ''
    group_name = ['']

    def get_name(lineinput):
        name = lineinput.text()
        if not name:
            name = _self.get_unused_name(lineinput.placeholderText())
        return name

    def do_group(name):
        if form.check_no_group.isChecked():
            return
        if not group_name[0]:
            group_name[0] = _self.get_unused_name('run', 1)
        _self.group(group_name[0], name)

    if form.do_prepare.isChecked():
        from . import creating as pc

        prep_name = get_name(form.prep_name)
        method = form.prep_method.currentText()
        warnings = ''

        get_res_charge = lambda resn: ((-1) if resn in ('GLU', 'ASP') else
                (1) if resn in ('ARG', 'LYS') else 0)

        def get_cb_pseudocharge(resn, name):
            if name != 'CB':
                return 0
            return get_res_charge(resn)

        try:
            state = _self.get_selection_state(selection)
        except BaseException as e:
            state = -1

        if method == 'pdb2pqr':
            warnings = pc.pdb2pqr_cli(prep_name, selection,
                    options=form.pdb2pqr_args.text(),
                    quiet=0,
                    state=state,
                    preserve=form.check_preserve.isChecked(),
                    fixrna=form.pdb2pqr_fixrna.isChecked(),
                    _proclist=form._proclist,
                    exe=form.pdb2pqr_exe.text())
            if form.pdb2pqr_ignore_warnings.isChecked():
                warnings = ''
        elif method == 'protein_assign_charges_and_radii':
            _self.create(prep_name, selection)
            _self.util.protein_assign_charges_and_radii(prep_name)
        elif method.startswith('prepwizard'):
            pc.prepwizard(prep_name, selection,
                    options=form.prepwizard_args.text(),
                    _proclist=form._proclist,
                    quiet=0,
                    preserve=form.check_preserve.isChecked())
            _self.alter(prep_name, 'elec_radius = vdw')
        elif method.startswith('use formal'):
            _self.create(prep_name, selection)
            _self.alter(prep_name,
                        '(elec_radius, partial_charge) = (vdw, formal_charge)')
        elif method.startswith('use CB'):
            _self.create(prep_name, selection)
            _self.alter(
                prep_name,
                '(elec_radius, partial_charge) = (vdw, getpc(resn, name))',
                space={'getpc': get_cb_pseudocharge})
        elif method.startswith('use CA'):
            _self.create(prep_name, '(%s) & name CA' % (selection))
            _self.alter(
                prep_name,
                '(elec_radius, vdw, partial_charge) = (3.0, 3.0, getpc(resn))',
                space={'getpc': get_res_charge})
        elif method.startswith('use vdw'):
            _self.create(prep_name, selection)
            _self.alter(prep_name, 'elec_radius = vdw')
        else:
            raise ValueError('unknown method: ' + method)

        selection = prep_name
        form.input_sele.addItem(prep_name)
        do_group(prep_name)

        if warnings:
            @form._callInMainThread
            def result():
                msgbox = QMessageBox(QMessageBox.Question, 'Continue?',
                    method + ' emmitted warnings, do you want to continue?',
                    QMessageBox.Yes | QMessageBox.No , form._dialog)
                msgbox.setDetailedText(warnings)
                return msgbox.exec_()

            if result == QMessageBox.No:
                raise SilentAbort

    if form.do_apbs.isChecked():
        map_name = get_name(form.apbs_map)
        template = form.apbs_template.toPlainText()
        electrostatics.map_new_apbs(
            map_name,
            selection,
            grid=form.apbs_grid.value(),
            focus=form.focus_sele.text(),
            quiet=0,
            preserve=form.check_preserve.isChecked(),
            exe=form.apbs_exe.text(),
            _template=template,
            _proclist=form._proclist,
            assign=0)
        form.surf_map.addItem(map_name)
        do_group(map_name)

    if form.do_surface.isChecked():
        if not map_name:
            map_name = form.surf_map.currentText()
            if not map_name:
                raise ValueError('missing map')

        if not prep_name:
            prep_name = _self.get_object_list(selection)[0]

        ramp_name = get_name(form.surf_ramp)

        v = form.surf_range.value()
        _self.ramp_new(ramp_name, map_name, [-v, 0, v])
        do_group(ramp_name)

        sas = 'Accessible' in form.surf_type.currentText()
        _self.set('surface_ramp_above_mode', not sas, prep_name)
        _self.set('surface_solvent', sas, prep_name)
        _self.set('surface_color', ramp_name, prep_name)
        _self.show('surface', selection)


def dialog(_self=None):
    if _self is None:
        from pymol import cmd as _self

    dialog = QtWidgets.QDialog()
    uifile = os.path.join(os.path.dirname(__file__), 'apbs.ui')
    form = loadUi(uifile, dialog)
    form._dialog = dialog
    form._proclist = []

    def set_apbs_in(contents):
        form.apbs_template.setPlainText(contents.strip())

    # hide options widgets
    form.optarea_prep.setVisible(False)
    form.optarea_apbs.setVisible(False)
    form.optarea_surf.setVisible(False)
    form.optarea_other.setVisible(False)

    # pre-fill form with likely data
    names = _self.get_object_list()
    names += ['(' + n + ')' for n in _self.get_names('public_selections')]
    if names:
        form.input_sele.clear()
        form.input_sele.addItems([('polymer & ' + name)
                                  if _self.count_atoms('polymer & ' + name) > 0
                                  else name for name in names])
    form.surf_map.addItems(_self.get_names_of_type('object:map'))
    set_apbs_in(electrostatics.template_apbs_in)

    # executables
    from shutil import which
    form.apbs_exe.setText(electrostatics.find_apbs_exe() or 'apbs')
    form.pdb2pqr_exe.setText(
            which('pdb2pqr') or
            which('pdb2pqr30') or
            # acellera::htmd-pdb2pqr provides pdb2pqr_cli
            which('pdb2pqr_cli') or 'pdb2pqr')


    # for async panels
    form._callInMainThread = MainThreadCaller()
    run_impl_async = AsyncFunc(run_impl)

    # "Run" button callback
    def run():
        form.tabWidget.setEnabled(False)
        form.button_ok.clicked.disconnect()
        form.button_ok.clicked.connect(abort)
        form.button_ok.setText('Abort')

        form._capture = StdOutCapture()

        # detach from main thread
        run_impl_async(form, _self)

    # "Run" button "finally" actions (main thread)
    @run_impl_async.finished.connect
    def run_finally(args):
        _, exception = args

        form._proclist[:] = []
        stdout = form._capture.release()
        print(stdout)

        form.button_ok.setText('Run')
        form.button_ok.clicked.disconnect()
        form.button_ok.clicked.connect(run)
        form.button_ok.setEnabled(True)
        form.tabWidget.setEnabled(True)

        if exception is not None:
            handle_exception(exception, stdout)
            return

        quit_msg = "Finished with Success. Close the APBS dialog?"
        if QMessageBox.Yes == QMessageBox.question(
                form._dialog, 'Finished', quit_msg, QMessageBox.Yes,
                QMessageBox.No):
            form._dialog.close()

    def handle_exception(e, stdout):
        if isinstance(e, SilentAbort):
            return

        msg = str(e) or 'unknown error'
        msgbox = QMessageBox(QMessageBox.Critical, 'Error',
                             msg, QMessageBox.Close, form._dialog)
        if stdout.strip():
            msgbox.setDetailedText(stdout)
        msgbox.exec_()

    # "Abort" button callback
    def abort():
        form.button_ok.setEnabled(False)
        while form._proclist:
            p = form._proclist.pop()
            try:
                p.terminate()
                p.returncode = -15  # SIGTERM
            except OSError as e:
                print(e)

    # selection checker
    check_sele_timer = QtCore.QTimer()
    check_sele_timer.setSingleShot(True)

    # grid auto-value
    form.apbs_grid_userchanged = False
    form.apbs_grid.setStyleSheet('background: #ff6')

    @form.apbs_grid.editingFinished.connect
    def _():
        form.apbs_grid_userchanged = True
        form.apbs_grid.setStyleSheet('')

    @check_sele_timer.timeout.connect
    def _():
        has_props = ['no', 'no']

        def callback(partial_charge, elec_radius):
            if partial_charge: has_props[0] = 'YES'
            if elec_radius > 0: has_props[1] = 'YES'

        n = _self.iterate(
            form.input_sele.currentText(),
            'callback(partial_charge, elec_radius)',
            space={'callback': callback})

        # grid auto-value (keep map size in the order of 200x200x200)
        if n > 1 and not form.apbs_grid_userchanged:
            e = _self.get_extent(form.input_sele.currentText())
            volume = (e[1][0] - e[0][0]) * (e[1][1] - e[0][1]) * (e[1][2] - e[0][2])
            grid = max(0.5, volume ** 0.333 / 200.0)
            form.apbs_grid.setValue(grid)

        if n < 1:
            label = 'Selection is invalid'
            color = '#f66'
        elif has_props == ['YES', 'YES']:
            label = 'No preparation necessary, selection has charges and radii'
            form.do_prepare.setChecked(False)
            color = '#6f6'
        else:
            label = 'Selection needs preparation (partial_charge: %s, elec_radius: %s)' % tuple(
                has_props)
            form.do_prepare.setChecked(True)
            color = '#fc6'

        form.label_sele_has.setText(label)
        form.label_sele_has.setStyleSheet('background: %s; padding: 5' % color)

    check_sele_timer.start(0)

    @form.apbs_exe_browse.clicked.connect
    def _():
        fnames = getOpenFileNames(None, filter='apbs (apbs*);;All Files (*)')[0]
        if fnames:
            form.apbs_exe.setText(fnames[0])

    @form.pdb2pqr_exe_browse.clicked.connect
    def _():
        fnames = getOpenFileNames(None, filter='pdb2pqr (pdb2pqr*);;All Files (*)')[0]
        if fnames:
            form.pdb2pqr_exe.setText(fnames[0])

    # hook up events
    form.input_sele.currentIndexChanged.connect(
        lambda: check_sele_timer.start(0))
    form.input_sele.editTextChanged.connect(
        lambda: check_sele_timer.start(1000))

    form.button_ok.clicked.connect(run)

    # "Register" opens a web browser
    @form.button_register.clicked.connect
    def _():
        import webbrowser
        webbrowser.open("http://www.poissonboltzmann.org/")

    @form.button_load.clicked.connect
    def _():
        fnames = getOpenFileNames(
            None, filter='APBS Input (*.in);;All Files (*)')[0]
        if fnames:
            contents = load_apbs_in(form, fnames[0])
            set_apbs_in(contents)

    @form.button_reset.clicked.connect
    def _():
        set_apbs_in(electrostatics.template_apbs_in)

    form._dialog.show()
    form._dialog.resize(500, 600)


def load_apbs_in(form, filename, contents=''):
    import shlex
    from pymol import cmd, importing

    wdir = os.path.dirname(filename)

    if not contents:
        contents = cmd.file_read(filename)

    if not isinstance(contents, str):
        contents = contents.decode()

    sectionkeys = ('read', 'elec', 'apolar', 'print')
    section = ''
    insert_write_pot = True

    lines = []

    for line in contents.splitlines():
        a = shlex.split(line)
        key = a[0].lower() if a else ''

        if not section:
            if key in sectionkeys:
                section = key
        elif key == 'end':
            if section == 'elec' and insert_write_pot:
                lines.append('write pot dx "{mapfile}"')
            section = ''
        elif section == 'read':
            if len(a) > 2 and key in ('charge', 'kappa', 'mol', 'parm', 'pot'):
                filename = os.path.join(wdir, a[2])
                if os.path.exists(filename):
                    format = a[1].lower()
                    if key == 'mol' and format in ('pqr', 'pdb'):
                        # load into PyMOL and update selection dropdown
                        oname = importing.filename_to_objectname(a[2])
                        oname = cmd.get_unused_name(oname, 0)
                        cmd.load(filename, oname, format=format)
                        form.input_sele.addItem(oname)
                        form.input_sele.setEditText(oname)

                    # absolute path in input file
                    a[2] = '"' + filename + '"'
                    line = ' '.join(a)
                else:
                    QMessageBox.warning(
                        form._dialog, "Warning",
                        f'Warning: File "{filename}" does not exist')

        elif section == 'elec':
            if key == 'write':
                if a[1:4] == ['pot', 'dx', "{mapfile}"]:
                    insert_write_pot = False

        lines.append(line)

    return '\n'.join(lines)


def __init_plugin__(app=None):
    from pymol.plugins import addmenuitemqt as addmenuitem
    addmenuitem('APBS Electrostatics', dialog)
