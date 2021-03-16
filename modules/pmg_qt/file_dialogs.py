import os
import sys
import pymol

from pymol.Qt import QtGui, QtCore, QtWidgets
from pymol.Qt.utils import getSaveFileNameWithExt, AsyncFunc
from pymol.Qt.utils import PopupOnException

import urllib.request as urllib


def _get_cms_traj_file(fname):
    """For a .cms file, get the filename of the corresponding trajectory file.
    Return None if no such trajecotry file is found.
    """
    candidates = []

    # 2) By filename conventions
    stem = fname[:-8] if fname.endswith("-out.cms") else fname[:-4]
    candidates += [
        (stem + '_trj', 'clickme.dtr'),
        (stem + '.xtc',),
    ]

    for components in candidates:
        traj = os.path.join(*components)
        if os.path.exists(traj):
            return traj

    return None


def load_dialog(parent, fname, **kwargs):
    '''
    Load a file into PyMOL. May show a file format specific options
    dialog (e.g. trajectory loading dialog). Registers the filename
    in the "recent files" history.
    '''
    if '://' not in fname:
        parent.initialdir = os.path.dirname(fname)

    parent.recent_filenames_add(fname)

    format = pymol.importing.filename_to_format(fname)[2]

    if fname[-4:] in ['.dcd', '.dtr', '.xtc', '.trr']:
        load_traj_dialog(parent, fname)
    elif format in ('aln', 'fasta'):
        load_aln_dialog(parent, fname, format)
    elif format == 'mae':
        load_mae_dialog(parent, fname)
    elif format in ('ccp4', 'map'):
        load_map_dialog(parent, fname, 'ccp4')
    elif format == 'brix':
        load_map_dialog(parent, fname, 'o')
    elif format == 'mtz':
        load_mtz_dialog(parent, fname)
    else:
        if format in ('pse', 'psw') and not ask_partial(parent, kwargs, fname):
            return

        if format in ('pml', 'py', 'pym'):
            parent.cmd.cd(parent.initialdir, quiet=0)

        try:
            parent.cmd.load(fname, quiet=0, **kwargs)
        except BaseException as e:
            QtWidgets.QMessageBox.critical(parent, "Error", str(e))
            return

        # auto-load desmond trajectory
        if fname.endswith('.cms'):
            traj = _get_cms_traj_file(fname)
            if traj:
                load_traj_dialog(parent, traj)

    return True


def ask_partial(parent, kwargs, fname):
    if kwargs.get('partial', 0) or not parent.cmd.get_names():
        return True

    form = parent.load_form('askpartial')
    form.check_rename.setChecked(parent.cmd.get_setting_boolean(
        'auto_rename_duplicate_objects'))

    if not form._dialog.exec_():
        return False

    if form.check_partial.isChecked():
        kwargs['partial'] = 1
        parent.cmd.set('auto_rename_duplicate_objects',
                form.check_rename.isChecked(), quiet=0)
    elif form.check_new.isChecked():
        parent.new_window([fname])
        return False

    return True


def load_traj_dialog(parent, filename):
    '''Open a trajectory loading dialog'''
    names = parent.cmd.get_object_list()

    if not names:
        msg = "To load a trajectory, you first need to load a molecular object"  #noqa
        QtWidgets.QMessageBox.warning(parent, "Warning", msg)
        return

    form = parent.load_form('load_traj')
    form.input_object.addItems(names)
    form.input_object.setCurrentIndex(form.input_object.count() - 1)

    def get_command(*args):
        command = ''
        if form.input_dbm3.isChecked():
            command += 'set defer_builds_mode, 3\n'
        command += ('load_traj \\\n'
                   '    %s, \\\n'
                   '    %s, %d, \\\n'
                   '    start=%d, stop=%d, interval=%d' % (
                       filename,
                       form.input_object.currentText(),
                       form.input_state.value(),
                       form.input_start.value(),
                       form.input_stop.value(),
                       form.input_interval.value()))
        return command

    def update_output_command(*args):
        form.output_command.setText(get_command())

    def run():
        parent.cmd.do(get_command())
        form._dialog.close()

    # hook up events
    form.input_object.currentIndexChanged.connect(update_output_command)
    form.input_state.valueChanged.connect(update_output_command)
    form.input_start.valueChanged.connect(update_output_command)
    form.input_stop.valueChanged.connect(update_output_command)
    form.input_interval.valueChanged.connect(update_output_command)
    form.input_dbm3.toggled.connect(update_output_command)
    form.button_ok.clicked.connect(run)

    update_output_command()
    form._dialog.setModal(True)
    form._dialog.show()


def load_mtz_dialog(parent, filename):
    from pymol import headering

    _fileData = headering.MTZHeader(filename)

    FCols = _fileData.getColumnsOfType("F") + \
            _fileData.getColumnsOfType("G")
    PCols = _fileData.getColumnsOfType("P")
    WCols = _fileData.getColumnsOfType("W") + \
            _fileData.getColumnsOfType("Q")
    _2FC, _2PC, _looksLike = _fileData.guessCols("2FoFc")
    _FC, _PC, _looksLike = _fileData.guessCols("FoFc")

    form = parent.load_form('load_mtz')
    form.input_amplitudes.addItems(FCols)
    form.input_phases.addItems(PCols)
    form.input_weights.addItem("")
    form.input_weights.addItems(WCols)
    form.input_prefix.setFocus(),

    for col in [_2FC, _FC]:
        if col in FCols:
            form.input_amplitudes.setCurrentIndex(FCols.index(col))
            break
    for col in [_2PC, _PC]:
        if col in PCols:
            form.input_phases.setCurrentIndex(PCols.index(col))
            break

    if _fileData.reso_min is not None:
        form.input_reso_min.setValue(_fileData.reso_min)
    if _fileData.reso_max is not None:
        form.input_reso_max.setValue(_fileData.reso_max)

    def run():
        try:
            parent.cmd.load_mtz(filename,
                    form.input_prefix.text(),
                    form.input_amplitudes.currentText(),
                    form.input_phases.currentText(),
                    form.input_weights.currentText(),
                    form.input_reso_min.value(),
                    form.input_reso_max.value(),
                    quiet=0)
        except BaseException as e:
            QtWidgets.QMessageBox.critical(parent, "Error", str(e))

    form._dialog.accepted.connect(run)
    form._dialog.setModal(True)
    form._dialog.show()


def load_aln_dialog(parent, filename, format):
    _self = parent.cmd

    import numpy
    import difflib
    import pymol.seqalign as seqalign

    try:
        # for fasta format, this only succeeds if all sequences have the
        # same length, raises ValueError otherwise
        alignment = seqalign.aln_magic_read(filename)

        # a single sequence is not an aligment
        if format == "fasta" and len(alignment) < 2:
            raise ValueError
    except ValueError:
        # fasta files which don't contain alignments will be loaded as
        # extended structures (fab command) instead
        _self.load(filename)
        return

    # alignment record ids and PyMOL model names
    ids = [rec.id for rec in alignment]
    ids_remain = list(ids)
    models = _self.get_object_list()
    models_remain = list(models)
    mapping = {}

    N = len(ids)
    M = len(models)

    # ids -> models similarity matrix
    similarity = numpy.zeros((N, M))
    for i in range(N):
        for j in range(M):
            similarity[i, j] = difflib.SequenceMatcher(None,
                    ids[i], models[j], False).ratio()

    # guess mapping
    for _ in range(min(N, M)):
        i, j = numpy.unravel_index(similarity.argmax(), similarity.shape)
        mapping[ids_remain.pop(i)] = models_remain.pop(j)
        similarity = numpy.delete(similarity, i, axis=0)
        similarity = numpy.delete(similarity, j, axis=1)

    form = parent.load_form('load_aln')
    comboboxes = {}

    # mapping GUI
    for row, rec_id in enumerate(ids, 1):
        label = QtWidgets.QLabel(rec_id, form._dialog)
        combobox = QtWidgets.QComboBox(form._dialog)
        combobox.addItem("")
        combobox.addItems(models)
        combobox.setCurrentText(mapping.get(rec_id, ""))
        form.layout_mapping.addWidget(label, row, 0)
        form.layout_mapping.addWidget(combobox, row, 1)
        comboboxes[rec_id] = combobox

    def run():
        mapping = dict((rec_id, combobox.currentText())
                for (rec_id, combobox) in comboboxes.items())
        seqalign.load_aln_multi(filename, mapping=mapping, _self=_self)
        form._dialog.close()

    def cancel():
        form._dialog.close()
        if format == 'fasta' and QtWidgets.QMessageBox.question(
                parent, "Load as structures?",
                "Load sequences as extended structures instead?"
        ) == QtWidgets.QMessageBox.Yes:
            _self.load(filename)

    # hook up events
    form.button_ok.clicked.connect(run)
    form.button_cancel.clicked.connect(cancel)

    form._dialog.setModal(True)
    form._dialog.show()


def load_mae_dialog(parent, filename):
    form = parent.load_form('load_mae')

    form.input_object_name.setPlaceholderText(
            pymol.importing.filename_to_objectname(filename))
    form.input_object_props.setText(parent.cmd.get('load_object_props_default') or '*')
    form.input_atom_props.setText(parent.cmd.get('load_atom_props_default') or '*')

    def get_command(*args):
        command = ('load \\\n    %s' % (filename))
        name = form.input_object_name.text()
        if name:
            command += ', \\\n    ' + name
        command += ', \\\n    mimic=' + ('1' if form.input_mimic.isChecked() else '0')
        command += (
                ', \\\n    object_props=%s'
                ', \\\n    atom_props=%s' % (
                    form.input_object_props.text(),
                    form.input_atom_props.text()))
        multiplex, discrete = [
            (-2, -1),
            (0, 0),
            (0, 1),
            (1, -1),
        ][form.input_multiplex.currentIndex()]
        if discrete != -1:
            command += ', \\\n    discrete={}'.format(discrete)
        if multiplex != -2:
            command += ', \\\n    multiplex={}'.format(multiplex)
        return command

    def update_output_command(*args):
        form.output_command.setText(get_command())

    def run():
        parent.cmd.do(get_command())
        form._dialog.close()

    # hook up events
    form.input_multiplex.currentIndexChanged.connect(update_output_command)
    form.input_mimic.stateChanged.connect(update_output_command)
    form.input_object_name.textChanged.connect(update_output_command)
    form.input_object_props.textChanged.connect(update_output_command)
    form.input_atom_props.textChanged.connect(update_output_command)
    form.button_ok.clicked.connect(run)

    update_output_command()
    form._dialog.setModal(True)
    form._dialog.show()


def load_map_dialog(parent, filename, format='ccp4'):
    form = parent.load_form('load_map')
    normalize_setting = 'normalize_' + format + '_maps'

    form.input_object_name.setText(
            pymol.importing.filename_to_objectname(filename))
    form.input_normalize.setChecked(parent.cmd.get_setting_int(normalize_setting) > 0)

    def get_command(*args):
        command = 'set %s, %d\n' % (normalize_setting,
                1 if form.input_normalize.isChecked() else 0)
        command += 'load ' + filename
        name = form.input_object_name.text()
        if name:
            command += ', \\\n    ' + name
        else:
            name = pymol.importing.filename_to_objectname(filename)

        selesuffix = ''

        level = round(form.input_level.value(), 4)
        sele = form.input_selection.currentText()
        if sele:
            buf = form.input_buffer.value()
            selesuffix += ', %s, %s' % (sele, buf)
            if form.check_carve.isChecked():
                selesuffix += ', carve=%s' % (buf)

        if form.check_volume.isChecked():
            name_volume = form.input_name_volume.text() or (name + '_volume')
            command += '\nvolume %s, %s, %s blue .5 %s yellow 0' % (
                    name_volume, name, level, level * 2)
            command += selesuffix

        if form.check_isomesh.isChecked():
            name_isomesh = form.input_name_isomesh.text() or (name + '_isomesh')
            command += '\nisomesh %s, %s, %s' % (name_isomesh, name, level)
            command += selesuffix

        if form.check_isosurface.isChecked():
            name_isosurface = form.input_name_isosurface.text() or (name + '_isosurface')
            command += '\nisosurface %s, %s, %s' % (name_isosurface, name, level)
            command += selesuffix

        return command

    def update_output_command(*args):
        form.output_command.setText(get_command())

    def run():
        parent.cmd.do(get_command())
        form._dialog.close()

    # hook up events
    form.input_normalize.stateChanged.connect(update_output_command)
    form.input_object_name.textChanged.connect(update_output_command)
    form.check_volume.stateChanged.connect(update_output_command)
    form.check_isomesh.stateChanged.connect(update_output_command)
    form.check_isosurface.stateChanged.connect(update_output_command)
    form.check_carve.stateChanged.connect(update_output_command)
    form.input_name_volume.textChanged.connect(update_output_command)
    form.input_name_isomesh.textChanged.connect(update_output_command)
    form.input_name_isosurface.textChanged.connect(update_output_command)
    form.input_selection.editTextChanged.connect(update_output_command)
    form.input_level.valueChanged.connect(update_output_command)
    form.input_buffer.valueChanged.connect(update_output_command)
    form.button_ok.clicked.connect(run)

    update_output_command()
    form._dialog.setModal(True)
    form._dialog.show()


def _get_assemblies(pdbid):
    # TODO move to another module
    import json
    pdbid = pdbid.lower()
    url = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/" + pdbid
    try:
        data = json.load(urllib.urlopen(url))
        assembies = data[pdbid][0]['assemblies']
        return [a['assembly_id'] for a in assembies]
    except LookupError:
        pass
    except Exception as e:
        print('_get_assemblies failed')
        print(e)
    return []


def _get_chains(pdbid):
    # TODO move to another module
    import json
    pdbid = pdbid.lower()
    url = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/polymer_coverage/" + pdbid
    try:
        data = json.load(urllib.urlopen(url))
        return [
            chain['chain_id']
            for molecule in data[pdbid]['molecules']
            for chain in molecule['chains']
        ]
    except Exception as e:
        print('_get_chains failed')
        print(e)
    return []


def file_fetch_pdb(parent):
    form = parent.load_form('fetch')
    form.input_assembly.setEditText(parent.cmd.get('assembly'))

    def get_command(*args):
        code = form.input_code.text()
        if len(code) != 4:
            return ''

        def get_name(w):
            name = w.text()
            if name:
                return ', ' + name
            return ''

        if form.input_check_pdb.isChecked():
            command = 'set assembly, "%s"\nfetch %s%s%s' % (
                    form.input_assembly.currentText(),
                    code, form.input_chain.currentText(),
                    get_name(form.input_name))
        else:
            command = ''

        if form.input_check_2fofc.isChecked():
            command += '\nfetch %s%s, type=2fofc' % (
                    code, get_name(form.input_name_2fofc))

        if form.input_check_fofc.isChecked():
            command += '\nfetch %s%s, type=fofc' % (
                    code, get_name(form.input_name_fofc))

        return command

    def update_output_command(*args):
        form.output_command.setText(get_command())

    def code_changed(code):
        for combo in [form.input_assembly, form.input_chain]:
            if combo.count() != 1:
                text = combo.currentText()
                combo.clear()
                combo.addItem('')
                combo.setEditText(text)
        if len(code) == 4:
            update_assemblies(code)
            update_chains(code)
        update_output_command()

    def run():
        if len(form.input_code.text()) != 4:
            QtWidgets.QMessageBox.warning(parent, "Error", "Need 4 letter PDB code")
            return
        parent.cmd.do(get_command())
        form._dialog.close()

    # async events
    update_assemblies = AsyncFunc(_get_assemblies, form.input_assembly.addItems)
    update_chains = AsyncFunc(_get_chains, form.input_chain.addItems)

    # hook up events
    form.input_code.textChanged.connect(code_changed)
    form.input_chain.editTextChanged.connect(update_output_command)
    form.input_assembly.editTextChanged.connect(update_output_command)
    form.input_name.textChanged.connect(update_output_command)
    form.input_name_2fofc.textChanged.connect(update_output_command)
    form.input_name_fofc.textChanged.connect(update_output_command)
    form.input_check_pdb.stateChanged.connect(update_output_command)
    form.input_check_2fofc.stateChanged.connect(update_output_command)
    form.input_check_fofc.stateChanged.connect(update_output_command)
    form.button_ok.clicked.connect(run)

    update_output_command()
    form._dialog.show()


def file_save(parent):
    form = parent.load_form('save_molecule')
    default_selection = form.input_selection.currentText()

    get_setting_int = parent.cmd.get_setting_int
    form.input_no_pdb_conect_nodup.setChecked(  not get_setting_int('pdb_conect_nodup'))
    form.input_pdb_conect_all.setChecked(           get_setting_int('pdb_conect_all'))
    form.input_no_ignore_pdb_segi.setChecked(   not get_setting_int('ignore_pdb_segi'))
    form.input_pdb_retain_ids.setChecked(           get_setting_int('pdb_retain_ids'))
    form.input_retain_order.setChecked(             get_setting_int('retain_order'))

    models = parent.cmd.get_object_list()
    selections = parent.cmd.get_names('public_selections')
    names = models + selections

    form.input_state.addItems(map(str, range(1, parent.cmd.count_states() + 1)))

    form.input_selection.addItems(names)
    form.input_selection.lineEdit().setPlaceholderText(default_selection)

    formats = [
        'PDBx/mmCIF (*.cif *.cif.gz)',
        'PDB (*.pdb *.pdb.gz)',
        'PQR (*.pqr)',
        'MOL2 (*.mol2)',
        'MDL SD (*.sdf *.mol)',
        'Maestro (*.mae)',
        'MacroModel (*.mmd *.mmod *.dat)',
        'ChemPy Pickle (*.pkl)',
        'XYZ (*.xyz)',
        'MMTF (*.mmtf)',
        'By Extension (*.*)',
    ]

    @PopupOnException.decorator
    def run(*_):
        selection = form.input_selection.currentText() or default_selection
        state = int(form.input_state.currentText().split()[0])

        parent.cmd.set('pdb_conect_nodup',    not form.input_no_pdb_conect_nodup.isChecked())
        parent.cmd.set('pdb_conect_all',          form.input_pdb_conect_all.isChecked())
        parent.cmd.set('ignore_pdb_segi',     not form.input_no_ignore_pdb_segi.isChecked())
        parent.cmd.set('pdb_retain_ids',          form.input_pdb_retain_ids.isChecked())
        parent.cmd.set('retain_order',            form.input_retain_order.isChecked())

        if form.input_multi_state.isChecked():
            fmt = form.input_multi_state_fmt.text()
        elif form.input_multi_object.isChecked():
            fmt = form.input_multi_object_fmt.text()
        else:
            fmt = ''

        if fmt and form.input_multi_prompt.isChecked():
            fss = parent.cmd.multifilenamegen(fmt, selection, state)
        else:
            fss = [(fmt, selection, state)]

        for fname, selection, state in fss:
            fname = getSaveFileNameWithExt(parent,
                'Save Molecule As...',
                os.path.join(parent.initialdir, fname),
                filter=';;'.join(formats))

            if not fname:
                return

            parent.initialdir = os.path.dirname(fname)

            if form.input_multisave.isChecked():
                parent.cmd.multisave(fname, selection, state, quiet=0)
            elif '{' in os.path.basename(fname):
                parent.cmd.multifilesave(fname, selection, state, quiet=0)
            else:
                parent.cmd.save(fname, selection, state, quiet=0)
                parent.recent_filenames_add(fname)

        form._dialog.close()

    form.input_multi_state.pressed.connect(
        lambda: form.input_state.setCurrentIndex(1))

    form.button_ok.clicked.connect(run)
    form._dialog.show()


def file_save_png(parent):  #noqa
    if parent.dialog_png is not None:
        parent.dialog_png.show()
        return

    form = parent.load_form('png')
    parent.dialog_png = form._dialog

    def run():
        from pymol import exporting
        fname = getSaveFileNameWithExt(parent, 'Save As...', parent.initialdir,
                                filter='PNG File (*.png)')
        if not fname:
            return

        parent.initialdir = os.path.dirname(fname)

        rendering = form.input_rendering.currentIndex()
        ray = 0
        width, height, dpi = 0, 0, -1

        '''
        dpi = float(form.input_dpi.currentText())

        width = exporting._unit2px(
            form.input_width.value(), dpi,
            form.input_width_unit.currentText())
        height = exporting._unit2px(
            form.input_height.value(), dpi,
            form.input_height_unit.currentText())
        '''

        form._dialog.hide()

        if rendering == 1:
            parent.cmd.do('draw %d, %d' % (width, height))
            width = 0
            height = 0
        elif rendering == 2:
            parent.cmd.do('set opaque_background, 1')
            ray = 1
        elif rendering == 3:
            parent.cmd.do('set opaque_background, 0')
            ray = 1

        parent.cmd.sync()
        parent.cmd.do(
            'png %s, %d, %d, %d, ray=%d' %
            (fname, width, height, dpi, ray))

    '''
    def units_changed():
        dpi = float(form.input_dpi.currentText())
        width_unit = form.input_width_unit.currentText()
        height_unit = form.input_width_unit.currentText()
        if dpi < 1 and (width_unit != 'px' or height_unit != 'px'):
            form.input_dpi.setCurrentIndex(1)
            dpi = float(form.input_dpi.currentText())

    # initial values
    form.input_width.setValue(pymol.cmd.get_viewport()[0])

    dpi_index = 0
    dpi_values = [-1, 150, 300]
    dpi = pymol.cmd.get_setting_int('image_dots_per_inch')

    if dpi > 0:
        try:
            dpi_index = dpi_values.index(dpi)
        except ValueError:
            dpi_values.append(dpi)
            dpi_values.sort()
            dpi_index = dpi_values.index(dpi)

    for dpi in dpi_values:
        form.input_dpi.addItem(str(dpi))
    form.input_dpi.setCurrentIndex(dpi_index)

    # hook up events
    form.input_width_unit.currentIndexChanged.connect(units_changed)
    form.input_height_unit.currentIndexChanged.connect(units_changed)
    '''
    form.button_ok.clicked.connect(run)

    form._dialog.show()


def file_save_mpeg(parent, _preselect=None):
    form = parent.load_form('movieexport')

    filters = {
        'png': 'Numbered PNG Files (*.png)',
        'mp4': 'MPEG 4 movie file (*.mp4)',
        'mpg': 'MPEG 1 movie file (*.mpg *.mpeg)',
        'mov': 'QuickTime (*.mov)',
        'gif': 'Animated GIF (*.gif)',
    }

    support = {
        '':             {'mp4': 0, 'mpg': 0, 'mov': 0, 'gif': 0},
        'ffmpeg':       {'mp4': 1, 'mpg': 1, 'mov': 1, 'gif': 1},
        'mpeg_encode':  {'mp4': 0, 'mpg': 1, 'mov': 0, 'gif': 0},
        'convert':      {'mp4': 0, 'mpg': 0, 'mov': 0, 'gif': 1},
    }

    from pymol.movie import find_exe as has_exe

    def update_encoder_options():
        encoder = form.input_encoder.currentText()

        for fmt, enabled in support[encoder].items():
            w = getattr(form, 'format_' + fmt)
            w.setEnabled(enabled)
            if not enabled and w.isChecked():
                if encoder == 'mpeg_encode':
                    form.format_mpg.setChecked(True)
                elif encoder == 'convert':
                    form.format_gif.setChecked(True)
                else:
                    form.format_png.setChecked(True)

        form.input_quality.setEnabled(encoder not in ("", "convert"))

        if encoder and not has_exe(encoder):
            msg = "Encoder '%s' is not installed." % encoder
            pkg = None
            url = None
            if not pkg:
                QtWidgets.QMessageBox.warning(parent, "Warning", msg)
            else:
                from pymol.Qt import utils
                utils.conda_ask_install(pkg[1], pkg[0], msg, url=url)

    if _preselect == 'png':
        form.format_png.setChecked(True)
        form.group_format.hide()
    else:
        for i in range(1, form.input_encoder.count()):
            encoder = form.input_encoder.itemText(i)
            if has_exe(encoder):
                form.input_encoder.setCurrentIndex(i)
                break

        if _preselect == 'mov':
            encoder = 'ffmpeg'
            form.input_encoder.setCurrentIndex(1)
            form.format_mov.setChecked(True)
        elif encoder == 'ffmpeg':
            form.format_mp4.setChecked(True)
        else:
            form.format_mpg.setChecked(True)

    form._dialog.adjustSize()

    @PopupOnException.decorator
    def run(*_):
        for fmt in filters:
            w = getattr(form, 'format_' + fmt)
            if w.isChecked():
                break
        fname = getSaveFileNameWithExt(parent,
                'Save As...', parent.initialdir,
                filter=filters[fmt])
        if not fname:
            return

        parent.initialdir = os.path.dirname(fname)

        if fmt == 'png':
            parent.cmd.mpng(fname,
                    width=form.input_width.value(),
                    height=form.input_height.value(),
                    mode=2 if form.input_ray.isChecked() else 1,
                    quiet=0, modal=-1)
        else:
            mode = 'ray' if form.input_ray.isChecked() else 'draw'
            encoder = form.input_encoder.currentText()

            parent.cmd.movie.produce(fname,
                    width=form.input_width.value(),
                    height=form.input_height.value(),
                    quality=form.input_quality.value(),
                    mode=mode, encoder=encoder,
                    quiet=0)

    def set_resolution(height):
        w, h = form.input_width.value(), max(1, form.input_height.value())
        aspect = (w / h) if (w and h) else 9999
        width = int(round(min(aspect, 16. / 9.) * height / 2.)) * 2
        form.input_width.setValue(width)
        form.input_height.setValue(height)

    # initial values
    viewport = parent.cmd.get_viewport()
    form.input_width.setValue(viewport[0])
    form.input_height.setValue(viewport[1])
    form.input_quality.setValue(parent.cmd.get_setting_int('movie_quality'))
    if parent.cmd.get_setting_int('ray_trace_frames'):
        form.input_ray.setChecked(True)
    update_encoder_options()

    # hook up events
    form.button_ok.clicked.connect(run)
    form.input_encoder.currentIndexChanged.connect(update_encoder_options)

    for height in (720, 480, 360):
        button = getattr(form, 'button_%dp' % height)
        button.pressed.connect(lambda h=height: set_resolution(h))

    form._dialog.show()


def _file_save_object(self, otype, formats, noobjectsmsg):
    names = self.cmd.get_names_of_type(otype)

    if not names:
        QtWidgets.QMessageBox.warning(self, "Warning", noobjectsmsg)
        return

    form = self.load_form('save_object')
    form.input_name.addItems(names)
    form._dialog.setWindowTitle('Save ' + otype)

    def run():
        name = form.input_name.currentText()

        fname = getSaveFileNameWithExt(self,
            'Save As...',
            self.initialdir,
            filter=';;'.join(formats))

        if not fname:
            return

        self.cmd.save(fname, name, -1, quiet=0)
        form._dialog.close()

    form.button_ok.clicked.connect(run)
    form._dialog.show()


def file_save_map(self):
    return _file_save_object(self, 'object:map', ['CCP4 (*.ccp4 *.map)'],
            'No map objects loaded')


def file_save_aln(self):
    url = "http://pymolwiki.org/index.php/Align#Alignment_Objects"
    return _file_save_object(self, 'object:alignment', ['clustalw (*.aln)'],
            'No alignment objects loaded\n\n'
            'Hint: create alignment objects with "align" and\n'
            '"super" using the "object=..." argument.')
