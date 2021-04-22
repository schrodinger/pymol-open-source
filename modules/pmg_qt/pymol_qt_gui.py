"""
Contains main class for PyMOL QT GUI
"""


from collections import defaultdict
import os
import re
import sys

import pymol
import pymol._gui
from pymol import colorprinting, save_shortcut

from pymol.Qt import QtGui, QtCore, QtWidgets
from pymol.Qt.utils import (getSaveFileNameWithExt, UpdateLock, WidgetMenu,
        MainThreadCaller,
        PopupOnException,
        connectFontContextMenu, getMonospaceFont)

from .pymol_gl_widget import PyMOLGLWidget
from . import keymapping

from pmg_qt import properties_dialog, file_dialogs

Qt = QtCore.Qt
QFileDialog = QtWidgets.QFileDialog
getOpenFileNames = QFileDialog.getOpenFileNames


class PyMOLQtGUI(QtWidgets.QMainWindow, pymol._gui.PyMOLDesktopGUI):
    '''
    PyMOL QMainWindow GUI
    '''

    from pmg_qt.file_dialogs import (
            load_dialog,
            load_mae_dialog,
            file_fetch_pdb,
            file_save_png,
            file_save_mpeg,
            file_save_map,
            file_save_aln,
            file_save
    )

    _ext_window_visible = True
    _initialdir = ''

    def keyPressEvent(self, ev):
        args = keymapping.keyPressEventToPyMOLButtonArgs(ev)

        if args is not None:
            self.pymolwidget.pymol.button(*args)

    def closeEvent(self, event):
        self.cmd.quit()

    # for thread-safe viewport command
    viewportsignal = QtCore.Signal(int, int)

    def pymolviewport(self, w, h):
        cw, ch = self.cmd.get_viewport()
        pw = self.pymolwidget
        scale = pw.fb_scale

        # maintain aspect ratio
        if h < 1:
            if w < 1:
                pw.pymol.reshape(int(scale * pw.width()),
                                 int(scale * pw.height()), True)
                return
            h = (w * ch) / cw
        if w < 1:
            w = (h * cw) / ch

        win_size = self.size()
        delta = QtCore.QSize(w - cw, h - ch) / scale

        # window resize
        self.resize(delta + win_size)

    def get_view(self):
        self.cmd.get_view(2, quiet=0)
        QtWidgets.QApplication.clipboard().setText(self.cmd.get_view(3))
        print(" get_view: matrix copied to clipboard.")

    def __init__(self):  # noqa
        QtWidgets.QMainWindow.__init__(self)
        self.setDockOptions(QtWidgets.QMainWindow.AllowTabbedDocks |
                            QtWidgets.QMainWindow.AllowNestedDocks)

        # resize Window before it is shown
        options = pymol.invocation.options
        self.resize(
            options.win_x + (220 if options.internal_gui else 0),
            options.win_y + (246 if options.external_gui else 18))

        # for thread-safe viewport command
        self.viewportsignal.connect(self.pymolviewport)

        # reusable dialogs
        self.dialog_png = None
        self.advanced_settings_dialog = None
        self.props_panel = None
        self.builder = None
        self.shortcut_menu_filter_dialog = None

        # setting index -> callable
        self.setting_callbacks = defaultdict(list)

        # "session_file" setting in window title
        self.setting_callbacks[440].append(
            lambda v: self.setWindowTitle("PyMOL (" + os.path.basename(v) + ")")
        )

        # "External" Command Line and Loggin Widget
        self._setup_history()
        self.lineedit = CommandLineEdit()
        self.lineedit.setObjectName("command_line")
        self.browser = QtWidgets.QPlainTextEdit()
        self.browser.setObjectName("feedback_browser")
        self.browser.setReadOnly(True)

        # convenience: clicking into feedback browser gives focus to command
        # line. Drawback: Copying with CTRL+C doesn't work in feedback
        # browser -> clear focus proxy while text selected
        self.browser.setFocusProxy(self.lineedit)

        @self.browser.copyAvailable.connect
        def _(yes):
            self.browser.setFocusProxy(None if yes else self.lineedit)
            self.browser.setFocus()

        # Font
        self.browser.setFont(getMonospaceFont())
        connectFontContextMenu(self.browser)

        lineeditlayout = QtWidgets.QHBoxLayout()
        command_label = QtWidgets.QLabel("PyMOL>")
        command_label.setObjectName("command_label")
        lineeditlayout.addWidget(command_label)
        lineeditlayout.addWidget(self.lineedit)
        self.lineedit.setToolTip('''Command Input Area

Get the list of commands by hitting <TAB>

Get the list of arguments for one command with a question mark:
PyMOL> color ?

Read the online help for a command with "help":
PyMOL> help color

Get autocompletion for many arguments by hitting <TAB>
PyMOL> color ye<TAB>    (will autocomplete "yellow")
''')

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.browser)
        layout.addLayout(lineeditlayout)

        quickbuttonslayout = QtWidgets.QVBoxLayout()
        quickbuttonslayout.setSpacing(2)

        extguilayout = QtWidgets.QBoxLayout(QtWidgets.QBoxLayout.LeftToRight)
        extguilayout.setContentsMargins(2, 2, 2, 2)
        extguilayout.addLayout(layout)
        extguilayout.addLayout(quickbuttonslayout)

        class ExtGuiFrame(QtWidgets.QFrame):
            def mouseDoubleClickEvent(_, event):
                self.toggle_ext_window_dockable(True)

            _size_hint = QtCore.QSize(options.win_x, options.ext_y)

            def sizeHint(self):
                return self._size_hint

        dockWidgetContents = ExtGuiFrame(self)
        dockWidgetContents.setLayout(extguilayout)
        dockWidgetContents.setObjectName("extgui")

        self.ext_window = \
            dockWidget = QtWidgets.QDockWidget(self)
        dockWidget.setWindowTitle("External GUI")
        dockWidget.setWidget(dockWidgetContents)
        if options.external_gui:
            dockWidget.setTitleBarWidget(QtWidgets.QWidget())
        else:
            dockWidget.hide()

        self.addDockWidget(Qt.TopDockWidgetArea, dockWidget)

        # rearrange vertically if docking left or right
        @dockWidget.dockLocationChanged.connect
        def _(area):
            if area == Qt.LeftDockWidgetArea or area == Qt.RightDockWidgetArea:
                extguilayout.setDirection(QtWidgets.QBoxLayout.BottomToTop)
                quickbuttonslayout.takeAt(quickbuttons_stretch_index)
            else:
                extguilayout.setDirection(QtWidgets.QBoxLayout.LeftToRight)
                if quickbuttons_stretch_index >= quickbuttonslayout.count():
                    quickbuttonslayout.addStretch()

        # OpenGL Widget
        self.pymolwidget = PyMOLGLWidget(self)
        self.setCentralWidget(self.pymolwidget)

        cmd = self.cmd = self.pymolwidget.cmd

        '''
        # command completion
        completer = QtWidgets.QCompleter(cmd.kwhash.keywords, self)
        self.lineedit.setCompleter(completer)
        '''

        # overload <Tab> action
        self.lineedit.installEventFilter(self)
        self.pymolwidget.installEventFilter(self)

        # Quick Buttons
        for row in [
            [
                ('Reset', cmd.reset),
                ('Zoom', lambda: cmd.zoom(animate=1.0)),
                ('Orient', lambda: cmd.orient(animate=1.0)),

                # render dialog will be constructed when the menu is shown
                # for the first time. This way it's populated with the current
                # viewport and settings. Also defers parsing of the ui file.
                ('Draw/Ray', WidgetMenu(self).setSetupUi(self.render_dialog)),
            ],
            [
                ('Unpick', cmd.unpick),
                ('Deselect', cmd.deselect),
                ('Rock', cmd.rock),
                ('Get View', self.get_view),
            ],
            [
                ('|<', cmd.rewind),
                ('<', cmd.backward),
                ('Stop', cmd.mstop),
                ('Play', cmd.mplay),
                ('>', cmd.forward),
                ('>|', cmd.ending),
                ('MClear', cmd.mclear),
            ],
            [
                ('Builder', self.open_builder_panel),
                ('Properties', self.open_props_dialog),
                ('Rebuild', cmd.rebuild),
            ],
        ]:
            hbox = QtWidgets.QHBoxLayout()
            hbox.setSpacing(2)

            for name, callback in row:
                btn = QtWidgets.QPushButton(name)
                btn.setProperty("quickbutton", True)
                btn.setAttribute(Qt.WA_LayoutUsesWidgetRect) # OS X workaround
                hbox.addWidget(btn)

                if callback is None:
                    btn.setEnabled(False)
                elif isinstance(callback, QtWidgets.QMenu):
                    btn.setMenu(callback)
                else:
                    btn.released.connect(callback)

            quickbuttonslayout.addLayout(hbox)

        # progress bar
        hbox = QtWidgets.QHBoxLayout()
        self.progressbar = QtWidgets.QProgressBar()
        self.progressbar.setSizePolicy(
                QtWidgets.QSizePolicy.Minimum,
                QtWidgets.QSizePolicy.Minimum)
        hbox.addWidget(self.progressbar)
        self.abortbutton = QtWidgets.QPushButton('Abort')
        self.abortbutton.setStyleSheet("background: #FF0000; color: #FFFFFF")
        self.abortbutton.released.connect(cmd.interrupt)
        hbox.addWidget(self.abortbutton)
        quickbuttonslayout.addLayout(hbox)

        quickbuttonslayout.addStretch()
        quickbuttons_stretch_index = quickbuttonslayout.count() - 1

        # menu top level
        self.menubar = menubar = self.menuBar()

        # action groups
        actiongroups = {}

        def _addmenu(data, menu):
            '''Fill a menu from "data"'''
            menu.setTearOffEnabled(True)
            menu.setWindowTitle(menu.title())  # needed for Windows
            for item in data:
                if item[0] == 'separator':
                    menu.addSeparator()
                elif item[0] == 'menu':
                    _addmenu(item[2], menu.addMenu(item[1].replace('&', '&&')))
                elif item[0] == 'command':
                    command = item[2]
                    if command is None:
                        print('warning: skipping', item)
                    else:
                        if isinstance(command, str):
                            command = lambda c=command: cmd.do(c)
                        menu.addAction(item[1], command)
                elif item[0] == 'check':
                    if len(item) > 4:
                        menu.addAction(
                            SettingAction(self, cmd, item[2], item[1],
                                          item[3], item[4]))
                    else:
                        menu.addAction(
                            SettingAction(self, cmd, item[2], item[1]))
                elif item[0] == 'radio':
                    label, name, value = item[1:4]
                    try:
                        group, type_, values = actiongroups[item[2]]
                    except KeyError:
                        group = QtWidgets.QActionGroup(self)
                        type_, values = cmd.get_setting_tuple(name)
                        actiongroups[item[2]] = group, type_, values
                    action = QtWidgets.QAction(label, self)
                    action.triggered.connect(lambda _=0, args=(name, value):
                                             cmd.set(*args, log=1, quiet=0))

                    self.setting_callbacks[cmd.setting._get_index(
                        name)].append(
                            lambda v, V=value, a=action: a.setChecked(v == V))

                    group.addAction(action)
                    menu.addAction(action)
                    action.setCheckable(True)
                    if values[0] == value:
                        action.setChecked(True)
                elif item[0] == 'open_recent_menu':
                    self.open_recent_menu = menu.addMenu('Open Recent...')
                else:
                    print('error:', item)

        # recent files menu
        self.open_recent_menu = None

        # for plugins
        self.menudict = {'': menubar}

        # menu
        for _, label, data in self.get_menudata(cmd):
            assert _ == 'menu'
            menu = menubar.addMenu(label)
            self.menudict[label] = menu
            _addmenu(data, menu)

        # hack for macOS to hide "Edit > Start Dictation"
        # https://bugreports.qt.io/browse/QTBUG-43217
        if pymol.IS_MACOS:
            self.menudict['Edit'].setTitle('Edit_')
            QtCore.QTimer.singleShot(10, lambda:
                    self.menudict['Edit'].setTitle('Edit'))

        # recent files menu
        if self.open_recent_menu:
            @self.open_recent_menu.aboutToShow.connect
            def _():
                self.open_recent_menu.clear()
                for fname in self.recent_filenames:
                    self.open_recent_menu.addAction(
                            fname if len(fname) < 128 else '...' + fname[-120:],
                            lambda fname=fname: self.load_dialog(fname))

        # some experimental window control
        menu = self.menudict['Display'].addSeparator()
        menu = self.menudict['Display'].addMenu('External GUI')
        menu.addAction('Toggle floating', self.toggle_ext_window_dockable,
                       QtGui.QKeySequence('Ctrl+E'))
        ext_vis_action = self.ext_window.toggleViewAction()
        ext_vis_action.setText('Visible')
        menu.addAction(ext_vis_action)

        # extra key mappings (MacPyMOL compatible)
        QtWidgets.QShortcut(QtGui.QKeySequence('Ctrl+O'), self).activated.connect(self.file_open)
        QtWidgets.QShortcut(QtGui.QKeySequence('Ctrl+S'), self).activated.connect(self.session_save)

        # feedback
        self.feedback_timer = QtCore.QTimer()
        self.feedback_timer.setSingleShot(True)
        self.feedback_timer.timeout.connect(self.update_feedback)
        self.feedback_timer.start(100)

        # legacy plugin system
        self.menudict['Plugin'].addAction(
            'Initialize Plugin System', self.initializePlugins)

        # focus in command line
        if options.external_gui:
            self.lineedit.setFocus()
        else:
            self.pymolwidget.setFocus()

        # Apply PyMOL stylesheet
        try:
            with open(cmd.exp_path('$PYMOL_DATA/pmg_qt/styles/pymol.sty')) as f:
                style = f.read()
        except IOError:
            print('Could not read PyMOL stylesheet.')
            print('DEBUG: PYMOL_DATA=' + repr(os.getenv('PYMOL_DATA')))
            style = ""

        if style:
            self.setStyleSheet(style)

        # Load saved shortcuts on launch
        self.saved_shortcuts = pymol.save_shortcut.load_and_set(self.cmd)

    def lineeditKeyPressEventFilter(self, watched, event):
        key = event.key()
        if key == Qt.Key_Tab:
            self.complete()
        elif key == Qt.Key_Up:
            if event.modifiers() & Qt.ControlModifier:
                self.back_search()
            else:
                self.back()
        elif key == Qt.Key_Down:
            self.forward()
        elif key == Qt.Key_Return or key == Qt.Key_Enter:
            # filter out "Return" instead of binding lineedit.returnPressed,
            # because otherwise OrthoKey would capture it as well.
            self.doPrompt()
        else:
            return False
        return True

    def eventFilter(self, watched, event):
        '''
        Filter out <Tab> event to do tab-completion instead of move focus
        '''
        type_ = event.type()
        if type_ == QtCore.QEvent.KeyRelease:
            if event.key() == Qt.Key_Tab:
                # silently skip tab release
                return True
        elif type_ == QtCore.QEvent.KeyPress:
            if watched is self.lineedit:
                return self.lineeditKeyPressEventFilter(watched, event)
            elif event.key() == Qt.Key_Tab:
                self.keyPressEvent(event)
                return True
        return False

    def toggle_ext_window_dockable(self, neverfloat=False):
        '''
        Toggle whether the "external" GUI is dockable
        '''
        dockWidget = self.ext_window

        if dockWidget.titleBarWidget() is None:
            tbw = QtWidgets.QWidget()
        else:
            tbw = None

        dockWidget.setFloating(tbw is None and not neverfloat)
        dockWidget.setTitleBarWidget(tbw)
        dockWidget.show()

    def toggle_fullscreen(self, toggle=-1):
        '''
        Full screen
        '''
        is_fullscreen = self.windowState() == Qt.WindowFullScreen

        if toggle == -1:
            toggle = not is_fullscreen

        if not is_fullscreen:
            self._ext_window_visible = self.ext_window.isVisible()

        if toggle:
            self.menubar.hide()
            if not self.ext_window.isFloating():
                self.ext_window.hide()
            self.showFullScreen()
            self.pymolwidget.setFocus()
        else:
            self.menubar.show()
            if self._ext_window_visible:
                self.ext_window.show()
            self.showNormal()

    @property
    def initialdir(self):
        '''
        Be in sync with cd/pwd on the console until the first file has been
        browsed, then remember the last directory.
        '''
        return self._initialdir or os.getcwd()

    @initialdir.setter
    def initialdir(self, value):
        self._initialdir = value

    ##################
    # UI Forms
    ##################

    def load_form(self, name, dialog=None):
        '''Load a form from pmg_qt/forms/{name}.py'''
        import importlib
        if dialog is None:
            dialog = QtWidgets.QDialog(self)
            widget = dialog
        elif dialog == 'floating':
            widget = QtWidgets.QWidget(self)
        else:
            widget = dialog

        try:
            m = importlib.import_module('.forms.' + name, 'pmg_qt')
        except ImportError as e:
            if pymol.Qt.DEBUG:
                print('load_form import failed (%s)' % (e,))
            uifile = os.path.join(os.path.dirname(__file__), 'forms', '%s.ui' % name)
            form = pymol.Qt.utils.loadUi(uifile, widget)
        else:
            if hasattr(m, 'Ui_Form'):
                form = m.Ui_Form()
            else:
                form = m.Ui_Dialog()

            form.setupUi(widget)

        if dialog == 'floating':
            dialog = QtWidgets.QDockWidget(widget.windowTitle(), self)
            dialog.setFloating(True)
            dialog.setWidget(widget)
            dialog.resize(widget.size())

        form._dialog = dialog
        return form

    def edit_colors_dialog(self):
        form = self.load_form('colors')
        form.list_colors.setSortingEnabled(True)

        # populate list with named colors
        for color_index in self.cmd.get_color_indices():
            form.list_colors.addItem(color_index[0])

        # update spinboxes for given color
        def load_color(name):
            index = self.cmd.get_color_index(name)
            if index == -1:
                return
            rgb = self.cmd.get_color_tuple(index)
            form.input_R.setValue(rgb[0])
            form.input_G.setValue(rgb[1])
            form.input_B.setValue(rgb[2])

        # update spinbox from slider
        spinbox_lock = [False]
        def update_spinbox(spinbox, value):
            if not spinbox_lock[0]:
                spinbox.setValue(value / 100.)

        # update sliders and colored frame
        def update_gui(*args):
            spinbox_lock[0] = True
            R = form.input_R.value()
            G = form.input_G.value()
            B = form.input_B.value()
            form.slider_R.setValue(R * 100)
            form.slider_G.setValue(G * 100)
            form.slider_B.setValue(B * 100)
            form.frame_color.setStyleSheet(
                "background-color: rgb(%d,%d,%d)" % (
                    R * 0xFF, G * 0xFF, B * 0xFF))
            spinbox_lock[0] = False

        def run():
            name  = form.input_name.text()
            R = form.input_R.value()
            G = form.input_G.value()
            B = form.input_B.value()

            self.cmd.do('set_color %s, [%.2f, %.2f, %.2f]\nrecolor' %
                        (name, R, G, B))

            # if new color, insert and make current row
            if not form.list_colors.findItems(name, Qt.MatchExactly):
                form.list_colors.addItem(name)
                form.list_colors.setCurrentItem(
                    form.list_colors.findItems(name, Qt.MatchExactly)[0])

        # hook up events
        form.slider_R.valueChanged.connect(lambda v: update_spinbox(form.input_R, v))
        form.slider_G.valueChanged.connect(lambda v: update_spinbox(form.input_G, v))
        form.slider_B.valueChanged.connect(lambda v: update_spinbox(form.input_B, v))
        form.input_R.valueChanged.connect(update_gui)
        form.input_G.valueChanged.connect(update_gui)
        form.input_B.valueChanged.connect(update_gui)
        form.input_name.textChanged.connect(load_color)
        form.list_colors.currentTextChanged.connect(form.input_name.setText)
        form.button_apply.clicked.connect(run)

        form._dialog.show()

    def open_builder_panel(self):
        from pmg_qt.builder import BuilderPanelDocked
        from pymol import plugins

        app = plugins.get_pmgapp()
        if not self.builder:
            self.builder = BuilderPanelDocked(self, app)
            self.addDockWidget(Qt.TopDockWidgetArea, self.builder)

        self.builder.show()
        self.builder.raise_()

    def open_props_dialog(self):
        from .properties_dialog import PropsDialog

        if not self.props_panel:
            self.props_panel = PropsDialog(self)

        self.props_panel.get_dialog().show()
        self.props_panel.get_dialog().raise_()

    def edit_pymolrc(self):
        from . import TextEditor
        from pymol import plugins
        TextEditor.edit_pymolrc(plugins.get_pmgapp())

    ##################
    # Menu callbacks
    ##################

    def file_open(self):
        fnames = getOpenFileNames(self, 'Open file', self.initialdir)[0]
        partial = 0
        for fname in fnames:
            if not self.load_dialog(fname, partial=partial):
                break
            partial = 1

    def session_save(self):
        fname = self.cmd.get('session_file')
        fname = self.cmd.as_pathstr(fname)
        return self.session_save_as(fname)

    @PopupOnException.decorator
    def session_save_as(self, fname=''):
        formats = [
            'PyMOL Session File (*.pse *.pze *.pse.gz)',
            'PyMOL Show File (*.psw *.pzw *.psw.gz)',
        ]
        if not fname:
            fname = getSaveFileNameWithExt(
                self,
                'Save Session As...',
                self.initialdir,
                filter=';;'.join(formats))
        if fname:
            self.initialdir = os.path.dirname(fname)
            self.cmd.save(fname, format='pse', quiet=0)
            self.recent_filenames_add(fname)

    def render_dialog(self, widget=None):
        form = self.load_form('render', widget)
        lock = UpdateLock([ZeroDivisionError])

        def get_factor():
            units = form.input_units.currentText()
            factor = 1.0 if units == 'inch' else 2.54
            return factor / float(form.input_dpi.currentText())

        @lock.skipIfCircular
        def update_units(*args):
            width = form.input_width.value()
            height = form.input_height.value()
            factor = get_factor()
            form.input_width_units.setValue(width * factor)
            form.input_height_units.setValue(height * factor)

        @lock.skipIfCircular
        def update_pixels(*args):
            width = form.input_width_units.value()
            height = form.input_height_units.value()
            factor = get_factor()
            form.input_width.setValue(width / factor)
            form.input_height.setValue(height / factor)

        @lock.skipIfCircular
        def update_width(*args):
            if form.aspectratio > 0:
                width = form.input_height.value() * form.aspectratio
                form.input_width.setValue(int(width))
                form.input_width_units.setValue(width * get_factor())

        @lock.skipIfCircular
        def update_height(*args):
            if form.aspectratio > 0:
                height = form.input_width.value() / form.aspectratio
                form.input_height.setValue(int(height))
                form.input_height_units.setValue(height * get_factor())

        def update_aspectratio(checked=True):
            if checked:
                try:
                    form.aspectratio = (
                            float(form.input_width.value()) /
                            float(form.input_height.value()))
                except ZeroDivisionError:
                    form.button_lock.setChecked(False)
            else:
                form.aspectratio = 0

        def update_from_viewport():
            w, h = self.cmd.get_viewport()
            form.aspectratio = 0
            form.input_width.setValue(w)
            form.input_height.setValue(h)
            update_aspectratio(form.button_lock.isChecked())

        def run_draw(ray=False):
            width = form.input_width.value()
            height = form.input_height.value()
            if ray:
                self.cmd.set('opaque_background',
                        not form.input_transparent.isChecked())
                self.cmd.do('ray %d, %d, async=1' % (width, height))
            else:
                self.cmd.do('draw %d, %d' % (width, height))
            form.stack.setCurrentIndex(1)

        def run_ray():
            run_draw(ray=True)

        def run_save():
            fname = getSaveFileNameWithExt(self, 'Save As...', self.initialdir,
                    filter='PNG File (*.png)')
            if not fname:
                return
            self.initialdir = os.path.dirname(fname)
            self.cmd.png(fname, prior=1, dpi=form.input_dpi.currentText())

        def run_copy_clipboard():
            with PopupOnException():
                _copy_image(self.cmd, False, form.input_dpi.currentText())

        dpi = self.cmd.get_setting_int('image_dots_per_inch')
        if dpi > 0:
            form.input_dpi.setEditText(str(dpi))
        form.input_dpi.setValidator(QtGui.QIntValidator())

        form.input_units.currentIndexChanged.connect(update_units)
        form.input_dpi.editTextChanged.connect(update_pixels)
        form.input_width.valueChanged.connect(update_units)
        form.input_height.valueChanged.connect(update_units)
        form.input_width_units.valueChanged.connect(update_pixels)
        form.input_height_units.valueChanged.connect(update_pixels)

        # set values before connecting mutual width<->height updates
        update_from_viewport()

        form.input_width.valueChanged.connect(update_height)
        form.input_height.valueChanged.connect(update_width)
        form.input_width_units.valueChanged.connect(update_height)
        form.input_height_units.valueChanged.connect(update_width)
        form.button_lock.toggled.connect(update_aspectratio)

        form.button_draw.clicked.connect(run_draw)
        form.button_ray.clicked.connect(run_ray)
        form.button_current.clicked.connect(update_from_viewport)
        form.button_back.clicked.connect(lambda: form.stack.setCurrentIndex(0))
        form.button_clip.clicked.connect(run_copy_clipboard)
        form.button_save.clicked.connect(run_save)

        if widget is None:
            form._dialog.show()

    @PopupOnException.decorator
    def _file_save(self, filter, format):
        fname = getSaveFileNameWithExt(
            self,
            'Save As...',
            self.initialdir,
            filter=filter)
        if fname:
            self.cmd.save(fname, format=format, quiet=0)

    def file_save_wrl(self):
        self._file_save('VRML 2 WRL File (*.wrl)', 'wrl')

    def file_save_dae(self):
        self._file_save('COLLADA File (*.dae)', 'dae')

    def file_save_pov(self):
        self._file_save('POV File (*.pov)', 'pov')

    def file_save_mpng(self):
        self.file_save_mpeg('png')

    def file_save_mov(self):
        self.file_save_mpeg('mov')

    def file_save_stl(self):
        self._file_save('STL File (*.stl)', 'stl')

    def file_save_gltf(self):
        self._file_save('GLTF File (*.gltf)', 'gltf')

    LOG_FORMATS = [
        'PyMOL Script (*.pml)',
        'Python Script (*.py *.pym)',
        'All (*)',
    ]

    def log_open(self, fname='', mode='w'):
        if not fname:
            fname = getSaveFileNameWithExt(self, 'Open Logfile...', self.initialdir,
                                    filter=';;'.join(self.LOG_FORMATS))
        if fname:
            self.initialdir = os.path.dirname(fname)
            self.cmd.log_open(fname, mode)

    def log_append(self):
        return self.log_open(mode='a')

    def log_resume(self):
        fname = getSaveFileNameWithExt(self, 'Open Logfile...', self.initialdir,
                                filter=';;'.join(self.LOG_FORMATS))
        if fname:
            self.initialdir = os.path.dirname(fname)
            self.cmd.resume(fname)

    def file_run(self):
        formats = [
            'All Runnable (*.pml *.py *.pym)',
            'PyMOL Command Script (*.pml)',
            'PyMOL Command Script (*.txt)',
            'Python Script (*.py *.pym)',
            'Python Script (*.txt)',
            'All Files(*)',
        ]
        fnames, selectedfilter = getOpenFileNames(
            self, 'Open file', self.initialdir, filter=';;'.join(formats))
        is_py = selectedfilter.startswith('Python')

        with PopupOnException():
            for fname in fnames:
                self.initialdir = os.path.dirname(fname)
                self.cmd.cd(self.initialdir, quiet=0)
                # detect: .py, .pym, .pyc, .pyo, .py.txt
                if is_py or re.search(r'\.py(|m|c|o|\.txt)$', fname, re.I):
                    self.cmd.run(fname)
                else:
                    self.cmd.do("@" + fname)

    def cd_dialog(self):
        dname = QFileDialog.getExistingDirectory(
            self, "Change Working Directory", self.initialdir)
        self.cmd.cd(dname or '.', quiet=0)

    def confirm_quit(self):
        QtWidgets.qApp.quit()

    def settings_edit_all_dialog(self):
        from .advanced_settings_gui import PyMOLAdvancedSettings
        if self.advanced_settings_dialog is None:
            self.advanced_settings_dialog = PyMOLAdvancedSettings(self,
                                                                  self.cmd)
        self.advanced_settings_dialog.show()

    def shortcut_menu_edit_dialog(self):
        from .shortcut_menu_gui import PyMOLShortcutMenu
        if self.shortcut_menu_filter_dialog is None:
            self.shortcut_menu_filter_dialog = PyMOLShortcutMenu(self, self.saved_shortcuts, self.cmd)
        self.shortcut_menu_filter_dialog.show()

    def show_about(self):
        msg = [
            'The PyMOL Molecular Graphics System\n',
            'Version %s' % (self.cmd.get_version()[0]),
            u'Copyright (C) Schr\xF6dinger, LLC.',
            'All rights reserved.\n',
            'License information:',
        ]

        msg.append('Open-Source Build')

        msg += [
            '',
            'For more information:',
            'https://pymol.org',
            'sales@schrodinger.com',
        ]
        QtWidgets.QMessageBox.about(self, "About PyMOL", '\n'.join(msg))

    #################
    # GUI callbacks
    #################

    def command_get(self):
        return self.lineedit.text()

    def command_set(self, v):
        return self.lineedit.setText(v)

    def command_set_cursor(self, i):
        return self.lineedit.setCursorPosition(i)

    def update_progress(self):
        progress = self.cmd.get_progress()
        if progress >= 0:
            self.progressbar.setValue(progress * 100)
            self.progressbar.show()
            self.abortbutton.show()
        else:
            self.progressbar.hide()
            self.abortbutton.hide()

    def update_feedback(self):
        self.update_progress()

        feedback = self.cmd._get_feedback()
        if feedback:
            html = colorprinting.text2html('\n'.join(feedback))
            self.browser.appendHtml(html)

            scrollbar = self.browser.verticalScrollBar()
            scrollbar.setValue(scrollbar.maximum())

        for setting in self.cmd.get_setting_updates() or ():
            if setting in self.setting_callbacks:
                current_value = self.cmd.get_setting_tuple(setting)[1][0]
                for callback in self.setting_callbacks[setting]:
                    callback(current_value)

        self.feedback_timer.start(500)

    def doPrompt(self):
        self.doTypedCommand(self.command_get())
        self.pymolwidget._pymolProcess()
        self.lineedit.clear()
        self.feedback_timer.start(0)

    ##########################
    # legacy plugin system
    ##########################

    @PopupOnException.decorator
    def initializePlugins(self):
        from pymol import plugins
        from . import mimic_tk

        self.menudict['Plugin'].clear()

        app = plugins.get_pmgapp()

        plugins.legacysupport.addPluginManagerMenuItem()

        # Redirect to Legacy submenu
        self.menudict['PluginQt'] = self.menudict['Plugin']
        self.menudict['Plugin'] = self.menudict['PluginQt'].addMenu('Legacy Plugins')
        self.menudict['Plugin'].setTearOffEnabled(True)
        self.menudict['PluginQt'].addSeparator()

        plugins.HAVE_QT = True
        plugins.initialize(app)

    def createlegacypmgapp(self):
        from . import mimic_pmg_tk as mimic
        pmgapp = mimic.PMGApp()
        pmgapp.menuBar = mimic.PmwMenuBar(self.menudict)
        return pmgapp

    def window_cmd(self, action, x, y, w, h):
        if action == 0: # hide
            self.hide()
        elif action == 1: # show
            self.show()
        elif action == 2: # position
            self.move(x, y)
        elif action == 3: # size (first two arguments)
            self.resize(x, y)
        elif action == 4: # box
            self.move(x, y)
            self.resize(w, h)
        elif action == 5: # maximize
            self.showMaximized()
        elif action == 6: # fit
            if hasattr(QtGui, 'QWindow') and self.windowHandle().visibility() in (
                    QtGui.QWindow.Maximized, QtGui.QWindow.FullScreen):
                return
            a = QtWidgets.QApplication.desktop().availableGeometry(self)
            g = self.geometry()
            f = self.frameGeometry()
            w = min(f.width(), a.width())
            h = min(f.height(), a.height())
            x = max(min(f.x(), a.right() - w), a.x())
            y = max(min(f.y(), a.bottom() - h), a.y())
            self.setGeometry(
                x - f.x() + g.x(),
                y - f.y() + g.y(),
                w - f.width() + g.width(),
                h - f.height() + g.height(),
            )
        elif action == 7: # focus
            self.setFocus(Qt.OtherFocusReason)
        elif action == 8: # defocus
            self.clearFocus()


def commandoverloaddecorator(func):
    name = func.__name__
    func.__doc__ = getattr(pymol.cmd, name).__doc__
    setattr(pymol.cmd, name, func)
    pymol.cmd.extend(func)
    return func


def SettingAction(parent, cmd, name, label='', true_value=1, false_value=0,
                  command=None):
    '''
    Menu toggle action for a PyMOL setting

    parent: parent QObject
    cmd: PyMOL instance
    name: setting name
    label: menu item text
    '''
    if not label:
        label = name

    index = cmd.setting._get_index(name)
    type_, values = cmd.get_setting_tuple(index)
    action = QtWidgets.QAction(label, parent)

    if not command:
        command = lambda: cmd.set(
            index,
            true_value if action.isChecked() else false_value,
            log=1,
            quiet=0)

    parent.setting_callbacks[index].append(
        lambda v: action.setChecked(v != false_value))

    if type_ in (
            1,  # bool
            2,  # int
            3,  # float
            5,  # color
            6,  # str
    ):
        action.setCheckable(True)
        if values[0] == true_value:
            action.setChecked(True)
    else:
        print('TODO', type_, name)

    action.triggered.connect(command)
    return action

window = None


class CommandLineEdit(QtWidgets.QLineEdit):
    '''
    Line edit widget with instant text insert on drag-enter
    '''
    _saved_pos = -1

    def dragMoveEvent(self, event):
        pass

    def dropEvent(self, event):
        if event.mimeData().hasText():
            event.acceptProposedAction()

    def dragEnterEvent(self, event):
        if not event.mimeData().hasText():
            self._saved_pos = -1
            return

        event.acceptProposedAction()

        urls = event.mimeData().urls()
        if urls and urls[0].isLocalFile():
            droppedtext = urls[0].toLocalFile()
        else:
            droppedtext = event.mimeData().text()

        pos = self.cursorPosition()
        text = self.text()
        self._saved_pos = pos
        self._saved_text = text

        self.setText(text[:pos] + droppedtext + text[pos:])
        self.setSelection(pos, len(droppedtext))

    def dragLeaveEvent(self, event):
        if self._saved_pos != -1:
            self.setText(self._saved_text)
            self.setCursorPosition(self._saved_pos)


class PyMOLApplication(QtWidgets.QApplication):
    '''
    Catch drop events on app icon
    '''
    # FileOpen event is only activated after the first
    # application state change, otherwise sys.argv would be
    # handled by Qt, we don't want that.

    def handle_file_open(self, ev):
        if ev.type() == QtCore.QEvent.ApplicationActivate:
            self.handle_file_open = self.handle_file_open_active
        return False

    def handle_file_open_active(self, ev):
        if ev.type() != QtCore.QEvent.FileOpen:
            return False

        # When double clicking a file in Finder, open it in a new instance
        if not pymol.invocation.options.reuse_helper and pymol.cmd.get_names():
            window.new_window([ev.file()])
            return True

        # pymol -I -U
        if pymol.invocation.options.auto_reinitialize:
            pymol.cmd.reinitialize()

        # PyMOL Show
        if ev.file().endswith('.psw'):
            pymol.cmd.set('presentation')
            pymol.cmd.set('internal_gui', 0)
            pymol.cmd.set('internal_feedback', 0)
            pymol.cmd.full_screen('on')

        window.load_dialog(ev.file())
        return True

    def event(self, ev):
        if self.handle_file_open(ev):
            return True
        return super(PyMOLApplication, self).event(ev)


# like pymol.internal._copy_image
def _copy_image(_self=pymol.cmd, quiet=1, dpi=-1):
    import tempfile
    fname = tempfile.mktemp('.png')

    if not _self.png(fname, prior=1, dpi=dpi):
        print("no prior image")
        return

    try:
        qim = QtGui.QImage(fname)
        QtWidgets.QApplication.clipboard().setImage(qim)
    finally:
        os.unlink(fname)

    if not quiet:
        print(" Image copied to clipboard")


def make_pymol_qicon():
    icons_dir = os.path.expandvars('$PYMOL_DATA/pymol/icons')
    return QtGui.QIcon(os.path.join(icons_dir, 'icon2.svg'))


def execapp():
    '''
    Run PyMOL as a Qt application
    '''
    global window
    global pymol

    # don't let exceptions stop PyMOL
    import traceback
    sys.excepthook = traceback.print_exception

    # use QT_OPENGL=desktop (auto-detection may fail on Windows)
    if hasattr(Qt, 'AA_UseDesktopOpenGL') and pymol.IS_WINDOWS:
        QtCore.QCoreApplication.setAttribute(Qt.AA_UseDesktopOpenGL)

    # enable 4K scaling on Windows and Linux
    if hasattr(Qt, 'AA_EnableHighDpiScaling') and not any(
            v in os.environ
            for v in ['QT_SCALE_FACTOR', 'QT_SCREEN_SCALE_FACTORS']):
        QtCore.QCoreApplication.setAttribute(Qt.AA_EnableHighDpiScaling)

    # fix Windows taskbar icon
    if pymol.IS_WINDOWS:
        import ctypes
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(
                u'com.schrodinger.pymol')

    app = PyMOLApplication(['PyMOL'])
    app.setWindowIcon(make_pymol_qicon())

    window = PyMOLQtGUI()
    window.setWindowTitle("PyMOL")

    @commandoverloaddecorator
    def viewport(w=-1, h=-1, _self=None):
        window.viewportsignal.emit(int(w), int(h))

    @commandoverloaddecorator
    def full_screen(toggle=-1, _self=None):
        from pymol import viewing as v
        toggle = v.toggle_dict[v.toggle_sc.auto_err(str(toggle), 'toggle')]
        window.toggle_fullscreen(toggle)

    import pymol.gui
    pymol.gui.createlegacypmgapp = window.createlegacypmgapp

    pymol.cmd._copy_image = _copy_image
    pymol.cmd._call_in_gui_thread = MainThreadCaller()

    # Assume GUI thread, make OpenGL context current before calling func().
    def _call_with_opengl_context_gui_thread(func):
        with window.pymolwidget:
            return func()

    # Dispatch to GUI thread and make OpenGL context current before calling func().
    pymol.cmd._call_with_opengl_context = lambda func: pymol.cmd._call_in_gui_thread(
        lambda: _call_with_opengl_context_gui_thread(func))

    window.show()
    window.raise_()

    # window size according to -W -H options
    options = pymol.invocation.options
    if options.win_xy_set:
        scale = window.pymolwidget.fb_scale
        viewport(scale * options.win_x, scale * options.win_y)

    # load plugins
    if options.plugins:
        window.initializePlugins()

    app.exec_()
