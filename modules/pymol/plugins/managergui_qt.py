from pymol.gui import get_qtwindow as getPyMOLWindow
from pymol.plugins import installation, repository

from pymol.Qt import QtGui, QtCore
from pymol.Qt import QtWidgets

from . import pref_get

Qt = QtCore.Qt

def confirm_network_access():
    '''
    Popup dialog with network access notification (only once per session)
    '''
    self = confirm_network_access
    if self.ok < 0:
        check = QtWidgets.QMessageBox.information(None, 'Info',
            'Network download has been disabled, sorry!')
        return False
    if self.ok > 0:
        return True
    if QtWidgets.QMessageBox.question(None, 'Confirm',
        'PyMOL will now download executable code from the internet!'
        ' Proceed?') == QtWidgets.QMessageBox.Yes:
        self.ok = 1
    else:
        self.ok = 0

    return self.ok

# valid values: 1=never ask, 0=ask once per session, -1=network access disabled
confirm_network_access.ok = pref_get('network_access_ok', 0)

class PluginManager(QtCore.QObject):

    def __init__(self, cmd):

        super(PluginManager, self).__init__(None)

        window = getPyMOLWindow()

        self.form = window.load_form('pluginmanager')
        self.form.b_wiki.pressed.connect(self.fetchplugin)
        self.form.e_wiki.returnPressed.connect(self.fetchplugin)
        self.form.b_local.pressed.connect(self.installplugin)
        self.form.b_startup_all.pressed.connect(self.startup_all)
        self.form.b_startup_none.pressed.connect(self.startup_none)
        self.form.e_filter.textChanged.connect(self.filter_plugins)
        self.form.c_loaded.stateChanged.connect(self.filter_plugins)
        self.form.c_startup.stateChanged.connect(self.filter_plugins)
        self.form.l_repositories.itemClicked.connect(self.repo_changed)
        self.form.bb_path_add.pressed.connect(self.add_path)
        self.form.bb_path_up.pressed.connect(self.move_path_up)
        self.form.bb_path_down.pressed.connect(self.move_path_down)
        self.form.bb_path_remove.pressed.connect(self.remove_path)
        self.form.b_add_repo.pressed.connect(self.add_repository)
        self.form.b_remove_repo.pressed.connect(self.remove_repository)
        self.reload_plugins()
        self.form.l_repo_plugins.setSelectionMode(
            QtWidgets.QAbstractItemView.ExtendedSelection)
        self.form.b_install.pressed.connect(self.install_repo_plugins)
        self.form.b_info.pressed.connect(self.info_repo_plugin)
        self.populate_repositories()
        self.populate_startup_paths()
        self.populate_preferences()
        self.form._dialog.show()

    def populate_repositories(self):
        repos = [
            'http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/',
            'https://github.com/Pymol-Scripts/Pymol-script-repo',
            # for testing
            'http://www.thomas-holder.de/projects/pymol/repository/'
        ]
        for rep in repos:
            self.form.l_repositories.addItem(rep)

    def populate_startup_paths(self):
        from . import get_startup_path
        items = get_startup_path(True)
        self.form.slb_path.addItems(items)

    def update_startup_paths(self):
        from . import set_startup_path
        items = []
        for index in range(self.form.slb_path.count()):
            items.append(self.form.slb_path.item(index).text())
        set_startup_path(items)

    def populate_preferences(self):
        from . import preferences

        self.t_preferences_keys = list(preferences)

        w = self.form.t_preferences
        w.horizontalHeader().setStretchLastSection(True)
        w.setRowCount(len(self.t_preferences_keys))

        for row, key in enumerate(self.t_preferences_keys):
            item = QtWidgets.QTableWidgetItem(key)
            item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)
            w.setItem(row, 0, item)

            item = QtWidgets.QTableWidgetItem()
            value = preferences[key]
            if isinstance(value, bool):
                item.setCheckState(Qt.Checked if value else Qt.Unchecked)
            else:
                item.setFlags(item.flags() & ~Qt.ItemIsUserCheckable)
                item.setText(str(value))
            w.setItem(row, 1, item)

        w.itemChanged.connect(self.update_preferences_item)

    def update_preferences_item(self, item):
        from . import pref_get, pref_set
        key = self.t_preferences_keys[item.row()]

        if item.flags() & Qt.ItemIsUserCheckable:
            value = item.checkState() == Qt.Checked
        else:
            value = item.data(0)

            if hasattr(value, 'toString'):  # PyQt4
                value = value.toString()

            try:
                value = type(pref_get(key, ''))(value)
            except BaseException as e:
                print(e)

        pref_set(key, value)

    def install_repo_plugins(self):
        from .installation import installPluginFromFile, get_plugdir
        items = self.form.l_repo_plugins.selectedItems()
        if len(items) == 0:
            return
        plugdir = get_plugdir(None)
        if not plugdir:
            return
        import sys, tempfile, shutil, os
        from .legacysupport import tkMessageBox
        from .repository import guess
        tmpdir = tempfile.mkdtemp()
        try:
            for item in items:
                name = item.text()
                filename = self.repo_r.copy(name, tmpdir)
                installPluginFromFile(filename, None, plugdir)
        except:
            err = str(sys.exc_info()[1])
            tkMessageBox.showinfo('Error', 'Could not install plugin ' + name + '\n\n' + err)
        finally:
            shutil.rmtree(tmpdir)
        self.reload_plugins()

    def info_repo_plugin(self):
        from . import PluginInfo
        from .installation import get_name_and_ext, extract_zipfile, zip_extensions
        from .legacysupport import tkMessageBox
        items = self.form.l_repo_plugins.selectedItems()
        if len(items) == 0:
            return
        import tempfile, shutil, os
        tmpdir = tempfile.mkdtemp()
        tmpdirs = [tmpdir]
        try:
            name = items[0].text()
            filename = self.repo_r.copy(name, tmpdir)
            name, ext = get_name_and_ext(filename)
            if ext in zip_extensions:
                tmpdir, dirnames = extract_zipfile(filename, ext)
                tmpdirs.append(tmpdir)
                name = dirnames[-1]
                filename = os.path.join(os.path.join(tmpdir, *dirnames), '__init__.py')
            info = PluginInfo(name, filename)
            self.show_plugin_info_dialog(info)
        except:
            tkMessageBox.showinfo('Error', 'Could not get plugin info')
        finally:
            for tmpdir in tmpdirs:
                shutil.rmtree(tmpdir)

    def repo_changed(self, item):
        from .repository import guess
        url = item.text()
        self.repo_r = guess(url)
        plist = self.repo_r.list()
        self.form.l_repo_plugins.clear()
        for p in plist:
            self.form.l_repo_plugins.addItem(p)

    def add_repository(self):
        repo, result = QtWidgets.QInputDialog.getText(None,
            'Add repository', 'Enter repository URL')
        if result and repo.startswith('http'):
            self.form.l_repositories.addItem(repo)

    def remove_repository(self):
        items = self.form.l_repositories.selectedItems()
        if len(items) == 0:
            return
        row = self.form.l_repositories.row(items[0])
        self.form.l_repositories.takeItem(row)
        items = self.form.l_repositories.selectedItems()
        if len(items) > 0:
            self.repo_changed(items[0])

    def show(self):
        self.form._dialog.show()
        self.form._dialog.raise_()

    def w_startup_changed(self, state):
        item = self.sender().parent()._form
        info = self.plugin_info[item]
        info.autoload = state

    def reload_plugins(self, setting=None):
        from . import plugins
        window = getPyMOLWindow()

        self.clear_plugins()

        def add_plugin_item(info):
            item = window.load_form('pluginitem', QtWidgets.QFrame())
            item._widget = item._dialog
            item._widget.setFrameStyle(QtWidgets.QFrame.Sunken)
            item._widget.setFrameShape(QtWidgets.QFrame.Panel)

            item.w_title.setText(info.name)
            item.w_version.setText(info.get_version())
            item.w_startup.setChecked(info.autoload)
            item.w_startup.stateChanged.connect(self.w_startup_changed)
            item.w_enable.pressed.connect(info.load)
            item.w_enable.pressed.connect(self.disable_load_button)
            item.w_uninstall.pressed.connect(info.uninstall)
            item.w_uninstall.pressed.connect(self.reload_plugins)
            if info.loaded:
                item.w_enable.setEnabled(False)
            if hasattr(info.module, 'settings_dialog'):
                item.w_settings.setVisible(True)
            else:
                item.w_settings.setVisible(False)
            item.w_info.pressed.connect(self.show_info)
            if info.loadtime:
                item.w_loadtime.setText("Took %.3f seconds to load" % info.loadtime)
            else:
                item.w_loadtime.setText("Not loaded")

            item._widget._form = item
            self.plugin_info[item] = info

            self.form.f_installed_layout.addWidget(item._widget)

        for info in sorted(plugins.values(), key=lambda i: i.name.lower()):
            try:
                add_plugin_item(info)
            except Exception as e:
                print(e)

        self.filter_plugins()
        self.form.f_installed_layout.addStretch()

    def clear_plugins(self):
        self.plugin_info = {}

        layout = self.form.f_installed_layout
        while layout.count():
            child = layout.takeAt(0)
            if child.widget() is not None:
                child.widget().deleteLater()

    def filter_plugins(self):
        filter_loaded = self.form.c_loaded.isChecked()
        filter_startup = self.form.c_startup.isChecked()
        filter_name = self.form.e_filter.text().lower()

        for (item, info) in self.plugin_info.items():
            if filter_name not in info.name.lower():
                vis = False
            elif filter_loaded and not info.loaded:
                vis = False
            elif filter_startup and not info.autoload:
                vis = False
            else:
                vis = True

            item._widget.setVisible(vis)

    def startup_all(self):
        for item, info in self.plugin_info.items():
            item.w_startup.setChecked(True)
            info.autoload = True

    def startup_none(self):
        for item, info in self.plugin_info.items():
            item.w_startup.setChecked(False)
            info.autoload = False

    def installplugin(self):
        from .legacysupport import installPlugin
        installPlugin(self)
        self.reload_plugins()

    def fetchplugin(self):
        if not confirm_network_access():
            return
        from .installation import installPluginFromFile
        from .repository import fetchscript
        from pymol import CmdException
        url = self.form.e_wiki.text()
        if not len(url):
            return
        import tempfile, shutil
        tmpdir = tempfile.mkdtemp()
        try:
            filename = fetchscript(url, tmpdir, False)
        except BaseException as e:
            QtWidgets.QMessageBox.critical(None, 'Error',
                'Fetching Plugin failed.\n' + str(e))
            return

        if filename:
            installPluginFromFile(filename)
        shutil.rmtree(tmpdir)
        self.reload_plugins()

    def show_info(self):
        item = self.sender().parent()._form
        info = self.plugin_info[item]
        self.show_plugin_info_dialog(info)

    def disable_load_button(self):
        item = self.sender().parent()._form
        item.w_enable.setEnabled(False)

    def show_plugin_info_dialog(self, info):
        dialog = QtWidgets.QDialog(None)
        dialog.setWindowTitle('Plugin Information')
        layout = QtWidgets.QVBoxLayout()
        table = QtWidgets.QTableWidget(0, 2)
        table.verticalHeader().hide()
        table.horizontalHeader().hide()
        table.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(table)
        dialog.setLayout(layout)

        def add_line(label, text):
            row = table.rowCount()
            table.insertRow(table.rowCount())
            table_item = QtWidgets.QTableWidgetItem(label)
            table_item.setFlags(table_item.flags() &
                 ~(Qt.ItemIsEditable))
            table.setItem(row, 0, table_item)
            table_item = QtWidgets.QTableWidgetItem(text)
            table_item.setFlags(table_item.flags() &
                 ~(Qt.ItemIsEditable))
            table.setItem(row, 1, table_item)

        add_line('Name', info.name)
        if not info.is_temporary:
            add_line('Python Module Name', info.mod_name)
            add_line('Filename', info.filename)

        metadata = info.get_metadata()
        for label, value in metadata.items():
            add_line(label, value)

        if not info.is_temporary:
            if info.loaded:
                add_line('commands', ', '.join(info.commands))
        docstring = info.get_docstring() or 'No documentation available.'

        browser = QtWidgets.QTextBrowser()
        browser.setPlainText(docstring)
        layout.addWidget(browser)

        table.resizeColumnsToContents()

        dialog.resize(600, dialog.height())
        dialog.exec_()

    def add_path(self):
        from .installation import get_default_user_plugin_path as userpath
        d = QtWidgets.QFileDialog.getExistingDirectory(None,
            'Add plugin directory', userpath())
        if len(d) > 0:
            self.form.slb_path.addItem(d)
        self.update_startup_paths()

    def move_path_up(self):
        items = self.form.slb_path.selectedItems()
        if len(items) > 0:
            row = self.form.slb_path.row(items[0])
            if row > 0:
                item = self.form.slb_path.takeItem(row)
                self.form.slb_path.insertItem(row-1, item)
        self.update_startup_paths()

    def move_path_down(self):
        items = self.form.slb_path.selectedItems()
        if len(items) > 0:
            row = self.form.slb_path.row(items[0])
            if row < self.form.slb_path.rowCount() - 1:
                item = self.form.slb_path.takeItem(row)
                self.form.slb_path.insertItem(row+1, item)
        self.update_startup_paths()

    def remove_path(self):
        items = self.form.slb_path.selectedItems()
        if len(items) > 0:
            row = self.form.slb_path.row(items[0])
            self.form.slb_path.takeItem(row)
        self.update_startup_paths()
