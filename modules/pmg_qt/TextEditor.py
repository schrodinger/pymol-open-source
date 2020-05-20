'''
Simple Text Editor
'''

import os
import sys
import importlib

try:
    from pymol.Qt import QtCore, QtWidgets, QtGui
    from pymol.Qt.utils import connectFontContextMenu, getMonospaceFont
except ImportError:
    from PyQt5 import QtCore, QtWidgets, QtGui
    connectFontContextMenu = lambda *a: None
    getMonospaceFont = lambda: QtGui.QFont()


class TextEditor(QtWidgets.QMainWindow):

    def sizeHint(self):
        return QtCore.QSize(600, 400)

    def _write(self, handle):
        content = self._get()
        handle.write(content)
        self._savedcontent = content

    def _open(self, filename):
        self.filename = filename or ''
        if filename and os.path.exists(filename):
            with open(filename, 'r') as handle:
                content = handle.read()
        else:
            content = ''
        self._set(content)

        if filename.endswith('.py'):
            self.setSyntax('python')
        elif filename.endswith('.pml') or filename.endswith('pymolrc'):
            self.setSyntax('pml')
        else:
            self.setSyntax()

        if filename:
            self.setWindowTitle('%s (%s)' % (os.path.basename(filename),
                os.path.dirname(filename)))

    def setSyntax(self, filetype='plain'):
        self.syntaxactions[filetype].setChecked(True)

        if self.highlight is not None:
            self.highlight.setDocument(None)

        try:
            m = importlib.import_module(__name__.rsplit('.', 1)[0]
                    + '.syntax.' + filetype)
        except ImportError as e:
            return

        self.highlight = m.Highlighter(self.text.document())

    def _get(self):
        return self.text.toPlainText()

    def _set(self, content):
        self.text.setPlainText(content)
        self._savedcontent = self._get()

    def doSaveAs(self, *args):
        fname, selectedfilter = QtWidgets.QFileDialog.getSaveFileName(
            None, 'Save As...', os.path.dirname(self.filename))
        if not fname:
            return

        with open(fname, 'w') as handle:
            self._write(handle)
            self.filename = handle.name

    def doSave(self, *args):
        if not self.filename:
            return self.doSaveAs()
        with open(self.filename, 'w') as handle:
            self._write(handle)

    def doOpen(self, *args):
        if not self.check_ask_save():
            return

        fnames = QtWidgets.QFileDialog.getOpenFileNames(None, 'Open file')[0]
        if fnames:
            self._open(fnames[0])

    def check_ask_save(self):
        QMessageBox = QtWidgets.QMessageBox
        if self._get() != self._savedcontent:
            ok = QMessageBox.question(None, "Save?", "Save changes?",
                                      QMessageBox.Yes | QMessageBox.No |
                                      QMessageBox.Cancel, QMessageBox.Yes)
            if ok == QMessageBox.Yes:
                self.doSave()
            elif ok == QMessageBox.Cancel:
                return False
        return True

    def closeEvent(self, event):
        if not self.check_ask_save():
            event.ignore()
            return
        self.close()

    def __init__(self, parent=None, filename='', title='Text Editor'):
        super(TextEditor, self).__init__()

        self.highlight = None

        self.root = self
        self.root.setWindowTitle(title)

        menubar = self.root.menuBar()
        filemenu = menubar.addMenu("File")
        filemenu.addAction("Open", self.doOpen, QtGui.QKeySequence("Ctrl+O"))
        filemenu.addAction("Save", self.doSave, QtGui.QKeySequence("Ctrl+S"))
        filemenu.addAction("Save as ...", self.doSaveAs,
                           QtGui.QKeySequence("Ctrl+Shift+S"))

        syntaxmenu = menubar.addMenu("Syntax")
        syntaxgroup = QtWidgets.QActionGroup(self)
        self.syntaxactions = {}

        for label in ['Python', 'PML', 'Plain Text']:
            key = label.split()[0].lower()
            action = syntaxmenu.addAction(label,
                    lambda t=key: self.setSyntax(t))
            syntaxgroup.addAction(action)
            action.setCheckable(True)
            self.syntaxactions[key] = action

        self.text = QtWidgets.QPlainTextEdit()
        self.text.setFont(getMonospaceFont())
        self.root.setCentralWidget(self.text)

        connectFontContextMenu(self.text)

        self._open(filename)

        self.root.show()
        self.root.raise_()


def edit_pymolrc(app=None):
    try:
        import pymol
        pymolrc_list = pymol.invocation.options.pymolrc
    except ImportError:
        pymolrc_list = []

    if len(pymolrc_list or ()) < 2:
        _edit_pymolrc(app, pymolrc_list)
        return

    s, ok = QtWidgets.QInputDialog.getItem(
        None, 'Select pymolrc file', 'Active pymolrc files:', pymolrc_list)

    if ok:
        _edit_pymolrc(app, [s])


def _edit_pymolrc(app, _list=()):
    try:
        pymolrc = _list[0]
    except (TypeError, IndexError):
        if sys.platform.startswith('win'):
            pymolrc = os.path.expandvars(r'$HOMEDRIVE$HOMEPATH\pymolrc.pml')
        else:
            pymolrc = os.path.expandvars(r'$HOME/.pymolrc')

        pymolrc, ok = QtWidgets.QInputDialog.getText(
            None, 'Create new pymolrc?', 'Filename of new pymolrc',
            QtWidgets.QLineEdit.Normal, pymolrc)

        if not ok:
            return

    if pymolrc:
        TextEditor(None, pymolrc)


if __name__ == '__main__':
    try:
        filename = sys.argv[1]
    except:
        filename = ''
    app = QtWidgets.QApplication(['Test'])
    edit_pymolrc()
    app.exec_()
