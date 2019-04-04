'''
Injects PyQt replacements for the following modules into sys.modules:

tkMessageBox
tkFileDialog
tkSimpleDialog
'''

import sys
from pymol.Qt import QtWidgets


class _qtSimpleDialog:
    _NAME = 'tkSimpleDialog'
    _tkSimpleDialog = None

    @property
    def Dialog(self):
        cls = self.__class__
        if cls._tkSimpleDialog is None:
            del sys.modules[cls._NAME]
            __import__(cls._NAME)
            cls._tkSimpleDialog = sys.modules[cls._NAME]
            sys.modules[cls._NAME] = self
        return cls._tkSimpleDialog.Dialog

    def askinteger(_, title, prompt, **kw):
        r = QtWidgets.QInputDialog.getInt(None, title, prompt,
                kw.get('initialvalue', 0),
                kw.get('minvalue', -0x7fffffff),
                kw.get('maxvalue', 0x7fffffff))
        return r[0] if r[1] else None

    def askfloat(_, title, prompt, **kw):
        r = QtWidgets.QInputDialog.getDouble(None, title, prompt,
                float(kw.get('initialvalue', 0)),
                float(kw.get('minvalue', -0x7fffffff)),
                float(kw.get('maxvalue', 0x7fffffff)))
        return r[0] if r[1] else None

    def askstring(_, title, prompt, **kw):
        QLE = QtWidgets.QLineEdit
        r = QtWidgets.QInputDialog.getText(None, title, prompt,
                QLE.Password if kw.get('show') else QLE.Normal,
                kw.get('initialvalue', ''))
        return r[0] if r[1] else None


class _qtMessageBox:
    def __getattr__(self, name):
        QMB = QtWidgets.QMessageBox

        variants = {
            'askyesno': ('question', QMB.Yes, QMB.No),
            'askquestion': ('question', QMB.Yes, QMB.No),
            'askokcancel': ('question', QMB.Ok, QMB.Cancel),
            'askretrycancel': ('question', QMB.Retry, QMB.Cancel),
            'showinfo': ('information', QMB.Ok, QMB.NoButton),
            'showerror': ('critical', QMB.Ok, QMB.NoButton),
            'showwarning': ('warning', QMB.Ok, QMB.NoButton),
        }

        try:
            method, ok, others = variants[name]
        except KeyError:
            raise AttributeError(name)

        return lambda title, message, **kw: ok == getattr(QMB, method)(
                None, title, message, buttons=ok | others)


class _qtFileDialog:
    def _getfilter(self, filetypes):
        import collections
        extensions = collections.defaultdict(list)
        names = []
        for name, ext in filetypes:
            if ext.startswith('.'):
                ext = '*' + ext
            extensions[name].append(ext)
            if name not in names:
                names.append(name)
        return ';;'.join('%s (%s)' % (name, ' '.join(extensions[name]))
                for name in names)

    def askopenfilename(self, **options):
        if options.get('multiple'):
            func = QtWidgets.QFileDialog.getOpenFileNames
        else:
            func = QtWidgets.QFileDialog.getOpenFileName
        return func(None,
            options.get('title', ''),
            options.get('initialdir', ''),
            self._getfilter(options.get('filetypes', '')))[0]

    def askopenfilenames(self, **options):
        options['multiple'] = 1
        return self.askopenfilename(**options)

    def askopenfile(self, mode='r', **options):
        r = self.askopenfilename(**options)
        if options.get('multiple'):
            return [open(f, mode) for f in r]
        if not r:
            return None
        return open(r, mode)

    def askopenfiles(self, **options):
        options['multiple'] = 1
        return self.askopenfile(**options)

    def asksaveasfilename(self, **options):
        return QtWidgets.QFileDialog.getSaveFileName(None,
            options.get('title', ''),
            options.get('initialdir', ''),
            self._getfilter(options.get('filetypes', '')))[0]

    def asksaveasfile(self, mode='w', **options):
        r = self.asksaveasfilename(**options)
        if not r:
            return None
        return open(r, mode)

    def askdirectory(self, **options):
        return QtWidgets.QFileDialog.getExistingDirectory(None,
                options.get('title', ''),
                options.get('initialdir', ''))


qtMessageBox = _qtMessageBox()
qtFileDialog = _qtFileDialog()

# for all Python versions - allows plugin manager to import this with Python 3
# without importing "tkinter"
sys.modules['tkMessageBox'] = qtMessageBox
sys.modules['tkFileDialog'] = qtFileDialog

if sys.version_info[0] < 3:
    sys.modules['tkSimpleDialog'] = _qtSimpleDialog()
else:
    # injecting 'X.Y' into sys.modules without assigning the attribute
    # (import X;X.Y = ...) doesn't work. Use a meta_path solution instead.

    mapping = {
        'tkinter.messagebox': qtMessageBox,
        'tkinter.filedialog': qtFileDialog,
    }

    class MimicTkImporter:
        def load_module(self, fullname):
            # Loader method (deprecated but simple)
            m = mapping[fullname]
            if isinstance(m, str):
                import importlib
                m = importlib.import_module(m)
            sys.modules[fullname] = m
            return m

        def find_spec(self, fullname, path, target=None):
            # Python 3 Finder method (MetaPathFinder)
            if fullname in mapping:
                from importlib.machinery import ModuleSpec
                return ModuleSpec(fullname, self)
            return None

    sys.meta_path.insert(0, MimicTkImporter())
