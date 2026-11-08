# PyInstaller one-folder bundle. Build from the repository root with an
# installed (or editable) pymol in the current environment:
#
#     pyinstaller --noconfirm packaging/pymolish.spec
#
# Output: dist/pymolish/ with the pymolish executable. Data lives at
# pymolish/_internal/pymol/pymol_path so pymol.guess_pymol_path() finds it.
from PyInstaller.utils.hooks import collect_data_files, collect_submodules

# plugins under data/startup are loaded from pymol_path at runtime, not frozen
hidden = (collect_submodules('pymol', filter=lambda n: not n.startswith('pymol.plugins.startup.'))
          + collect_submodules('pmg_qt') + collect_submodules('chempy')
          + ['PySide6.QtUiTools'])

a = Analysis(
    ['launcher.py'],
    datas=[('../data', 'pymol/pymol_path/data'),
           ('../LICENSE', 'pymol/pymol_path')]
          + collect_data_files('pmg_qt') + collect_data_files('pymol'),
    hiddenimports=hidden,
    excludes=['tkinter', 'PySide6.QtWebEngineCore', 'PySide6.QtQml',
              'PySide6.QtQuick', 'PySide6.Qt3DCore', 'PySide6.QtMultimedia',
              'PySide6.QtCharts', 'PySide6.QtPdf', 'PySide6.QtDesigner'],
)
# Qt pieces PyMOL never loads; PySide6 hooks still collect them
# (DBus stays: the Linux xcb platform plugin links it)
PRUNE = ('Qt6Quick', 'Qt6Qml', 'Qt6Pdf', 'Qt6VirtualKeyboard', 'Qt6Network',
         'Qt6Designer', 'Qt/qml/', 'Qt/translations/', 'QtNetwork')
a.binaries = [b for b in a.binaries if not any(k in b[0] for k in PRUNE)]
a.datas = [d for d in a.datas if not any(k in d[0] for k in PRUNE)]

pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, exclude_binaries=True, name='pymolish', console=True)
coll = COLLECT(exe, a.binaries, a.datas, name='pymolish')
