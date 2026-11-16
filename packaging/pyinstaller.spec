# PyInstaller one-folder bundle. Build from the repository root with an
# installed (or editable) pymol in the current environment:
#
#     pyinstaller --noconfirm packaging/pymolish.spec
#
# Output: dist/pymolish/ with the pymolish executable. Data lives at
# pymolish/_internal/pymol/pymol_path so pymol.guess_pymol_path() finds it.
import glob
import os
from importlib.util import find_spec

from PyInstaller.utils.hooks import (collect_data_files, collect_submodules,
                                     copy_metadata)

# plugins under data/startup are loaded from pymol_path at runtime, not frozen
hidden = (collect_submodules('pymol', filter=lambda n: not n.startswith('pymol.plugins.startup.'))
          + collect_submodules('pmg_qt') + collect_submodules('chempy')
          + ['PySide6.QtUiTools'])

# pdbfixer and openmm (extra "fix", found through pymol/computing.py).
# openmm loads its platforms from OpenMM.libs/lib/plugins next to the package
# and reads force-field files; pdbfixer reads residue templates. Only the
# Reference and CPU platforms are shipped: the OpenCL, CUDA and HIP ones are
# LGPL and need drivers. License texts: licenses/OpenMM.txt and the metadata.
fix_datas, fix_binaries = [], []
if find_spec('pdbfixer'):
    fix_datas = (collect_data_files('openmm') + copy_metadata('pdbfixer')
                 + collect_data_files('pdbfixer', excludes=['tests/*']))
    libs = os.path.join(os.path.dirname(find_spec('openmm').origin), os.pardir,
                        'OpenMM.libs', 'lib')
    for sub in ('', 'plugins'):
        fix_binaries += [
            (f, os.path.join('OpenMM.libs', 'lib', sub))
            for f in glob.glob(os.path.join(libs, sub, '*'))
            if os.path.isfile(f) and not any(
                gpu in f.lower() for gpu in ('opencl', 'cuda', 'hip'))]

a = Analysis(
    ['launcher.py'],
    datas=[('../data', 'pymol/pymol_path/data'),
           ('../LICENSE', 'pymol/pymol_path'),
           ('../LICENSE-MIT', 'pymol/pymol_path'),
           # Qt is LGPL: a bundle that ships it must ship the text
           ('licenses', 'pymol/pymol_path/licenses')]
          + collect_data_files('pmg_qt') + collect_data_files('pymol')
          # rdkit (extra "chem", found through the import in pymol/smiles.py)
          # and the Pillow it pulls in are BSD-style: ship their license texts
          + (copy_metadata('rdkit') + copy_metadata('pillow')
             if find_spec('rdkit') else [])
          + fix_datas,
    binaries=fix_binaries,
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
exe = EXE(pyz, a.scripts, exclude_binaries=True, name='pymolish', console=True,
          # PyInstaller converts the PNG with Pillow (dev extra)
          icon='../data/pymol/icons/pymolish.png')
coll = COLLECT(exe, a.binaries, a.datas, name='pymolish')
