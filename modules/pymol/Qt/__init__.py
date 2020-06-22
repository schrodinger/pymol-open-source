"""
Wrapper for PyMOL scripts to get PySide or PyQt

Useful link for PySide/PyQt4 differences:
https://deptinfo-ensip.univ-poitiers.fr/ENS/pyside-docs/pysideapi2.html

PyQt5/PyQt4 differences:
http://pyqt.sourceforge.net/Docs/PyQt5/pyqt4_differences.html
"""

DEBUG = False

PYQT_NAME = None
QtWidgets = None

try:
    from pymol._Qt_pre import *
except ImportError:
    if DEBUG:
        print('import _Qt_pre failed')

import os

qt_api = os.environ.get('QT_API', '')

if not PYQT_NAME and qt_api in ('', 'pyqt5'):
    try:
        from PyQt5 import QtGui, QtCore, QtOpenGL, QtWidgets
        PYQT_NAME = 'PyQt5'
    except ImportError:
        if DEBUG:
            print('import PyQt5 failed')

if not PYQT_NAME and qt_api in ('', 'pyside2'):
    try:
        from PySide2 import QtGui, QtCore, QtOpenGL, QtWidgets
        PYQT_NAME = 'PySide2'
    except ImportError:
        if DEBUG:
            print('import PySide2 failed')

if not PYQT_NAME and qt_api in ('', 'pyqt4'):
    try:
        try:
            import PyQt4.sip as sip
        except ImportError:
            import sip
        sip.setapi("QString", 2)
        from PyQt4 import QtGui, QtCore, QtOpenGL
        PYQT_NAME = 'PyQt4'
    except ImportError:
        if DEBUG:
            print('import PyQt4 failed')

if not PYQT_NAME and qt_api in ('', 'pyside'):
    try:
        from PySide import QtGui, QtCore, QtOpenGL
        PYQT_NAME = 'PySide'
    except ImportError:
        if DEBUG:
            print('import PySide failed')

if not PYQT_NAME:
    raise ImportError(__name__)

# qtpy compatibility
os.environ['QT_API'] = PYQT_NAME.lower()

if QtWidgets is None:
    QtWidgets = QtGui

if hasattr(QtCore, 'QAbstractProxyModel'):
    QtCoreModels = QtCore
else:
    QtCoreModels = QtGui

if PYQT_NAME == 'PyQt4':
    QFileDialog = QtWidgets.QFileDialog
    QFileDialog.getOpenFileName = QFileDialog.getOpenFileNameAndFilter
    QFileDialog.getOpenFileNames = QFileDialog.getOpenFileNamesAndFilter
    QFileDialog.getSaveFileName = QFileDialog.getSaveFileNameAndFilter
    del QFileDialog

if PYQT_NAME[:4] == 'PyQt':
    QtCore.Signal = QtCore.pyqtSignal
    QtCore.Slot = QtCore.pyqtSlot
else:
    QtCore.pyqtSignal = QtCore.Signal
    QtCore.pyqtSlot = QtCore.Slot
    QtCore.QT_VERSION_STR = QtCore.__version__
    QtCore.QT_VERSION = (
            0x10000 * QtCore.__version_info__[0] +
            0x00100 * QtCore.__version_info__[1] +
            0x00001 * QtCore.__version_info__[2])

del qt_api
del os
