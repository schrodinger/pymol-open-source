"""
Qt binding for PyMOL scripts and plugins. PySide6 only.

Plugins should import Qt from here, not from PySide6 directly.
"""

import os

from PySide6 import QtCore, QtGui, QtOpenGLWidgets, QtWidgets

PYQT_NAME = 'PySide6'
os.environ['QT_API'] = 'pyside6'  # qtpy compatibility

# ponytail: Qt5-style names kept for third-party plugins written against them
QtWidgets.QOpenGLWidget = QtOpenGLWidgets.QOpenGLWidget
QtWidgets.QActionGroup = QtGui.QActionGroup
QtWidgets.QAction = QtGui.QAction
QtWidgets.QShortcut = QtGui.QShortcut
QtCore.pyqtSignal = QtCore.Signal
QtCore.pyqtSlot = QtCore.Slot

del os
