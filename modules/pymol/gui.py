#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* Copyright (c) Schrodinger, LLC.
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information.
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-*
#-*
#-*
#Z* -------------------------------------------------------------------

# abstract (external or internal) gui control interface

import pymol
cmd = __import__("sys").modules["pymol.cmd"]

def get_pmgapp():
    '''
    Returns the PMGApp instance.
    '''
    if pymol._ext_gui is None:
        pymol._ext_gui = createlegacypmgapp()
    return pymol._ext_gui

def get_qtwindow():
    '''
    Returns the PyMOLQtGUI/QMainWindow instance, or None if not available.
    '''
    try:
        from pmg_qt.pymol_qt_gui import window
        return window
    except ImportError:
        return None

def createlegacypmgapp():
    import pymol.plugins.legacysupport as m
    return m.createlegacypmgapp()

# external gui control

def ext_hide(_self=cmd):
    print('ignoring gui.ext_hide')

def ext_show(_self=cmd):
    print('ignoring gui.ext_show')

# common actions

def save_as(_self=cmd):
    qtwindow = get_qtwindow()
    if qtwindow is not None:
        qtwindow.session_save_as()

def save_image(_self=cmd):
    qtwindow = get_qtwindow()
    if qtwindow is not None:
        qtwindow.file_save_png()
