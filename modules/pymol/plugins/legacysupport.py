'''
PyMOL Plugins Engine, Legacy Support

Provides the "PMGApp"-like object passed to a plugin's __init_plugin__(app)
and Qt-backed stand-ins for the tkMessageBox/tkFileDialog helpers that older
plugins used. Tk itself is not available.

(c) 2011-2012 Thomas Holder, PyMOL OS Fellow
License: BSD-2-Clause
'''

import os
import pymol
from . import startup

__all__ = [
    'startup',
    'get_pmgapp',
]


def get_pmgapp():
    '''
    Returns the PMGApp-like instance (Qt menu adapter, or a headless no-op).
    '''
    import pymol.gui
    return pymol.gui.get_pmgapp()


def installPlugin(parent=None):
    '''
    Open a file dialog to install a plugin. To install a directory, point to
    its __init__.py file.
    '''
    from pymol.Qt import QtWidgets
    from .installation import zip_extensions, installPluginFromFile

    filters = 'Plugins (*.py %s);;All Files (*)' % ' '.join(
        '*.' + ext for ext in zip_extensions)
    ofile = QtWidgets.QFileDialog.getOpenFileName(
        None, 'Install Plugin', os.getcwd(), filters)[0]
    if ofile:
        installPluginFromFile(ofile)


plugin_manager_panel = None


def addPluginManagerMenuItem():

    def plugin_manager():
        global plugin_manager_panel
        if not plugin_manager_panel:
            from .managergui_qt import PluginManager
            plugin_manager_panel = PluginManager(None)
        plugin_manager_panel.show()

    pymol.plugins.addmenuitem('Plugin Manager', plugin_manager)
    pymol.plugins.addmenuitem('-', None)


def createlegacypmgapp():
    '''
    Headless PMGApp stand-in: menu calls are no-ops, execute runs inline.
    '''
    app = pymol.Scratch_Storage()
    app.root = None
    app.menuBar = pymol.Scratch_Storage()
    app.menuBar.addmenuitem = \
    app.menuBar.deletemenuitems = \
    app.menuBar.addcascademenu = lambda *x, **y: None
    app.execute = lambda c: eval(c) if isinstance(c, str) else c()
    return app


class _QtDialogs:
    '''
    Lazy proxy to pmg_qt.mimic_tk.qtMessageBox / qtFileDialog. Without a
    running QApplication (headless) message boxes print and answer yes,
    file dialogs return an empty string.
    '''

    def __init__(self, attr, headless_result):
        self._attr = attr
        self._headless_result = headless_result

    def __getattr__(self, name):
        from pymol.Qt import QtWidgets
        if QtWidgets.QApplication.instance() is None:
            # ponytail: headless auto-confirms; no interactive console prompt
            def headless(title='', message='', **kw):
                if message:
                    print(' %s: %s' % (title, message))
                return self._headless_result
            return headless
        from pmg_qt import mimic_tk
        return getattr(getattr(mimic_tk, self._attr), name)


tkMessageBox = _QtDialogs('qtMessageBox', True)
tkFileDialog = _QtDialogs('qtFileDialog', '')
