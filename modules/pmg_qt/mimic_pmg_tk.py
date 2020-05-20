'''
Mimic the pmg_tk API for plugin legacy support
'''

import sys

tkinter = None

def tkinter_init():
    global tkinter

    if hasattr(_BaseWidget_setup, '_super'):
        raise RuntimeError('tkinter init failed')

    import tkinter

    _BaseWidget_setup._super = tkinter.BaseWidget._setup
    tkinter.BaseWidget._setup = _BaseWidget_setup


# monkeypatch the "_default_root" assignment of Tkinter
def _BaseWidget_setup(self, master, cnf):
    if not master and not tkinter._default_root:
        from pymol import plugins
        tkinter._default_root = plugins.get_tk_root()
    return _BaseWidget_setup._super(self, master, cnf)


class PmwMenuBar:
    def __init__(self, menudict):
        self._menudict = menudict

    def _get_menu(self, menuName):
        try:
            return self._menudict[menuName]
        except KeyError:
            print('Error: no such menu: ' + repr(menuName))
            return None

    def addmenu(self, menuName, *args, **kw):
        self.addcascademenu('', menuName)

    def deletemenuitems(self, menuName, start, end=None):
        menu = self._get_menu(menuName)
        if menu is None:
            return

        for a in menu.actions()[start - 1:(end or start)]:
            menu.removeAction(a)

    def addmenuitem(self, menuName, itemType, statusHelp='',
                    traverseSpec=None, **kw):
        menu = self._get_menu(menuName)
        if menu is None:
            return

        if itemType == 'separator':
            menu.addSeparator()
            return

        if itemType != 'command' or 'command' not in kw:
            return

        # Wrapper for exception safety. PyMOL would crash if an exception
        # is not caught!
        def wrapper(command=kw['command']):
            try:
                command()
            except BaseException as e:
                from pymol import colorprinting
                colorprinting.print_exc([__file__])

                from pymol.Qt import QtWidgets
                QtWidgets.QMessageBox.critical(None, 'Error', str(e))

        label = kw.get('label', statusHelp)
        menu.addAction(label, wrapper)

    def addcascademenu(self, parentMenuName, menuName, statusHelp='',
                       traverseSpec=None, **kw):
        menu = self._get_menu(parentMenuName)
        if menu is None:
            return

        if menuName in self._menudict:
            raise ValueError('menu ' + repr(menuName) + ' exists')

        label = kw.get('label', statusHelp) or menuName
        menu = menu.addMenu(label)
        menu.setTearOffEnabled(True)
        self._menudict[menuName] = menu


class PMGSkin(object):
    def __init__(self, pmgapp):
        self._pmgapp = pmgapp
        self._setting = None

    @property
    def setting(self):
        if self._setting is None:
            from pmg_tk.Setting import Setting
            self._setting = Setting(self._pmgapp)

        return self._setting


class tkapp_proxy(object):
    def __init__(self, proxied, pmgapp):
        self._proxied = proxied
        self._pmgapp = pmgapp

    def __getattr__(self, name):
        return getattr(self._proxied, name)

    def call(self, tkcmd, *args):
        # suspend our own updates for commands which enter the event loop
        pause = tkcmd in ('update', 'tkwait', 'vwait')

        if pause:
            self._pmgapp._tk_update_paused += 1

        try:
            r = self._proxied.call(tkcmd, *args)
        finally:
            if pause:
                self._pmgapp._tk_update_paused -= 1

        return r


class PMGApp(object):
    def __init__(self):
        import pymol
        self._root = None
        self.pymol = pymol
        self.skin = PMGSkin(self)
        self._tk_update_paused = 0

    @property
    def root(self):
        if self._root is None:
            from pymol.Qt import QtCore

            tkinter_init()

            # create Tk instance in this thread
            self._root = tkinter.Tk()
            self._root.tk = tkapp_proxy(self._root.tk, self)
            self._root.withdraw()

            # feed Tk event loop from this thread
            timer = QtCore.QTimer()
            @timer.timeout.connect
            def _():
                if not self._tk_update_paused:
                    self._root.update()
                timer.start()
            timer.setSingleShot(True)
            timer.start(50)

            # keep reference to timer
            self._tk_update_timer = timer

            import Pmw
            Pmw.initialise(self._root)

        return self._root

    def execute(self, c):
        return eval(c) if isinstance(c, str) else c()

    def my_show(self, w, c=1):
        w.show()


# vi:expandtab:sw=4
