'''
The "PMGApp" object handed to plugin __init_plugin__(app) functions when the
Qt GUI is running. app.menuBar maps the old Pmw menu calls onto Qt menus.
'''


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


class PMGApp:
    def __init__(self):
        import pymol
        self.pymol = pymol
        self.menuBar = None

    @property
    def root(self):
        raise RuntimeError(
            'The Tk GUI was removed in pymolish-studio 4.3. '
            'Build plugin windows with pymol.Qt widgets instead.')

    def execute(self, c):
        return eval(c) if isinstance(c, str) else c()

    def my_show(self, w, c=1):
        w.show()
