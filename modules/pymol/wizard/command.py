from itertools import zip_longest as izip_longest

from pymol import cmd, stored, wizard

class Command(wizard.Wizard):
    '''
    Generic wizard for any PyMOL command.
    '''

    async_ = 1
    ignored_args = ('quiet',)

    def __getstate__(self):
        d = super(Command, self).__getstate__()
        d['shortcut'] = {}
        return d

    @property
    def func(self):
        return self.cmd.keyword[self.command][0]

    def __init__(self, command=None, _self=cmd):

        self.input_arg = None
        self.input_value = None
        self.input_numeric = False

        wizard.Wizard.__init__(self, _self)

        self.stored_name = stored.get_unused_name('_wizard')
        setattr(stored, self.stored_name, self)
        self.varname = 'stored.' + self.stored_name

        self.set_command(command, False)

    def set_command(self, command, refresh=True):
        self.command = command

        if command is None:
            return

        import inspect

        sig = inspect.signature(self.func, follow_wrapped=False)
        self.args = []
        self.current = {}
        self.shortcut = {}

        for param, aa in izip_longest(sig.parameters.values(), self.cmd.auto_arg):
            if param is None:
                break

            if param.kind != inspect.Parameter.POSITIONAL_OR_KEYWORD:
                break

            arg = param.name
            if arg.startswith('_') or arg in self.ignored_args:
                continue
            self.args.append(arg)
            self.menu[arg] = [
                [2, arg.title(), ''],
                [0, '', ''],
                [1, 'Enter value...', '%s.set_input_arg("%s")' % (self.varname, arg)]
            ]
            try:
                self.shortcut[arg] = aa[command][0]
            except:
                pass

        for param in sig.parameters.values():
            if param.default == inspect.Parameter.empty:
                continue

            arg = param.name
            if not (arg in self.ignored_args or arg.startswith('_')):
                self.current[arg] = param.default

        if refresh:
            self.cmd.refresh_wizard()

    def set_menu_values(self, arg, values, first=1, last=-2):
        self.menu[arg][first:last] = [
            [1, str(value), '%s.set_current("%s", %s)' % (self.varname, arg, repr(value)) ]
            for value in values
        ]

    def get_menu(self, arg):
        if arg == '_commands':
            from pymol import keywords
            from pymol.helping import python_help
            return [
                [1, str(command), self.varname + f'.set_command({command!r})' ]
                for (command, value) in keywords.get_command_keywords().items()
                if not command.startswith('_') and value[0] != python_help
            ]

        if arg not in self.menu:
            return None
        try:
            sc = self.shortcut[arg]()
            self.set_menu_values(arg, sc.keywords[:20])
        except KeyError:
            pass
        return self.menu[arg]

    def set_current(self, k, v):
        self.current[k] = v
        self.cmd.refresh_wizard()

    def run(self):
        if self.async_:
            self.cmd.async_(self.func, **self.current)
        else:
            self.func(**self.current)

    def cleanup(self):
        delattr(stored, self.stored_name)

    def get_panel(self):
        if self.command is None:
            return [
                [3, 'Select a command...', '_commands'],
                [2, 'Done', 'cmd.set_wizard()'],
            ]

        if self.input_arg:
            return [
                [1, 'Input', ''],
                [2, 'Apply', self.varname + '.apply_input()'],
                [2, 'Cancel ', self.varname + '.set_input_arg()'],
            ]
        panel = [
            [1, self.command.title() + ' Wizard', ''],
        ] + [
            [3, '%s: %s' % (arg, self.current.get(arg, '')), arg]
            for arg in self.args
        ] + [
            [2, 'Run', self.varname + '.run()'],
            [2, 'Done', 'cmd.set_wizard()'],
        ]
        return panel

    def get_prompt(self):
        if not self.input_arg:
            return None
        return [r'Please enter value for \999' + self.input_arg +
                r'\---: \990' + self.input_value]

    def get_event_mask(self):
        if self.input_arg:
            return self.event_mask_key
        return 0

    def set_input_arg(self, arg=None, numeric=False):
        self.input_arg = arg
        if arg:
            v = self.current.get(arg, '')
            self.input_value = str(v) if isinstance(v, (str, int, float)) else ''
            self.input_numeric = numeric
        self.cmd.refresh_wizard()

    def apply_input(self):
        self.current[self.input_arg] = self.input_value
        self.set_input_arg()

    def do_key(self, k, x, y, m):
        if k in (8, 127): # delete
            self.input_value = self.input_value[:-1]
        elif k in (10, 13):
            self.apply_input()
        elif k == 27: # escape
            self.set_input_arg()
        elif 31 < k < 127 and not self.input_numeric or k == 46 or 47 < k < 58:
            self.input_value = self.input_value + chr(k)
        else:
            print(" Warning: invalid key")
        self.cmd.refresh_wizard()
        return 1
