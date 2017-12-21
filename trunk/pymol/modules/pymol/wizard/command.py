from __future__ import print_function

try:
    from itertools import izip_longest
except ImportError:
    from itertools import zip_longest as izip_longest

from pymol import cmd, stored, wizard

class Command(wizard.Wizard):
    '''
    Generic wizard for any PyMOL command.
    '''

    async_ = 1
    ignored_args = ('quiet',)

    def __init__(self, command, _self=cmd):
        import inspect

        self.input_arg = None
        self.input_value = None
        self.input_numeric = False

        wizard.Wizard.__init__(self, _self)

        self.command = command
        self.func = _self.keyword[command][0]
        spec = inspect.getargspec(self.func)
        self.args = []
        self.current = {}
        self.shortcut = {}

        self.stored_name = stored.get_unused_name('_wizard')
        setattr(stored, self.stored_name, self)
        self.varname = 'stored.' + self.stored_name

        for arg, aa in izip_longest(spec.args, _self.auto_arg):
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

        for (arg, value) in zip(reversed(spec.args), reversed(spec.defaults)):
            if not (arg in self.ignored_args or arg.startswith('_')):
                self.current[arg] = value

    def set_menu_values(self, arg, values, first=1, last=-2):
        self.menu[arg][first:last] = [
            [1, str(value), '%s.set_current("%s", %s)' % (self.varname, arg, repr(value)) ]
            for value in values
        ]

    def get_menu(self, arg):
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
