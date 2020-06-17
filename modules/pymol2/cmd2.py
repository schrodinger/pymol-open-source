from pymol import cmd as global_cmd
import pymol
import inspect
import itertools
import weakref


class Cmd:
    '''
    Proxy for the pymol.cmd module for multiple instances support.
    '''

    def __init__(self, _pymol, _COb):
        self._weakref = weakref.ref(self)
        self._weakrefproxy = weakref.proxy(self)

        # store parent

        self._pymol = weakref.proxy(_pymol)

        # store C object for easy access

        self._COb = _COb

        # private data

        self.color_sc = None
        self.reaper = None

        # deferred initiailization

        global_cmd._deferred_init_pymol_internals(_pymol)

        # PRIVATE FUNCTIONS (requiring '_self' as a keyword argument)

        # locking.py

        self.reaper = None

        if 1:
            # use own locks (for performance)
            self.lock_api = _pymol.lock_api
            self.lock_api_data = _pymol.lock_api_data
            self.lock_api_glut = _pymol.lock_api_glut
            self.lock_api_status = _pymol.lock_api_status
        else:
            # use global locks (for debugging)
            self.lock_api = global_cmd._pymol.lock_api
            self.lock_api_data = global_cmd._pymol.lock_api_data
            self.lock_api_glut = global_cmd._pymol.lock_api_glut
            self.lock_api_status = global_cmd._pymol.lock_api_status

        self.lock_api_allow_flush = 1
        self.lockcm = global_cmd.LockCM(self)

        # now we create the command langauge

        from pymol import keywords

        self.keyword = keywords.get_command_keywords()
        self.kw_list = list(self.keyword)

        keywords.fix_list(self.kw_list)
        self.kwhash = self.Shortcut(self.kw_list)
        keywords.fix_dict(self.keyword)

        self.help_only = keywords.get_help_only_keywords()
        self.help_sc = self.Shortcut(
            itertools.chain(self.keyword, self.help_only))


        self.selection_sc = lambda sc=self.Shortcut,gn=self.get_names:sc(gn('public')+['all'])
        self.object_sc = lambda sc=self.Shortcut,gn=self.get_names:sc(gn('objects'))
        self.map_sc = lambda sc=self.Shortcut,gnot=self.get_names_of_type:sc(gnot('object:map'))
        self.contour_sc =  lambda sc=self.Shortcut,gnot=self.get_names_of_type:sc(gnot('object:mesh')+gnot('object:surface'))
        self.group_sc = lambda sc=self.Shortcut,gnot=self.get_names_of_type:sc(gnot('object:group'))

        self.fb_action_sc = pymol.feedingback.fb_action_sc
        self.fb_module_sc = pymol.feedingback.fb_module_sc
        self.fb_mask_sc = pymol.feedingback.fb_mask_sc

        self.auto_arg = pymol.completing.get_auto_arg_list(self)
        self.color_sc = None

        # keyboard configuration

        from pymol import keyboard

        self.key_mappings = keyboard.get_default_keys(self)

# PUBLIC API METHODS which expect "self" as the first argument

    def __getattr__(self, key):
        v = getattr(global_cmd, key)
        # determine whether `v` is a callable and has a `_self` argument or
        # a keywords dictionary
        i = -1
        try:
            argspec = inspect.getfullargspec(v)
            if '_self' not in argspec.kwonlyargs:
                try:
                    i = argspec.args.index('_self')
                except ValueError:
                    if argspec.varkw is None:
                        raise TypeError
        except TypeError:
            setattr(self, key, v)
            return v
        # don't bind a circular reference into the wrapper
        # namespace, use a weak reference instead
        cmdref = self._weakref
        def wrapper(*a, **k):
            if i == -1 or len(a) <= i:
                k['_self'] = cmdref()
            return v(*a, **k)
        wrapper.__name__ = key
        setattr(self, key, wrapper)
        return wrapper
