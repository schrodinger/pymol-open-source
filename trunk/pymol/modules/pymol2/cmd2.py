from pymol import cmd as global_cmd
import pymol
import inspect

#most recently created Cmd (for now)
cmd = None

class Cmd:
    '''
    Proxy for the pymol.cmd module for multiple instances support.
    '''

    def __init__(self, _pymol, _COb):
        global cmd
        cmd = self
        # store parent
        
        self._pymol = _pymol

        # store C object for easy access
    
        self._COb = _COb
        
        # private data
    
        self.color_sc = None
        self.reaper = None
        
        # CONSTANTS (pymol/constants.py)

        self.fb_debug = global_cmd.fb_debug # this cannot be right...
        
        # deferred initiailization
        
        global_cmd._deferred_init_pymol_internals(_pymol)
        
        # PRIVATE FUNCTIONS (requiring '_self' as a keyword argument)
        
        # locking.py

        self.reaper = None

        if 1:
            # use own locks (for performance)
            self.lock_api = _pymol.lock_api
            self.lock_api_c = _pymol.lock_api_c
            self.lock_api_data = _pymol.lock_api_data            
            self.lock_api_glut = _pymol.lock_api_glut
            self.lock_api_status = _pymol.lock_api_status
        else:
            # use global locks (for debugging)
            self.lock_api = global_cmd._pymol.lock_api
            self.lock_api_c = global_cmd._pymol.lock_api_c
            self.lock_api_data = global_cmd._pymol.lock_api_data            
            self.lock_api_glut = global_cmd._pymol.lock_api_glut
            self.lock_api_status = global_cmd._pymol.lock_api_status

        self.lock_api_allow_flush = 1

        # now we create the command langauge

        from pymol import keywords

        self.keyword = keywords.get_command_keywords()
        self.kw_list = self.keyword.keys()

        keywords.fix_list(self.kw_list)
        self.kwhash = self.Shortcut(self.kw_list)
        keywords.fix_dict(self.keyword)

        self.help_only = keywords.get_help_only_keywords()
        self.help_sc = self.Shortcut(self.keyword.keys()+self.help_only.keys())

        
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
        
        self.special = keyboard.get_special(self)

        self.shft_special = keyboard.get_shft_special(self)        
        self.alt_special = keyboard.get_alt_special(self)        
        self.ctrl_special = keyboard.get_ctrl_special(self)
        self.ctsh_special = keyboard.get_ctsh_special(self)

        self.ctrl = keyboard.get_ctrl(self)        
        self.alt = keyboard.get_alt(self)
        
# PUBLIC API METHODS which expect "self" as the first argument

    def __getattr__(self, key):
        v = getattr(global_cmd, key)
        try:
            i = inspect.getargspec(v).args.index('_self')
        except:
            setattr(self, key, v)
            return v
        def wrapper(*a, **k):
            if len(a) <= i:
                k['_self'] = self
            return v(*a, **k)
        wrapper.__name__ = key
        setattr(self, key, wrapper)
        return wrapper
