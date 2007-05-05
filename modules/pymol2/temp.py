#        self.color = lambda s, *a, **k: ( setitem(k,'_self',s),  apply(global_cmd.color, a, k))[1]
    def old_color(self, color, selection="(all)", quiet=1, flags=0):
        from pymol import cmd
        from pymol.cmd import _cmd,lock,unlock,Shortcut,QuietException,_raising, \
          _feedback,fb_module,fb_mask, \
          repres,repres_sc, is_string, is_list, is_ok, is_error, \
          toggle_dict,toggle_sc,stereo_dict,stereo_sc, \
          palette_dict ,palette_sc, window_dict, window_sc, \
          safe_list_eval, lock_without_glut, DEFAULT_ERROR, DEFAULT_SUCCESS

        '''
DESCRIPTION

    "color" changes the color of objects or atoms.

USAGE

    color color [, selection]

ARGUMENTS

    color = a color or ramp name

    selection = a selection-expression or name-pattern
    
PYMOL API

    cmd.color(string color, string selection, int quiet)

EXAMPLE 

    color cyan
    color yellow, chain A
        '''
        # preprocess selection
        selection = selector.process(selection)
        color = cmd._interpret_color(str(color))
        #
        r = DEFAULT_ERROR
        try:
            cmd.lock()
            r = _cmd.color(self._c_self,str(color),str(selection),int(flags),int(quiet))
        finally:
            cmd.unlock(r)
        if _raising(r): raise QuietException
        return 
