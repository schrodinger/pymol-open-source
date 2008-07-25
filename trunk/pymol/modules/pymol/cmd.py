
#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
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

# cmd.py 
# Python interface module for PyMol
#
# **This is the only module which should be/need be imported by 
# ** PyMol API Based Programs

# NEW CALL RETURN CONVENTIONS for _cmd.so C-layer
#

# (1) Calls into C (_cmd) should return results/status and print
#     errors and feedback (according to mask) BUT NEVER RAISE EXCEPTIONS
#     from within the C code itself.

# (2) Effective with version 0.99, standard Python return conventions
# apply, but haven't yet been fully implemented.  In summary:

#     Unless explicitly specified in the function:

#     ==> Success with no information should return None

#     ==> Failure should return a negative number as follows:
#        -1 = a general, unspecified failure
          
#     Upon an error, exceptions will be raised by the Python wrapper
#     layer if the "raise_exceptions" setting is on.

#     ==> Boolean queries should return 1 for true/yes and 0 for false/no.

#     ==> Count queries should return 0 or a positive number

# (3) If _cmd produces a specific return result, be sure to include an
#     error result as one of the possibilities outside the range of the
#     expected return value.  For example, a negative distance
#
# (4) cmd.py API wrappers can then raise exceptions and return values.
#
# NOTE: Output tweaking via the "quiet" parameter of API functions.
#
# Many PyMOL API functions have a "quiet" parameter which is used to
# customize output depending on when and where the routine was called.
#
# As defined, quiet should be 1.  Called from an external module, output
# to the console should be minimal.
#
# However, when a command is run through the parser (often manually) the
# user expects a little more feedback.  The parser will automatically
# set "quiet" to zero
#
# In rare cases, certain nonserious error or warning output should
# also be suppressed.  Set "quiet" to 2 for this behavior.


def _deferred_init_pymol_internals(_pymol):
    # set up some global session tasks

    if viewing.session_restore_views not in _pymol._session_restore_tasks:
        _pymol._session_restore_tasks.append(viewing.session_restore_views)

    if viewing.session_save_views not in _pymol._session_save_tasks:
        _pymol._session_save_tasks.append(viewing.session_save_views)

    if viewing.session_restore_scenes not in _pymol._session_restore_tasks:
        _pymol._session_restore_tasks.append(viewing.session_restore_scenes)

    if viewing.session_save_scenes not in _pymol._session_save_tasks:
        _pymol._session_save_tasks.append(viewing.session_save_scenes)

    # take care of some deferred initialization

    _pymol._view_dict_sc = Shortcut({})
    _pymol._scene_dict_sc = Shortcut({})

    # 
if __name__=='pymol.cmd':

    import traceback
    import sys

    try:
        
        import re
        import _cmd
        import string
        import thread
        import threading
        import pymol
        import os
        import parsing
        import __main__
        import time

        _pymol = pymol
        
        from shortcut import Shortcut

        from chempy import io

        #######################################################################
        # symbols for early export
        #######################################################################

        from constants import *
        from constants import _load2str

        fb_debug = sys.stderr # can redirect python debugging output elsewhere if desred...

        #--------------------------------------------------------------------
        # convenient type and result checking

        from checking import *
        from checking import _raising
        
        #-------------------------------------------------------------------
        # path expansion, including our fixes for Win32

        def _nt_expandvars(path): # allow for //share/folder$/file
            path = nt_hidden_path_re.sub(r"$$\\",path)
            return os.path.expandvars(path)
        
        if "nt" in sys.builtin_module_names:
            _expandvars = _nt_expandvars
        else:
            _expandvars = os.path.expandvars
        
        def exp_path(path):
            return _expandvars(os.path.expanduser(path))
        
        #--------------------------------------------------------------------
        # locks and threading

        reaper = None

        # the following locks are used by both C and Python to insure that no more than
        # one active thread enters PyMOL at a given time. 
        
        lock_api = pymol.lock_api
        lock_api_c = pymol.lock_api_c
        lock_api_status = pymol.lock_api_status
        lock_api_glut = pymol.lock_api_glut
        lock_api_data = pymol.lock_api_data
        lock_api_allow_flush = 1
        
        from locking import *

        #--------------------------------------------------------------------
        # status monitoring
        
        from monitoring import *
        
        #--------------------------------------------------------------------
        # Feedback

        from feedingback import *
        from feedingback import _feedback

        #--------------------------------------------------------------------
        # internal API routines

        import internal

        _adjust_coord = internal._adjust_coord
        _alt = internal._alt
        _coordset_update_spawn = internal._coordset_update_spawn
        _coordset_update_thread = internal._coordset_update_thread
        _copy_image = internal._copy_image
        _ctrl = internal._ctrl
        _do = internal._do
        _dump_floats = internal._dump_floats
        _dump_ufloats = internal._dump_ufloats
        _fake_drag = internal._fake_drag
        _get_color_sc = internal._get_color_sc
        _get_feedback = internal._get_feedback
        _interpret_color = internal._interpret_color
        _invalidate_color_sc = internal._invalidate_color_sc
        _load = internal._load
        _mpng = internal._mpng
        _object_update_spawn = internal._object_update_spawn
        _object_update_thread = internal._object_update_thread
        _png = internal._png
        _quit = internal._quit
        _ray_anti_spawn = internal._ray_anti_spawn
        _ray_hash_spawn = internal._ray_hash_spawn
        _ray_spawn = internal._ray_spawn
        _refresh = internal._refresh
        _sgi_stereo = internal._sgi_stereo
        _special = internal._special
        _validate_color_sc = internal._validate_color_sc
        _cache_get = internal._cache_get
        _cache_set = internal._cache_set
        _cache_clear = internal._cache_clear
        _cache_purge = internal._cache_purge
        _cache_mark = internal._cache_mark
        _sdof = internal._sdof
        
        # when adding, remember to also edit cmd2.py

        get_feedback = _get_feedback # legacy
        
        #######################################################################
        # now import modules which depend on the above
        #######################################################################

        import editor

        #######################################################################
        # cmd module functions...
        #######################################################################
            
        # for extending the language

        from commanding import extend, alias, dummy

        # for documentation etc

        from helping import python_help
                
        def write_html_ref(file):
            lst = globals()
            f=open(file,'w')
            head = 'H2'
            f.write("<HTML><BODY><H1>Reference</H1>")
            kees = lst.keys()
            kees.sort()
            for a in kees:
                if hasattr(lst[a],'__doc__'):
                    if (a[0:1]!='_' and
                        (a not in ['string','thread',
                                   'setup_global_locks',
                                   'real_system', 'sys','imp','glob','vl','time',
                                   'threading', 'repres','re','python_help','os',
                                   'fb_debug', 'fb_dict','ctrl','auto_arg','alt','a',
                                   'help_only', 'special','stereo_dict','toggle_dict',
                                   'palette_dict', 'types' ])):
                        doc = lst[a].__doc__
                        if is_string(doc):
                            if len(doc):
                                doc = string.strip(doc)
                                doc = string.replace(doc,"<","&lt;")
                                f.write("<HR SIZE=1><%s>"%head+a+"</%s>\n"%head)
                                f.write("<PRE>"+string.strip(doc)+"\n\n</PRE>")
            f.write("</BODY></HTML>")
            f.close()

        
        #####################################################################
        # Here is where the PyMOL Command Language and API are built.
        #####################################################################

        # first we need to import a set of symbols into the local namespace

        from api import *

        # deferred initialization

        _deferred_init_pymol_internals(pymol)
        
        # now we create the command langauge
        
        import keywords
        keyword = keywords.get_command_keywords()
        kw_list = keyword.keys()
        
        keywords.fix_list(kw_list)
        kwhash = Shortcut(kw_list)
        keywords.fix_dict(keyword)
        
        # informational or API-only functions which don't exist in the
        # PyMOL command language namespace

        help_only = keywords.get_help_only_keywords()
        help_sc = Shortcut(keyword.keys()+help_only.keys())

        # keyboard configuration
        
        import keyboard
        
        special = keyboard.get_special()

        shft_special = keyboard.get_shft_special()        
        alt_special = keyboard.get_alt_special()        
        ctrl_special = keyboard.get_ctrl_special()
        ctsh_special = keyboard.get_ctsh_special()

        ctrl = keyboard.get_ctrl()        
        alt = keyboard.get_alt()

        selection_sc = lambda sc=Shortcut,gn=get_names:sc(gn('public')+['all'])
        object_sc = lambda sc=Shortcut,gn=get_names:sc(gn('objects'))
        map_sc = lambda sc=Shortcut,gnot=get_names_of_type:sc(gnot('object:map'))
        contour_sc =  lambda sc=Shortcut,gnot=get_names_of_type:sc(
            gnot('object:mesh')+gnot('object:surface'))
        
        # Table for argument autocompletion

        import completing
        
        auto_arg = completing.get_auto_arg_list()

        color_sc = None

    except:
        print "Error: unable to initalize the pymol.cmd module"
        traceback.print_exc()
        sys.exit(0)
        
else:
    from pymol.cmd import *
    

