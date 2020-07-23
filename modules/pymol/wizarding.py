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

if True:

    import pymol
    import sys
    cmd = __import__("sys").modules["pymol.cmd"]
    from .cmd import _cmd,lock,unlock,Shortcut,QuietException,_raising, \
          _feedback,fb_module,fb_mask, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error

    import pickle as cPickle
    import traceback

    class WizardError(Exception):
        pass

    def _wizard(name,arg,kwd,replace,_self=cmd):
        r = DEFAULT_ERROR
        from . import wizard
        try:
            full_name = 'pymol.wizard.'+name
            __import__(full_name)
        except ImportError:
            print("Error: Sorry, couldn't import the '"+name+"' wizard.")
        else:
            mod_obj = sys.modules[full_name]
            if mod_obj:
                oname = name.capitalize()
                r = DEFAULT_SUCCESS
                if hasattr(mod_obj,oname):
                    kwd['_self']=_self
                    try:
                        wiz = getattr(mod_obj,oname)(*arg, **kwd)
                    except TypeError as e:
                        # e.g. missing argument
                        raise pymol.CmdException(str(e))
                    except WizardError as e:
                        from pymol.wizard.message import Message
                        wiz = Message("Error: %s" % str(e), _self=_self)
                    if wiz:
                        _self.set_wizard(wiz,replace)
                        _self.do("_ refresh_wizard")
                else:
                    print("Error: Sorry, couldn't find the '"+oname+"' class.")
            else:
                print("Error: Sorry, couldn't import the '"+name+"' wizard.")
        return r

    def wizard(name=None, *arg, _self=cmd, **kwd):
        '''
DESCRIPTION

    "wizard" launches on of the built-in wizards.  There are special
    Python scripts which work with PyMOL in order to obtain direct user
    interaction and easily peform complicated tasks.

USAGE

    wizard name

PYMOL API

    cmd.wizard(string name)

EXAMPLE

    wizard distance  # launches the distance measurement wizard
    '''
        r = DEFAULT_ERROR
        if name is None:
            _self.set_wizard()
            r = DEFAULT_SUCCESS
        else:
            name = str(name)
            if name.lower() == 'distance': # legacy compatibility
                name = 'measurement'
            r = _wizard(name,arg,kwd,0,_self=_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def replace_wizard(name=None, *arg, _self=cmd, **kwd):
        '''
DESCRIPTION

    "replace_wizard" is an unsupported internal command.
    
    '''
        r = DEFAULT_ERROR
        if name is None:
            _self.set_wizard()
            r = DEFAULT_SUCCESS
        else:
            r = _wizard(name,arg,kwd,1,_self=_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def set_wizard(wizard=None,replace=0,_self=cmd): # INTERNAL
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.set_wizard(_self._COb,wizard,replace)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def set_wizard_stack(stack=[],_self=cmd): # INTERNAL
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.set_wizard_stack(_self._COb,stack)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def refresh_wizard(_self=cmd): # INTERNAL
        '''
DESCRIPTION

    "refresh_wizard" is in unsupported internal command.
    
    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.refresh_wizard(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def dirty_wizard(_self=cmd): # INTERNAL
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.dirty_wizard(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def get_wizard(_self=cmd): # INTERNAL
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_wizard(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def get_wizard_stack(_self=cmd): # INTERNAL
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_wizard_stack(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def session_save_wizard(session,_self=cmd):
        # double-pickle so that session file is class-independent
        stack = _self.get_wizard_stack()
        session['wizard']=cPickle.dumps(stack,1)
        return 1

    def session_restore_wizard(session,_self=cmd):
        if session is not None:
            version = session.get('version', 0)
            if 'wizard' in session:
                from chempy.io import pkl
                try:
                    wizards = pkl.fromString(session['wizard'])
                    for wiz in wizards:
                        wiz.cmd = _self
                        wiz.migrate_session(version)
                    _self.set_wizard_stack(wizards)
                except Exception as e:
                    print(e)
                    print("Session-Warning: unable to restore wizard.")
        return 1
