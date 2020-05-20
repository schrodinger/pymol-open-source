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

    import os
    import pymol
    import string
    from . import parsing
    import threading
    cmd = __import__("sys").modules["pymol.cmd"]
    import traceback

    from glob import glob
    from .cmd import _cmd,lock,unlock,Shortcut,QuietException, \
          _feedback,fb_module,fb_mask, exp_path, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error

    from os import getcwd as getcwdu

    def cd(dir="~",complain=1,quiet=1):
        '''
DESCRIPTION

    "cd" changes the current working directory.

USAGE

    cd <path>

SEE ALSO

    pwd, ls, system
        '''
        dir = exp_path(dir)
        try:
            os.chdir(dir)  # raises on error
            if not quiet:
                print(" cd: now in %s" % getcwdu())
        except BaseException as e:
            if complain:
                raise pymol.CmdException(str(e))
        return DEFAULT_SUCCESS

    def pwd():
        '''
DESCRIPTION

    Print current working directory.

USAGE

    pwd

SEE ALSO

    cd, ls, system
        '''
        print(getcwdu())
        return DEFAULT_SUCCESS

    def ls(pattern=None):
        '''
DESCRIPTION

    List contents of the current working directory.

USAGE

    ls [pattern]
    dir [pattern]

EXAMPLES

    ls
    ls *.pml

SEE ALSO

    cd, pwd, system   
        '''
        if pattern is None:
            pattern = u"*"
        else:
            pattern = exp_path(pattern)
        if '*' not in pattern:
            lst = glob(os.path.join(pattern, '*'))
        else:
            lst = []
        if not len(lst):
            lst = glob(pattern)
        if len(lst):
            lst.sort()
            lst = parsing.list_to_str_list(lst)
            for a in lst:
                print(a)
        else:
            print(" ls: Nothing found.  Is that a valid path?")
        return DEFAULT_SUCCESS

    def system(command, async_=0, _self=cmd, **kwargs):
        '''
DESCRIPTION

    "system" executes a command in a subshell under Unix or Windows.

USAGE

    system command 

PYMOL API

    cmd.system(string command,int async=0)

NOTES

    async can only be specified from the Python level (not the command language)

    if async is 0 (default), then the result code from "system" is returned in r

    if async is 1, then the command is run in a separate thread whose object is
    returned

SEE ALSO

    ls, cd, pwd
        '''
        async_ = int(kwargs.pop('async', async_))

        if kwargs:
            raise pymol.CmdException('unknown argument: ' + ', '.join(kwargs))

        r = None
        if async_:
            r = threading.Thread(target=_cmd.system,args=(str(command),1))
            r.start()
        else:
            r = _cmd.system(_self._COb,str(command),0)
        return r # special meaning

    def paste(_self=cmd): # INTERNAL
        r=DEFAULT_SUCCESS
        lst = []
        if hasattr(pymol,"machine_get_clipboard"):
            lst = pymol.machine_get_clipboard()
        if len(lst):
            new_lst = []
            for a in lst:
                while len(a):
                    if ord(a[-1])>32:
                        break
                    else:
                        a=a[:-1]
                # if nothing in the queue, this special string is printed; so
                # we ignore it
                if len(a):
                    if a=="""PRIMARY selection doesn't exist or form "STRING" not defined""":
                        new_list = []
                    else:
                        new_lst.append(a)
            r = _cmd.paste(_self._COb,new_lst)
        if _raising(r,_self): raise pymol.CmdException
        return r
