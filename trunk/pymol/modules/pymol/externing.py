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

if __name__=='pymol.externing':
    
    import os
    import pymol
    import string
    import parsing
    import threading

    from glob import glob
    from cmd import _cmd,lock,unlock,Shortcut,QuietException, \
          _feedback,fb_module,fb_mask, exp_path, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error        

    def cd(dir):
        '''
DESCRIPTION

    "cd" changes the current working directory.

USAGE

    cd <path>

SEE ALSO

    pwd, ls, system
        '''
        dir = exp_path(dir)
        os.chdir(dir)  # raises on error
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
        print os.getcwd()
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
        if pattern==None:
            pattern = "*"
        else:
            pattern = exp_path(pattern)
        if string.find("*",pattern)<0:
            lst = glob(pattern+"/*")
        else:
            lst = []
        if not len(lst):
            lst = glob(pattern)
        if len(lst):
            lst.sort()
            lst = parsing.list_to_str_list(lst)
            for a in lst:
                print a
        else:
            print " ls: Nothing found.  Is that a valid path?"
        return DEFAULT_SUCCESS

    def system(command,async=0):
        '''
DESCRIPTION

    "system" executes a command in a subshell under Unix or Windows.

USAGE

    system command 

PYMOL API

    cmd.system(string command,int sync=1)

NOTES

    async can only be specified from the Python level (not the command language)

    if async is 0 (default), then the result code from "system" is returned in r

    if async is 1, then the command is run in a separate thread whose object is
    returned

SEE ALSO

    ls, cd, pwd
        '''
        r = None
        if async:
            r = threading.Thread(target=_cmd.system,args=(str(command),1))
            r.start()
        else:
            r = _cmd.system(str(command),0)
        return r # special meaning

    def paste(): # INTERNAL
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
                if len(a):
                    new_lst.append(a)
            r = _cmd.paste(new_lst)
        if _raising(r): raise pymol.CmdException
        return r 

