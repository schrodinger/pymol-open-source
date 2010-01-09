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

if __name__=='pymol.commanding':

    import thread
    import string
    import re
    import os
    import time
    import threading
    import traceback
    import parsing
    import cmd
    import pymol
    
    from cmd import _cmd, Shortcut, QuietException, \
          fb_module, fb_mask, is_list, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, is_ok, is_error, is_string

    def resume(filename, _self=cmd):
        '''
        
DESCRIPTION

    "resume" executes a log file and opens it for recording of
    additional commands.

USAGE

    resume filename

SEE ALSO

    log, log_close

    '''
        pymol=_self._pymol        
        r = DEFAULT_ERROR
        if os.path.exists(filename):
            if(re.search(r"\.py$|\.PY$|\.pym$|.PYM$",filename)):
                r = _self.do("run %s"%filename)
            else:
                r = _self.do("@%s"%filename)
        if is_ok(r):
            r = _self.do("log_open %s,a"%filename)
        if _self._raising(r,_self): raise pymol.CmdException
        return r


    class QueueFile:

        def __init__(self,queue):
            self.queue = queue

        def write(self,command):
            self.queue.put(command)

        def flush(self):
            pass

        def close(self):
            del self.queue
 
    def log_open(filename='log.pml', mode='w', _self=cmd):
        '''
DESCRIPTION

    "log_open" opens a log file for writing.

USAGE

    log_open filename

SEE ALSO

    log, log_close
    
        '''
        pymol=_self._pymol
        if not is_string(filename): # we're logging to Queue, not a file.
            pymol._log_file = QueueFile(filename)
            _self.set("logging",1,quiet=1)
        else:
            try:
                try:
                    if hasattr(pymol,"_log_file"):
                        if pymol._log_file!=None:
                            pymol._log_file.close()
                            del pymol._log_file
                except:
                    pass
                pymol._log_file = open(filename,mode)
                if _self._feedback(fb_module.cmd,fb_mask.details): # redundant
                    if mode!='a':
                        print " Cmd: logging to '%s'."%filename
                    else:
                        print " Cmd: appending to '%s'."%filename            
                if mode=='a':
                    pymol._log_file.write("\n") # always start on a new line
                if(re.search(r"\.py$|\.PY$|\.pym$|\.PYM$",filename)):
                    _self.set("logging",2,quiet=1)
                else:
                    _self.set("logging",1,quiet=1)
            except:
                print"Error: unable to open log file '%s'"%filename
                pymol._log_file = None
                _self.set("logging",0,quiet=1)
                traceback.print_exc()
                raise QuietException

    def log(text, alt_text=None, _self=cmd):
        '''
DESCRIPTION

    "log" writes a command to the log file (if one is open).

SEE ALSO

    log_open, log_close
    
        '''
        pymol=_self._pymol        
        cmd=_self
        if hasattr(pymol,"_log_file"):
            if pymol._log_file!=None:
                mode = _self.get_setting_legacy("logging")
                if mode:
                    if mode==1:
                        pymol._log_file.write(text)
                    elif mode==2:
                        if alt_text!=None:
                            pymol._log_file.write(alt_text)
                        else:
                            pymol._log_file.write("cmd.do('''%s''')\n"%string.strip(text))
                    pymol._log_file.flush()

    def log_close(_self=cmd):
        '''
DESCRIPTION

    "log_close" closes the current log file (if one is open).

USAGE

    log_close

SEE ALSO

    log, log_open
    
        '''
        pymol=_self._pymol        
        cmd=_self
        if hasattr(pymol,"_log_file"):
            if pymol._log_file!=None:
                pymol._log_file.close()
                del pymol._log_file
                _self.set("logging",0,quiet=1)
                if _self._feedback(fb_module.cmd,fb_mask.details): # redundant
                    print " Cmd: log closed."

    def cls(_self=cmd): 
        '''
DESCRIPTION

    "cls" clears the output buffer.

USAGE

    cls
    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.cls(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def splash(mode=0, _self=cmd):
        pymol=_self._pymol        
        cmd=_self
        '''
DESCRIPTION

    "splash" shows the splash screen information.

USAGE

    splash
        '''
        r = DEFAULT_ERROR
        mode = int(mode)
        if mode == 2: # query
            try:
                _self.lock(_self)
                r = _cmd.splash(_self._COb,1)
            finally:
                _self.unlock(0,_self)
        elif mode == 1: # just show PNG
            show_splash = 1
            try:
                _self.lock(_self)
                show_splash = _cmd.splash(_self._COb,1)
            finally:
                _self.unlock(0,_self)
            r = DEFAULT_SUCCESS
            if show_splash==1: # generic / open-source
                png_path = _self.exp_path("$PYMOL_PATH/data/pymol/splash.png")
            elif show_splash==2: # evaluation builds
                png_path = _self.exp_path("$PYMOL_PATH/data/pymol/epymol.png")
            else: # incentive builds
                png_path = _self.exp_path("$PYMOL_PATH/data/pymol/ipymol.png")
            if os.path.exists(png_path):
                _self.do("_ cmd.load_png('%s',0,quiet=1)"%png_path)
        else:
            if _self.get_setting_legacy("internal_feedback")>0.1:
                _self.set("text","1",quiet=1)
            print
            try:
                _self.lock(_self)
                r = _cmd.splash(_self._COb,0)
            finally:
                _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    reinit_code = {
        'everything' : 0,
        'settings' : 1,
        'store_defaults' : 2,
        'original_settings' : 3,
        'purge_defaults' : 4,
    }
    
    reinit_sc = Shortcut(reinit_code.keys())

    def reinitialize(what='everything', object='', _self=cmd):
        '''
DESCRIPTION

    "reinitialize" reinitializes the program by deleting all objects
    and restoring the default program settings.

USAGE

    reinitialize

        '''
        r = DEFAULT_ERROR
        what = reinit_code[reinit_sc.auto_err(str(what),'option')]
        try:
            _self.lock(_self)
            r = _cmd.reinitialize(_self._COb,int(what),str(object))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def sync(timeout=1.0,poll=0.05,_self=cmd):
        '''
DESCRIPTION

    "sync" is an API-only function which waits until all current
    commmands have been executed before returning.  A timeout
    can be used to insure that this command eventually returns.

PYMOL API

    cmd.sync(float timeout=1.0,float poll=0.05)

SEE ALSO

    frame
    '''
        now = time.time()
        timeout = float(timeout)
        poll = float(poll)
        # first, make sure there aren't any commands waiting...
        if _cmd.wait_queue(_self._COb): # commands waiting to be executed?
            while 1:
                if not _cmd.wait_queue(_self._COb):
                    break
                e = threading.Event() # using this for portable delay
                e.wait(poll)
                del e
                if (timeout>=0.0) and ((time.time()-now)>timeout):
                    break
        if _cmd.wait_deferred(_self._COb): 
            # deferred tasks waiting for a display event?
            if thread.get_ident() == pymol.glutThread:
                _self.refresh()
            else:
                while 1:
                    if not _cmd.wait_queue(_self._COb):
                        break
                    e = threading.Event() # using this for portable delay
                    e.wait(poll)
                    del e
                    if (timeout>=0.0) and ((time.time()-now)>timeout):
                        break
        # then make sure we can grab the API
        while 1:
            if _self.lock_attempt(_self):
                _self.unlock(_self)
                break
            e = threading.Event() # using this for portable delay
            e.wait(poll)
            del e
            if (timeout>=0.0) and ((time.time()-now)>timeout):
                break
                    
    def do(commands,log=1,echo=1,flush=0,_self=cmd):
        # WARNING: don't call this routine if you already have the API lock
        # use cmd._do instead
        '''
DESCRIPTION

    "do" makes it possible for python programs to issue simple PyMOL
    commands as if they were entered on the command line.

PYMOL API

    cmd.do( commands )

USAGE (PYTHON)

    from pymol import cmd
    cmd.do("load file.pdb")
        '''
        r = DEFAULT_ERROR
        log = int(log)
        if is_list(commands):
            cmmd_list = commands
        else:
            cmmd_list = [ commands ]
        n_cmmd = len(cmmd_list)
        if n_cmmd>1: # if processing a list of commands, defer updates
            defer = _self.get_setting_legacy("defer_updates")
            _self.set('defer_updates',1)
        for cmmd in cmmd_list:
            lst = string.split(string.replace(cmmd,chr(13),chr(10)),chr(10))
            if len(lst)<2:
                for a in lst:
                    if(len(a)):
                        try:
                            _self.lock(_self)
                            r = _cmd.do(_self._COb,a,log,echo)
                        finally:
                            _self.unlock(r,_self)
                    else:
                        r = DEFAULT_SUCCESS
            else:
                try:
                    _self.lock(_self)
                    do_flush = flush or ((thread.get_ident() == _self._pymol.glutThread)
                                         and _self.lock_api_allow_flush)
                    for a in lst:
                        if len(a):
                            r = _cmd.do(_self._COb,a,log,echo)
                            if do_flush:
                                _self.unlock(r,_self) # flushes
                                _self.lock(_self)
                        else:
                            r = DEFAULT_SUCCESS
                finally:
                    _self.unlock(r,_self)
        if n_cmmd>1:
            _self.set('defer_updates',defer)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def quit(_self=cmd):
        '''
DESCRIPTION

    "quit" terminates the program. 

USAGE

    quit

PYMOL API

    cmd.quit()
        '''
        if thread.get_ident() == pymol.glutThread:
            _self._quit(_self)
        else:
            try:
                _self.lock(_self)
                _cmd.do(_self._COb,"_ time.sleep(0.100);cmd._quit()",0,0)
                # allow time for a graceful exit from the calling thread
                try:
                    thread.exit()
                except SystemExit:
                    pass
            finally:
                _self.unlock(-1,_self=_self)
        return None

    def delete(name,_self=cmd):
        '''
DESCRIPTION

    "delete" removes an object or a selection. 

USAGE

    delete name  
    delete all   # deletes all objects

    name = name of object or selection

PYMOL API

    cmd.delete (string name = object-or-selection-name )

SEE ALSO

    remove
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.delete(_self._COb,str(name))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException      
        return r

    def extend(name,function,_self=cmd):
        
        '''
DESCRIPTION

    "extend" is an API-only function which binds a new external
    function as a command into the PyMOL scripting language.

PYMOL API

    cmd.extend(string name,function function)

PYTHON EXAMPLE

    def foo(moo=2): print moo
    cmd.extend('foo',foo)

    The following would now work within PyMOL:

    PyMOL>foo
    2
    PyMOL>foo 3
    3
    PyMOL>foo moo=5
    5
    PyMOL>foo ?
    Usage: foo [ moo ]

NOTES

    For security reasons, new PyMOL commands created using "extend" are
    not saved or restored in sessions.

SEE ALSO

    alias, api
            '''
        _self.keyword[name] = [function, 0,0,',',parsing.STRICT]
        _self.kwhash.append(name)
        _self.help_sc.append(name)

        # for aliasing compound commands to a single keyword

    def alias(name, command, _self=cmd): 
        '''
DESCRIPTION

    "alias" binds routinely-used command inputs to a new command
    keyword.

USAGE

    alias name, command

ARGUMENTS

    name = string: new keyword
    
    command = string: literal input with commands separated by semicolons.
    
EXAMPLE

    alias my_scene, hide; show ribbon, polymer; show sticks, organic; show nonbonded, solvent
    
    my_scene

NOTES

    For security reasons, aliased commands are not saved or restored
    in sessions.

SEE ALSO

    extend, api
            '''
        _self.keyword[name] = [eval("lambda :do('''%s ''')"%command.replace("'''","")), 
                               0,0,',',parsing.STRICT]
        _self.kwhash.append(name)

    def dummy(*arg):
        '''
DESCRIPTION

    This is a dummy function which returns None.
            '''
        return None
