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

from pymol.shortcut import Shortcut
from pymol.parser import __file__ as _parser_filename

if True:

    import sys
    if True:
        import _thread as thread
        import urllib.request as urllib2
        from io import FileIO as file

    import inspect
    import shlex
    import builtins
    from enum import Enum
    if sys.version_info >= (3, 11):
        from enum import StrEnum
    from functools import wraps
    from pathlib import Path
    from textwrap import dedent
    from typing import get_args, Union, Any, get_origin, get_type_hints
    from types import UnionType

    import re
    import os
    import time
    import threading
    import traceback
    from . import colorprinting
    from . import parsing
    cmd = sys.modules["pymol.cmd"]
    import pymol

    from .cmd import _cmd, QuietException, \
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

    re_fetch = re.compile(r'(\bfetch\s+\w[^;\r\n\'"]+)')

    class LogFile(file):
        def _append_async0(self, m):
            s = m.group()
            if 'async' in s:
                return s
            return s.rstrip(',') + ', async=0'
        def write(self, s):
            s = re_fetch.sub(self._append_async0, s)
            file.write(self, s.encode())

    def log_open(filename='log.pml', mode='w', _self=cmd):
        '''
DESCRIPTION

    "log_open" opens a log file for writing.

USAGE

    log_open [ filename [, mode ]]

ARGUMENTS

    filename = str: file to write to (.pml or .py) {default: log.pml}

    mode = w/a: "w" to open an empty log file, "a" to append {default: w}

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
                        if pymol._log_file is not None:
                            pymol._log_file.close()
                            del pymol._log_file
                except:
                    pass
                pymol._log_file = LogFile(filename,mode)
                if _self._feedback(fb_module.cmd,fb_mask.details): # redundant
                    if mode!='a':
                        print(" Cmd: logging to '%s'."%filename)
                    else:
                        print(" Cmd: appending to '%s'."%filename)
                if mode=='a':
                    pymol._log_file.write("\n") # always start on a new line
                if(re.search(r"\.py$|\.PY$|\.pym$|\.PYM$",filename)):
                    _self.set("logging",2,quiet=1)
                else:
                    _self.set("logging",1,quiet=1)
            except:
                print("Error: unable to open log file '%s'"%filename)
                pymol._log_file = None
                _self.set("logging",0,quiet=1)
                traceback.print_exc()
                raise QuietException

    def log(text, alt_text=None, _self=cmd):
        '''
DESCRIPTION

    "log" writes a command to the log file (if one is open).

    `text` and/or `alt_text` must include the terminating line feed.

ARGUMENTS

    text = str: PyMOL command (optional if alt_text is given)
    alt_text = str: Python expression (optional)

SEE ALSO

    log_open, log_close
    
        '''
        # See also equivalent C impelemtation: PLog

        log_file = getattr(_self._pymol, "_log_file", None)

        if log_file is None:
            return

        mode = _self.get_setting_int("logging")

        if mode == 1:
            # .pml
            if not text and alt_text:
                text = '/' + alt_text
        elif mode == 2:
            # .py
            if alt_text:
                text = alt_text
            elif text.startswith('/'):
                text = text[1:]
            else:
                text = "cmd.do(" + repr(text.strip()) + ")\n"
        else:
            return

        if text:
            log_file.write(text)
            log_file.flush()

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
            if pymol._log_file is not None:
                pymol._log_file.close()
                del pymol._log_file
                _self.set("logging",0,quiet=1)
                if _self._feedback(fb_module.cmd,fb_mask.details): # redundant
                    print(" Cmd: log closed.")

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

    def _load_splash_image(filename, url, _self=cmd):
        import tempfile
        import struct

        tmp_filename = ""
        contents = None

        if url:
            try:
                handle = urllib2.urlopen(url)
                contents = handle.read()
                handle.close()

                # png magic number
                if contents[:4] != b'\x89\x50\x4e\x47':
                    raise IOError

                shape = struct.unpack('>II', contents[16:24])

                tmp_filename = tempfile.mktemp('.png')
                with open(tmp_filename, 'wb') as handle:
                    handle.write(contents)

                filename = tmp_filename
            except IOError:
                pass

        if os.path.exists(filename) and not _self.get_names():
            # hide text splash
            print()

            # fit window to image
            try:
                if not contents:
                    contents = open(filename, 'rb').read(24)
                shape = struct.unpack('>II', contents[16:24])
                scale = _self.get_setting_int('display_scale_factor')
                _self.viewport(shape[0] * scale, shape[1] * scale)
            except Exception as e:
                print(e)

            # load image
            _self.load_png(filename, 0, quiet=1)

        if tmp_filename:
            os.unlink(tmp_filename)

    def splash(mode=0, _self=cmd):
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
            png_url = ""
            if show_splash==1: # generic / open-source
                png_path = _self.exp_path("$PYMOL_DATA/pymol/splash.png")
            elif show_splash==2: # evaluation builds
                png_path = _self.exp_path("$PYMOL_DATA/pymol/epymol.png")
            elif show_splash==3: # edu builds
                png_path = _self.exp_path("$PYMOL_DATA/pymol/splash_edu.png")
                png_url = "http://pymol.org/splash/splash_edu_2.png"
            else: # incentive builds
                png_path = _self.exp_path("$PYMOL_DATA/pymol/ipymol.png")

            t = threading.Thread(target=_load_splash_image, args=(png_path, png_url, _self))
            t.setDaemon(1)
            t.start()
        else:
            if _self.get_setting_int("internal_feedback") > 0:
                _self.set("text","1",quiet=1)
            print()
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
        for t in async_threads:
            t.join()

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
                    colorprinting.warning('cmd.sync() timed out (wait_queue)')
                    break
        if _cmd.wait_deferred(_self._COb):
            # deferred tasks waiting for a display event?
            if _self.is_gui_thread():
                _self.refresh()
            else:
                while 1:
                    if not _cmd.wait_queue(_self._COb):
                        break
                    e = threading.Event() # using this for portable delay
                    e.wait(poll)
                    del e
                    if (timeout>=0.0) and ((time.time()-now)>timeout):
                        colorprinting.warning('cmd.sync() timed out (wait_deferred)')
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
                colorprinting.warning('cmd.sync() timed out (lock_attempt)')
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
        log = int(log)
        cmmd_list = commands if is_list(commands) else [commands]
        cmmd_list = [a for cmmd in cmmd_list for a in cmmd.splitlines() if a]
        n_cmmd = len(cmmd_list)
        if n_cmmd>1: # if processing a list of commands, defer updates
            defer = _self.get_setting_int("defer_updates")
            _self.set('defer_updates',1)
        if flush or (_self.is_gui_thread() and _self.lock_api_allow_flush):
            for a in cmmd_list:
                with _self.lockcm:
                    _cmd.do(_self._COb,a,log,echo)
        else:
            with _self.lockcm:
                for a in cmmd_list:
                    _cmd.do(_self._COb,a,log,echo)
        if n_cmmd>1:
            _self.set('defer_updates',defer)

    def quit(code=0, _self=cmd):
        '''
DESCRIPTION

    "quit" terminates the program. 

USAGE

    quit [code]

ARGUMENTS

    code = int: exit the application with status "code" {default: 0}

PYMOL API

    cmd.quit(int code)
        '''
        code = int(code)
        if _self.is_gui_thread():
            _self._quit(code, _self)
        else:
            try:
                _self.lock(_self)
                _cmd.do(_self._COb,"_ time.sleep(0.100);cmd._quit(%d)" % (code),0,0)
                # allow time for a graceful exit from the calling thread
                try:
                    thread.exit()
                except SystemExit:
                    pass
            finally:
                _self.unlock(-1,_self=_self)
        return None

    def delete(name, *, _self=cmd):
        '''
DESCRIPTION

    "delete" removes objects and named selections

USAGE

    delete name

ARGUMENTS

    name = name(s) of object(s) or selection(s), supports wildcards (*)

EXAMPLES

    delete measure*     # delete all objects which names start with "measure"
    delete all          # delete all objects and selections

PYMOL API

    cmd.delete (string name = object-or-selection-name )

SEE ALSO

    remove
    delete_states
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.delete(_self._COb,str(name))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def delete_states(name: str, states: str, *, _self=cmd):
        '''
DESCRIPTION

USAGE

    delete_states name, states

ARGUMENTS

    name = name(s) of object(s), supports wildcards (*)

    states = string: space separated list of state numbers or ranges

EXAMPLES

    delete_state 1nmr, 1-5     # delete states 1 to 5 from 1nmr
    delete_state *, 1-3 10-40  # deletes states 1 to 3 and 10 to 40 from all applicable objects

PYMOL API

    cmd.delete_states(string name = object-name, string states = states string)

SEE ALSO

    delete
        '''

        def get_state_list(states_str):
            output=[]
            input_str = re.sub(r"\s"," ", states_str)
            input_str = input_str.replace("-"," -")
            input_str = input_str.replace("- ","-")
            input_str = input_str.strip().split()

            last = -1
            for x in input_str:
                if not x.isdigit():
                    if x.startswith("-"):
                        direction = 1
                        current = last
                        last = int(x[1:]) - 1
                        if last < current:
                            direction = -1
                        while current != last:
                            current += direction
                            output.append(str(current))
                else:
                    val = int(x) - 1
                    output.append(str(val))
                    last = val
            return output

        with _self.lockcm:
            output = get_state_list(states)
            states_list = sorted(set(map(int, output)))
            return _cmd.delete_states(_self._COb, name, states_list)
    
        
    class ArgumentParsingError(ValueError):
        "Error on argument parsing."

        def __init__(self, arg_name, message):
            message = dedent(message).strip()
            if arg_name:
                s = f"Failed at parsing '{arg_name}'. {message}"
            else:
                s = message
            super().__init__(s)
    
    
    def _into_types(var, type, value):
        
        # Untyped string
        if type == Any:
            return value
                
        # Boolean flags
        elif type is bool:
            if isinstance(value, bool):
                return value
            trues = ["yes", "1", "true", "on", "y"]
            falses = ["no", "0", "false", "off", "n"]
            if value.lower() in trues:
                return True
            elif value.lower() in falses:
                return False
            else:
                raise ArgumentParsingError(
                    var,
                    f"Can't parse {value!r} as bool."
                    f" Supported true values are {', '.join(trues)}."
                    f" Supported false values are {', '.join(falses)}."
                )
        
        # Types from typing module
        elif origin := get_origin(type):

            if origin in {Union, UnionType}:
                funcs = get_args(type)
                for func in funcs:
                    try:
                        return _into_types(None, func, value)
                    except:
                        continue
                raise ArgumentParsingError(
                    var,
                    f"Can't parse {value!r} into {type}."
                    f" The parser tried each union type and none was suitable."
                )
            
            elif issubclass(origin, tuple):
                funcs = get_args(type)
                if funcs:
                    values = shlex.split(value)
                    if len(funcs) > 0 and len(funcs) != len(values):
                        raise ArgumentParsingError(
                            var,
                            f"Can't parse {value!r} into {type}."
                            f" The number of tuple arguments are incorrect."
                        )
                    try:
                        return tuple(_into_types(None, f, v) for f, v in zip(funcs, values))
                    except:
                        raise ArgumentParsingError(
                            var,
                            f"Can't parse {value!r} into {type}."
                            f" One or more tuple values are of incorrect types."
                        )
                else:
                    return tuple(shlex.split(value))

            elif issubclass(origin, list):
                funcs = get_args(type)
                if len(funcs) == 1:
                    func = funcs[0]
                    return [_into_types(None, func, a) for a in shlex.split(value)]
                return shlex.split(value)
        
        elif sys.version_info >= (3, 11) and issubclass(type, StrEnum):
            try:
                return type(value)
            except:
                names = [e.value for e in list(type)]
                raise ArgumentParsingError(
                    var,
                    f"Invalid value for {type.__name__}."
                    f" Accepted values are {', '.join(names)}."
                )

        # Specific types must go before other generic types
        #   isinstance(type, builtins.type) comes after
        elif issubclass(type, Enum):
            value = type.__members__.get(value)
            if value is None:
                raise ArgumentParsingError(
                    var,
                    f"Invalid value for {type.__name__}."
                    f" Accepted values are {', '.join(type.__members__)}."
                )
            return value
        
        # Generic types must accept str as single argument to __init__(s)
        elif isinstance(type, builtins.type):
            try:
                return type(value)
            except Exception as exc:
                raise ArgumentParsingError(
                    var,
                    f"Invalid value {value!r} for custom type {type.__name__}."
                    f" The type must accept str as the solo argument to __init__(s)."
                ) from exc


    def new_command(name, function=None, _self=cmd):

        if function is None:
            name, function = name.__name__, name

        # docstring text, if present, should be dedented
        if function.__doc__ is not None:
            function.__doc__ = dedent(function.__doc__).strip()

        # Resolve strings into real class objects (PEP 563).
        try:
            resolved_hints = get_type_hints(
                function,
                globalns=sys.modules[function.__module__].__dict__
            )
        except Exception:
            resolved_hints = function.__annotations__

        # Analysing arguments
        sign = inspect.signature(function)
        
        # Inner function that will be callable every time the command is executed
        @wraps(function)
        def inner(*args, **kwargs):
            caller = sys._getframe(1).f_code.co_filename
            # It was called from command line or pml script, so parse arguments
            if caller == _parser_filename:
                # special _self argument
                kwargs.pop("_self", None)
                new_kwargs = {}
                for var, param in sign.parameters.items():
                    if var in kwargs:
                        value = kwargs[var]
                        # special 'quiet' argument
                        if var == 'quiet' and isinstance(value, int):
                            new_kwargs[var] = bool(value)
                        else:
                            actual_type = resolved_hints.get(var, param.annotation)
                            new_kwargs[var] = _into_types(var, actual_type, value)
                    else:
                        if param.default is sign.empty:
                            raise ArgumentParsingError(f"Unknow argument '{var}'.")
                defaults = {
                    k: v.default for k, v in sign.parameters.items()
                    if v.default is not sign.empty
                }
                final_kwargs = {
                    **defaults,
                    **new_kwargs
                }
                return function(**final_kwargs)

            # It was called from Python, so pass the arguments as is
            else:
                return function(*args, **kwargs)
        
        _self.keyword[name] = [inner, 0, 0, ',', parsing.STRICT]
        _self.kwhash.append(name)
        _self.help_sc.append(name)
        
        # Accessor to the original function so bypass the stack extraction.
        # The purpose is optimization (loops, for instance).
        inner.func = inner.__wrapped__ 
        return inner


    def extend(name, function=None, _self=cmd):

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
        if function is None:
            name, function = name.__name__, name
        _self.keyword[name] = [function, 0,0,',',parsing.STRICT]
        _self.kwhash.append(name)
        _self.help_sc.append(name)
        return function

        # for aliasing compound commands to a single keyword

    def extendaa(*arg, _self=cmd):
        '''
DESCRIPTION

    API-only function to decorate a function as a PyMOL command with
    argument auto-completion.

EXAMPLE

    @cmd.extendaa(cmd.auto_arg[0]['zoom'])
    def zoom_organic(selection='*'):
        cmd.zoom('organic & (%s)' % selection)
        '''
        auto_arg = _self.auto_arg
        def wrapper(func):
            name = func.__name__
            _self.extend(name, func)
            for (i, aa) in enumerate(arg):
                if i == len(auto_arg):
                    auto_arg.append({})
                if aa is not None:
                    auto_arg[i][name] = aa
            return func
        return wrapper

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

    cmd.extend, api
            '''
        _self.keyword[name] = [eval("lambda :do('''%s ''')"%command.replace("'''","")),
                               0,0,',',parsing.STRICT]
        _self.kwhash.append(name)

    async_threads = []

    def async_(func, *args, _self=cmd, **kwargs):
        '''
DESCRIPTION

    Run function threaded and show "please wait..." message.
        '''
        from .wizard.message import Message

        wiz = Message(['please wait ...'], dismiss=0, _self=_self)

        try:
            _self.set_wizard(wiz)
        except:
            wiz = None

        if isinstance(func, str):
            func = _self.keyword[func][0]

        def wrapper():
            async_threads.append(t)
            try:
                func(*args, **kwargs)
            except (pymol.CmdException, cmd.QuietException) as e:
                if e.args:
                    print(e)
            finally:
                if wiz is not None:
                    try:
                        _self.set_wizard_stack([w
                            for w in _self.get_wizard_stack() if w != wiz])
                    except:
                        _self.do('_ wizard')
                    else:
                        _self.refresh_wizard()

                async_threads.remove(t)

        t = threading.Thread(target=wrapper)
        t.setDaemon(1)
        t.start()
