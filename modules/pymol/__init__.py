'''
PyMOL Molecular Graphics System
Copyright (c) Schrodinger, Inc.

Supported ways to launch PyMOL:

  If $PYMOL_PATH is a non-default location, it must be set and exported
  before launching PyMOL.

  From a terminal:

    shell> python /path/to/pymol/__init__.py [args]

  If the 'pymol' module is in PYTHONPATH

    shell> python -m pymol [args]

  From a python main thread:

    >>> # blocks the interpreter
    >>> import pymol
    >>> pymol.launch()

  From a python main thread, spawning a new thread:

    >>> # with GUI
    >>> # THIS IS NOT SUPPORTED ON macOS
    >>> import pymol
    >>> pymol.finish_launching()

    >>> # without GUI
    >>> import pymol
    >>> pymol.finish_launching(['pymol', '-cq'])

'''

import os
import sys
import __main__

if __name__ == '__main__':

    # PyMOL launched as "python pymol/__init__.py"
    # or via execfile(".../pymol/__init__.py",...) from main
    # or as "python -m pymol.__init__"

    if 'pymol' not in sys.modules:
        # "python /abc/pymol/__init__.py" will add /abc/pymol to PYTHONPATH
        # (we don't want that), but not /abc and not the current directory (we
        # want those)

        pymol_base = os.path.dirname(os.path.realpath(__file__))
        site_packages = os.path.dirname(pymol_base)

        # remove /abc/pymol
        if pymol_base in sys.path:
            sys.path.remove(pymol_base)

        # add /abc
        if site_packages not in sys.path:
            sys.path.insert(0, site_packages)

        # add current directory
        if '' not in sys.path:
            sys.path.insert(0, '')

    # arguments default to sys.argv... but also support execfile(...)
    # from a terminal where the user could set pymol_argv
    args = getattr(__main__, "pymol_argv", None)

    # standard launch (consume main thread)
    import pymol
    sys.exit(pymol.launch(args))

IS_WINDOWS = sys.platform.startswith('win')
IS_MACOS = sys.platform.startswith('darwin')
IS_LINUX = sys.platform.startswith('linux')

import _thread as thread

import copy
import threading
import re
import time
import traceback
import math

from . import invocation
from . import colorprinting

def _init_internals(_pymol):

    # Create a temporary object "stored" in the PyMOL global namespace
    # for usage with evaluate based-commands such as alter

    _pymol.stored = Scratch_Storage()

    # Create a permanent object in the PyMOL global namespace
    # that will be picked and unpickled along with the session

    _pymol.session = Session_Storage()

    # This global will be non-None if logging is active
    # (global variable used for efficiency)

    _pymol._log_file = None

    # This global will be non-None if an external gui
    # exists. It mainly exists so that events which occur
    # in the Python thread can be handed off to the
    # external GUI thread through one or more FIFO Queues
    # (global variable used for efficiency)

    _pymol._ext_gui = None

    # lists of functions to call when saving and restoring pymol session objects

    _pymol._session_save_tasks = []
    _pymol._session_restore_tasks = []

    # cached results (as a list):
    # [ [size, (hash1, hash2, ... ), (inp1, inp2, ...), output],
    #   [size, (hash1, hash2, ... ), (inp1, inp2, ...), output],
    #   ... ]

    _pymol._cache = []

    # standard input reading thread

    _pymol._stdin_reader_thread = None

    # stored views

    _pymol._view_dict = {}
    _pymol._view_dict_sc = None

    # stored scenes

    _pymol._scene_quit_on_action = ''

    # get us a private invocation pseudo-module

    _pymol._invocation = Scratch_Storage()
    _pymol._invocation.options = copy.deepcopy(invocation.options)
    _pymol._invocation.get_user_config = invocation.get_user_config
    _pymol._invocation.parse_args = invocation.parse_args

    # these locks are to be shared by all PyMOL instances within a
    # single Python interpeter

    _pymol.lock_api = threading.RLock() # mutex for API calls from the outside
    _pymol.lock_api_status = threading.RLock() # mutex for PyMOL status info
    _pymol.lock_api_glut = threading.RLock() # mutex for GLUT avoidance
    _pymol.lock_api_data = threading.RLock() # mutex for internal data structures

def get_version_message(v=None):
    '''
    Get an informative product + version string
    '''
    if not v:
        v = _cmd.get_version()

    p = "PyMOL %s " % v[0]
    p += "Incentive Product" if invocation.options.incentive_product else \
         "Open-Source"

    if v[4]:
        p += ' (' + v[4][:10] + ')'

    if v[3]:
        p += ', ' + time.strftime('%Y-%m-%d', time.localtime(v[3]))

    return p

def guess_pymol_path():
    '''
    Guess PYMOL_PATH from typical locations and return it as string.
    '''
    init_file = os.path.abspath(__file__)

    pymol_path_candidates = [
        # $PYMOL_PATH == <site-packages>/pymol/pymol_path
        os.path.join(os.path.dirname(init_file), 'pymol_path'),

        # $PYMOL_PATH/modules/pymol/__init__.py
        re.sub(r"[\/\\]modules[\/\\]pymol[\/\\]__init__\.py[c]*$", "", init_file),

        # /usr/share/pymol
        os.path.join(sys.prefix, 'share', 'pymol'),

        # venv --system-site-packages (experimental)
        os.path.join(sys.base_prefix, 'share', 'pymol'),
    ]

    for pymol_path in pymol_path_candidates:
        if os.path.isdir(pymol_path):
            return pymol_path

    return '.'

def setup_environ():
    # guess PYMOL_PATH if unset
    if 'PYMOL_PATH' not in os.environ:
        os.environ['PYMOL_PATH'] = guess_pymol_path()

    # other PyMOL variables
    if 'PYMOL_DATA' not in os.environ:
        os.environ['PYMOL_DATA'] = os.path.join(os.environ['PYMOL_PATH'], 'data')
    if 'PYMOL_SCRIPTS' not in os.environ:
        os.environ['PYMOL_SCRIPTS'] = os.path.join(os.environ['PYMOL_PATH'], 'scripts')
    os.environ['TUT'] = os.path.join(os.environ['PYMOL_DATA'], 'tut')

    # set Tcl/Tk environment if we ship it in ext/lib
    pymol_path = os.environ['PYMOL_PATH']
    for varname, dirname in [
            ('TCL_LIBRARY', 'tcl8.5'),
            ('TK_LIBRARY', 'tk8.5')]:
        dirname = os.path.join(pymol_path, "ext", "lib", dirname)
        if os.path.isdir(dirname):
            os.environ[varname] = dirname

def exec_str(self, string):
    '''
    Execute string in "self" namespace (used from C)
    '''
    try:
        exec(string, self.__dict__, self.__dict__)
    except Exception:
        traceback.print_exc()
    return None

def exec_deferred(self):
    '''
    Execute the stuff from invocations.options.deferred
    '''
    try:
        from socket import error as socket_error
    except ImportError:
        socket_error = None
        print('import socket failed')

    cmd = self.cmd
    _pymol = cmd._pymol

    # read from stdin (-p)
    if self.invocation.options.read_stdin and not _pymol._stdin_reader_thread:
        try:
            t = _pymol._stdin_reader_thread = \
                    threading.Thread(target=cmd._parser.stdin_reader)
            t.setDaemon(1)
            t.start()
        except:
            traceback.print_exc()

    # do the deferred stuff
    try:
        if cmd.ready():
            cmd.config_mouse(quiet=1)
            for a in self.invocation.options.deferred:
                if a[0:4] == "_do_":
                    cmd.do(a[4:])
                else:
                    cmd.load(a, quiet=0)
    except CmdException as e:
        colorprinting.error(str(e))
        colorprinting.error(
            " Error: Argument processing aborted due to exception (above).")
    except socket_error:
        # this (should) only happen if we're opening a PWG file on startup
        # and the port is busy.  For now, simply bail...
        cmd.wizard("message",["Socket.error: ","",
                              "\\999Assigned socket in use.","",
                              "\\779Is PyMOL already launched?","",
                              "\\966Shutting down..."])
        cmd.refresh()
        cmd.do("time.sleep(2);cmd.quit()")

def adapt_to_hardware(self):
    '''
    optimize for (or workaround) specific hardware
    '''
    cmd = self.cmd

    vendor, renderer, version = cmd.get_renderer()

    # Quadro cards don't support GL_BACK in stereo contexts
    if vendor.startswith('NVIDIA'):
        if 'Quadro' in renderer:
            if invocation.options.show_splash:
                print(" Adapting to Quadro hardware.")
            cmd.set('stereo_double_pump_mono', 1)

    elif vendor.startswith('Mesa'):
        if renderer[0:18]=='Mesa GLX Indirect':
            pass

    elif vendor.startswith('ATI'):
        if renderer[0:17] == 'FireGL2 / FireGL3':  # obsolete ?
            if invocation.options.show_splash:
                print(" Adapting to FireGL hardware.")
            cmd.set('line_width', 2, quiet=1)

        if IS_WINDOWS:
            if sys.getwindowsversion()[0] > 5:
                # prevent color corruption by calling glFlush etc.
                cmd.set('ati_bugs', 1)

        if 'Radeon HD' in renderer:
            if invocation.options.show_splash:
                print(" Adjusting settings to improve performance for ATI cards.")

            if cmd.get_setting_int("use_shaders")==0:
                # limit frame rate to 30 fps to avoid ATI "jello"
                # where screen updates fall way behind the user.
                cmd.set("max_ups", 30)

    elif vendor.startswith('Microsoft'):
        if renderer[0:17] == 'GDI Generic':
            cmd.set('light_count', 1)
            cmd.set('spec_direct', 0.7)

    elif vendor.startswith("Intel"):
        if "Express" in renderer:
            if invocation.options.show_splash:
                print(" Disabling shaders for Intel Express graphics")
            cmd.set("use_shaders", 0)

    elif (' R300 ' in vendor # V: X.Org R300 Project, R: Gallium 0.4 on ATI RV370
            ):
        if invocation.options.show_splash:
            print(" Detected blacklisted graphics driver.  Disabling shaders.")
        cmd.set("use_shaders", 0)

    # find out how many processors we have, and adjust hash
    # table size to reflect available RAM

    try:
        import multiprocessing
        ncpu = multiprocessing.cpu_count()
        if ncpu > 1:
             cmd.set("max_threads", ncpu)
             if invocation.options.show_splash:
                  print(" Detected %d CPU cores."%ncpu, end=' ')
                  print(" Enabled multithreaded rendering.")
    except:
        pass

    # store our adapted state as default
    cmd.reinitialize("store")

def launch_gui(self):
    '''
    Launch if requested:
    - external GUI
    '''
    pymol_path = os.getenv('PYMOL_PATH', '')

    try:
        poll = IS_MACOS

        if self.invocation.options.external_gui == 3:
            if 'DISPLAY' not in os.environ:
                os.environ['DISPLAY'] = ':0.0'

        if self.invocation.options.external_gui in (1, 3):
            __import__(self.invocation.options.gui)
            sys.modules[self.invocation.options.gui].__init__(self, poll,
                    skin = self.invocation.options.skin)

            # import plugin system
            import pymol.plugins

    except:
        traceback.print_exc()

def prime_pymol():
    '''
    Set the current thread as the glutThread
    '''
    global glutThread

    if not glutThread:
        glutThread = thread.get_ident()


def _launch_no_gui():
    import pymol2

    p = pymol2.SingletonPyMOL()
    p.start()

    # TODO sufficient?
    while (p.idle() or p.getRedisplay() or
            invocation.options.keep_thread_alive or
            cmd.get_modal_draw() or
            cmd.get_setting_int('keep_alive') or
            cmd._pymol._stdin_reader_thread is not None):
        p.draw()

    # TODO needed?
    cmd.sync()
    p.stop()


def launch(args=None, block_input_hook=0):
    '''
    Run PyMOL with args

    Only returns if we are running pretend GLUT.
    '''
    if args is None:
        args = sys.argv
    invocation.parse_args(args)

    if invocation.options.gui == 'pmg_qt':
        if invocation.options.no_gui:
            return _launch_no_gui()
        elif invocation.options.testing:
            return pymol._cmd.test2()

        try:
            from pmg_qt import pymol_qt_gui
            return pymol_qt_gui.execapp()
        except ImportError as ex:
            print(f'Qt not available ({ex}), using GLUT/Tk interface')
            invocation.options.gui = 'pmg_tk'

    prime_pymol()
    _cmd.runpymol(None, block_input_hook)

def finish_launching(args=None):
    '''
    Start the PyMOL process in a thread

    THIS IS NOT SUPPORTED ON macOS
    '''
    global glutThreadObject

    if cmd._COb is not None:
        return

    import pymol

    # legacy
    if args is None:
        args = getattr(pymol, 'pymol_argv', None)
    if args is None:
        args = getattr(__main__, 'pymol_argv', sys.argv)

    if True:
        # run PyMOL in thread
        invocation.options.keep_thread_alive = 1
        cmd.reaper = threading.current_thread()
        glutThreadObject = threading.Thread(target=launch,
                args=(list(args), 1))
        glutThreadObject.start()

    e = threading.Event()

    # wait for the C library to initialize
    while cmd._COb is None:
        e.wait(0.01)

    # make sure symmetry module has time to start...
    while not hasattr(pymol, 'xray'):
        e.wait(0.01)

class CmdException(Exception):
    '''
    Exception type for PyMOL commands
    '''
    label = "Error"
    def __init__(self, message='', label=None):
        self.message = message
        if message:
            self.args = (message,)
        if label:
            self.label = label
    def __str__(self):
        return " %s: %s" % (self.label, self.message)

class IncentiveOnlyException(CmdException):
    '''
    Exception type for features that are not available in Open-Source PyMOL
    '''
    label = "Incentive-Only-Error"
    def __init__(self, message=''):
        if not message:
            try:
                funcname = sys._getframe(1).f_code.co_name
                message = '"%s" is not available in Open-Source PyMOL' % (funcname,)
            except:
                message = 'Not available in Open-Source PyMOL'
        message += '\n\n' \
                    '    Please visit http://pymol.org if you are interested in the\n' \
                    '    full featured "Incentive PyMOL" version.\n'
        super(IncentiveOnlyException, self).__init__(message)

class Scratch_Storage:
    '''
    Generic namespace
    '''
    def __reduce__(self):
        # for loading Python 3 (new-style class) pickle with Python 2
        return (self.__class__, (), self.__dict__)

    def get_unused_name(self, prefix='tmp'):
        '''
        Get an unused name from this namespace
        '''
        i = 1
        while True:
            name = prefix + str(i)
            if not hasattr(self, name):
                setattr(self, name, None)
                return name
            i += 1

class Session_Storage:
    '''
    Generic namespace
    '''
    def __reduce__(self):
        # for loading Python 3 (new-style class) pickle with Python 2
        return (self.__class__, (), self.__dict__)


def _colortype(cmd):
    # backwards compatible color index type for iterate, which used
    # to expose colors as RGB tuples
    get_color_tuple = cmd.get_color_tuple
    class Color(int):
        def __getitem__(self, i):
            return get_color_tuple(self)[i]
        def __len__(self):
            return 3

    return Color


######### VARIABLES ############################

glutThread = None

######### ENVIRONMENT ##########################

setup_environ()

# initialize instance-specific module/object internals
_init_internals(sys.modules[__name__])

# get X-window support (machine_get_clipboard)
if 'DISPLAY' in os.environ:
    from .xwin import *

########## C MODULE ############################

import pymol._cmd
_cmd = sys.modules['pymol._cmd']

get_capabilities = _cmd.get_capabilities

from . import cmd

cmd._COb = None

try:
    import epymol
except ImportError:
    pass

########## WORKAROUND TO PREVENT "import cmd" ##############################
# Previous versions of PyMOL did relative imports and thus allowd
# "import cmd" in pymol scripts to import the pymol.cmd module. To be more
# strict and for compatibility with python3 we use absolute imports now,
# which unfortunately will import an unrelated "cmd" module from the default
# python library, and even worse will corrupt the pymol namespace with it.
# The following causes an import error for "import cmd":

class _NoCmdFinder:
    def find_spec(self, fullname, path=None, target=None):
        if path is None and fullname == 'cmd':
            msg = 'use "from pymol import cmd" instead of "import cmd"'
            print('Warning: {}'.format(msg))
        return None
    find_module = find_spec

sys.meta_path.insert(0, _NoCmdFinder())

########## LEGACY PRINT STATEMENT FOR PYMOL COMMAND LINE ###################

if True:
    def _print_statement(*args, **_):
        '''Legacy Python-2-like print statement for the PyMOL command line'''
        kw = {}
        if args and args[0].startswith('>>'):
            kw['file'] = eval(args[0][2:])
            args = args[1:]
        if args and not args[-1]:
            kw['end'] = ' '
            args = args[:-1]
        args = [eval(a) for a in args]
        print(*args, **kw)

    cmd.extend('print', _print_statement)
