# How do I launch PyMOL?

# THE SUPPORTED WAY:

# "python pymol/__init__.py" in an environment in which $PYMOL_PATH
# points to the main PyMOL directory and $PYTHONPATH includes
# $PYMOL_PATH/modules or where the contents of $PYMOL_PATH/modules
# have been installed in a standard location such as
# /usr/lib/python2.1/site-packages

# THE UNSUPPORTED/EXPERIMENTAL WAY:

# Method 2: "import pymol" from within a Python program in an
# environment where $PYMOL_PATH points to the main PyMOL directory
# and $PYTHONPATH includes $PYMOL_PATH/modules or where the contents of
# $PYMOL_PATH/modules have been installed in a standard location such
# as /usr/lib/python2.1/site-packages
#
# NOTE: with method 2, you should call pymol.finish_launching()
# before using any PYMOL API functions

# NOTE: with both methods, you should be able to get away with not
# specifying PYMOL_PATH if there is a subdirectory pymol_path located
# in the "pymol" modules directory which points to the main
# pymol directory

# NOTE: If you attempt to "import pymol" in the interactive Python
# prompt, PYTHON may crash with a "GC object already tracked" error
# message.  This is not a bug in PyMOL, as it can be produced with a
# trivial tcl/tk program.  (Tentatively WORKED AROUND in 0.99beta18!)

import copy
import __main__
if __name__!='__main__':
    import invocation
    
# Global variable "__main__.pymol_launch" tracks how we're launching PyMOL:
#
# 0: old way, now obsolete (e.g. "python launch_pymol.py")
# 1: new way, consume main thread: (e.g. from "python pymol/__init__.py")
# 2: new way, spawn own thread: (e.g."import pymol; pymol.finish_launching()")
# 3: dry run -- just get PyMOL environment information
# 4: monolithic (embedded) PyMOL.  Prime, but don't start.
# 5: Python embedded launch from within the PyMOL API

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
    # The entry 'None' represents the PyMOL C-API function call

    _pymol._session_save_tasks = [ None ]
    _pymol._session_restore_tasks = [ None ]

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

    _pymol._scene_dict = {}
    _pymol._scene_dict_sc = None
    _pymol._scene_order = []
    _pymol._scene_counter = 1
    _pymol._scene_quit_on_action = ''

    # get us a private invocation pseudo-module

    _pymol._invocation = Scratch_Storage()
    _pymol._invocation.options = copy.deepcopy(invocation.options)
    _pymol._invocation.get_user_config = invocation.get_user_config
    _pymol._invocation.parse_args = invocation.parse_args

    # these locks are to be shared by all PyMOL instances within a
    # single Python interpeter
        
    _pymol.lock_api = threading.RLock() # mutex for API calls from the outside
    _pymol.lock_api_c = threading.RLock() # mutex for C management of python threads
    _pymol.lock_api_status = threading.RLock() # mutex for PyMOL status info
    _pymol.lock_api_glut = threading.RLock() # mutex for GLUT avoidance
    _pymol.lock_api_data = threading.RLock() # mutex for internal data structures
    
if hasattr(__main__,'pymol_launch'):
    pymol_launch = __main__.pymol_launch
else:
    pymol_launch = 2 # default is standard threaded launch (unless overridden)

if pymol_launch != 3: # if this isn't a dry run

    import thread 
    import threading 
    import os
    import sys
    import re
    import string 
    import time
    import traceback
    import math

    # try to set PYMOL_PATH if unset...

    if __name__!='__main__':
        if not os.environ.has_key("PYMOL_PATH"):
            try:
                pymol_file = __file__
                # first, see if we've got "site-packages/pymol/pymol_path"
                # which would the case from a DISTUTILS install
                if (pymol_file[0:1] not in [ '\\', '/' ]) and pymol_file[1:2]!=':': 
                    pymol_file = os.getcwd()+"/"+pymol_file # make path absolute

                pymol_path = re.sub(r"[\/\\][^\/\\]*$","/pymol_path",pymol_file)

                if os.path.isdir(pymol_path):
                    os.environ['PYMOL_PATH'] = pymol_path
                # that didn't work, so check ther reverse situation "/modules/pymol/__init__.py"
                # which would right for an RPM install or simply import pymol with PYTHONPATH set
                else:
                    pymol_path = re.sub(r"[\/\\]modules[\/\\]pymol[\/\\]__init__\.py[c]*$","",pymol_file)
                    if os.path.isdir(pymol_path):
                        os.environ['PYMOL_PATH'] = pymol_path
            except NameError:
                pass

    # now start the launch process...

    if __name__=='__main__':

        # PyMOL launched as "python pymol/__init__.py"
        # or via execfile(".../pymol/__init__.py",...) from main

        if not hasattr(__main__,"pymol_argv"):
            __main__.pymol_argv = sys.argv

        pymol_launch = -1 # non-threaded launch import flag
        # NOTE: overrides current value (if any)
        
        try_again = 0

        try:
            import pymol
        except ImportError:
            try_again = 1

        if try_again: # insert directory with this file into search path
            try:
                sys.exc_clear()
                sys.path.insert(0,os.path.dirname(os.path.dirname(__file__)))
            except:
                pass
            import pymol
        
        pymol_launch = 1 # consume main thread for use by PyMOL
    
    elif pymol_launch == -1: # just passing through

        if hasattr(__main__,"pymol_argv"):
            pymol_argv = __main__.pymol_argv
        else:
            pymol_argv = [ "pymol", "-q" ]

    elif pymol_launch == 2: # standard threaded launch

        if hasattr(__main__,"pymol_argv"):
            pymol_argv = __main__.pymol_argv
        else:
            # suppresses startup messages with "import pymol"      
            pymol_argv = [ "pymol", "-q" ] 

    elif pymol_launch==4:

        if hasattr(__main__,"pymol_argv"):
            pymol_argv = __main__.pymol_argv
        else:
            pymol_argv = [ "pymol"]

    # PyMOL __init__.py

    if (__name__=='pymol') and not globals().has_key('_once'):

        # don't ever redefine these symbols...
        
        _once = None

        # Python exception type for PyMOL commands
        
        class CmdException:
            def __init__(self,args=None):
                self.args = args
        
        class Scratch_Storage:
            pass

        class Session_Storage:
            pass

        # initialize instance-specific module/object internals

        _init_internals(sys.modules['pymol'])
        
        # special handling for win32

        if sys.platform=='win32':
            # include modules directory (if it isn't already and it exists)
            loc2 = os.environ['PYMOL_PATH']+'/modules'
            if os.path.exists(loc2):
                if loc2 not in sys.path:
                    sys.path.insert(0,loc2)
            # include installed numpy
            loc1 = os.environ['PYMOL_PATH']+'/modules/numeric'
            if os.path.exists(loc1):
                sys.path.insert(0,loc1)

        if '' not in sys.path: # make sure cwd is in path like normal Python
            sys.path.insert(0,'') 

        sys.setcheckinterval(1) # maximize responsiveness

        # auto-detect bundled FREEMOL (if present)
        
        if not os.environ.has_key("FREEMOL"):
            test_path = os.path.join(os.environ['PYMOL_PATH'],"freemol")
            if os.path.isdir(test_path):
                os.environ['FREEMOL'] = test_path
                
        # include FREEMOL's libpy in sys.path (if present)
        
        if os.environ.has_key("FREEMOL"):
            freemol_libpy = os.path.join(os.environ['FREEMOL'],"libpy")
            if os.path.isdir(freemol_libpy):
                if freemol_libpy not in sys.path:
                    sys.path.append(freemol_libpy)
                
        def exec_str(self,st):
            try:
                exec st in self.__dict__, self.__dict__
            except StandardError:
                traceback.print_exc()
            return None

        def exec_deferred(self):
            import socket
            cmd=self.cmd
            if self.invocation.options.read_stdin:
                try:
                    _stdin_reader_thread = threading.Thread(target=cmd._parser.stdin_reader)
                    _stdin_reader_thread.setDaemon(1)
                    _stdin_reader_thread.start()
                except:
                    traceback.print_exc()
            try:
                socket_error = 0
                if cmd.ready():
                    cmd.config_mouse(quiet=1)
                    for a in self.invocation.options.deferred:
                        if a[0:4]=="_do_":
                            cmd.do(a[4:])
                        elif re.search(r"pymol\.py$",a):
                            pass
                        elif re.search(r"\.py$|\.pym|\.pyc$",a,re.I):
                            cmd.do("_ run %s" % a)
                        elif cmd.file_ext_re.search(a):
                            cmd.load(a,quiet=0)
                        elif re.search(r"\.pml$",a,re.I):
                            cmd.do("_ @%s" % a)
                        else:
                            cmd.load(a,quiet=0)
            except CmdException:
                traceback.print_exc()
                print "Error: Argument processing aborted due to exception (above)."
            except socket.error:
                socket_error = 1

            if socket_error:
                # this (should) only happen if we're opening a PWG file on startup
                # and the port is busy.  For now, simply bail...
                cmd.wizard("message",["Socket.error: ","",
                                      "\\999Assigned socket in use.","",
                                      "\\779Is PyMOL already launched?","",
                                      "\\966Shutting down..."])
                cmd.refresh()
                cmd.do("time.sleep(2);cmd.quit()")

        def adapt_to_hardware(self):
            cmd=self.cmd

            # optimize for VISTA
            if sys.platform[0:3]=='win':
                if sys.getwindowsversion()[0]>5:
                    # improves performance 
                    cmd.set('texture_fonts') 
                    
            # optimize for (or workaround) specific hardware
            (vendor,renderer,version) = cmd.get_renderer()

            # Quadro cards don't support GL_BACK in stereo contexts
            if vendor[0:6]=='NVIDIA':
                if 'Quadro' in renderer:
                    if invocation.options.show_splash:
                        print " Adapting to Quadro hardware."
                    cmd.set('stereo_double_pump_mono',1)                    

            elif vendor[0:4]=='Mesa':
                if renderer[0:18]=='Mesa GLX Indirect':
                    pass

            elif vendor[0:9]=='Parallels':
                if renderer[0:8]=='Parallel':
                    pass
                    # this was critical for older Parallels
                    # but actually slows down current versions
                    # cmd.set('texture_fonts',1) 

            elif vendor[0:3]=='ATI':
                if renderer[0:17]=='FireGL2 / FireGL3':  # obsolete ?
                    if invocation.options.show_splash:
                        print " Adapting to FireGL hardware."
                    cmd.set('line_width','2',quiet=1)            

                if sys.platform[0:3]=='win':
                    if sys.getwindowsversion()[0]>5:
                        # prevent color corruption by calling glFlush etc.
                        cmd.set('ati_bugs',1) 
                        
                if 'Radeon HD' in renderer:
                    #print " Note: Radeon HD cards tend not to run PyMOL well."
                    #print " Use nVidia or Intel instead, if OpenGL glitches occur."
                    print " Adjusting settings to improve performance for ATI cards."

                    # use display lists to minimize use of OpenGL
                    # immediate mode rendering (unreasonably slow on
                    # Radeon HD cards!)
                    cmd.set("use_display_lists") 

                    # disable line smooothing to prevent various
                    # bizarre screen-update and drawing artifacts
                    cmd.unset("line_smooth")

                    # limit frame rate to 30 fps to avoid ATI "jello"
                    # where screen updates fall way behind the user.
                    cmd.set("max_ups",30) 

            elif vendor[0:9]=='Microsoft':
                if renderer[0:17]=='GDI Generic':
                    cmd.set('light_count',1)
                    cmd.set('spec_direct',0.7)

            # find out how many processors we have, and adjust hash
            # table size to reflect available RAM

            try:
                ncpu = 1
                if sys.platform=='darwin':
                    if os.path.exists("/usr/sbin/sysctl"):
                        f=os.popen("/usr/sbin/sysctl hw.ncpu hw.physmem")
                        l=f.readlines()
                        f.close()
                        if len(l):
                            for ll in l:
                                ll = string.split(string.strip(ll))
                                if ll[0][0:7]=='hw.ncpu':
                                    ncpu = int(ll[-1])
                                elif ll[0][0:10]=='hw.physmem':
                                    mem = int(ll[-1:][0])
                                    if mem>1000000000: # Gig or more
                                        cmd.set("hash_max",130)
                                    elif mem>500000000:
                                        cmd.set("hash_max",100)
                                    elif mem<256000000:
                                        cmd.set("hash_max",70)
                elif sys.platform[0:5]=='linux':
                    f=os.popen(
         "egrep -c '^processor[^A-Za-z0-9:]*: [0-9]' /proc/cpuinfo")
                    l=f.readlines()
                    f.close()
                    ncpu = int(l[0])
                elif sys.platform[0:4]=='irix':
                    f=os.popen("hinv | grep IP | grep Processor | grep HZ")
                    l=f.readlines()
                    f.close()
                    ncpu=int(string.split(string.strip(l[0]))[0])
                if ncpu>1:
                     cmd.set("max_threads",ncpu)
                     if invocation.options.show_splash:  
                          print " Detected %d CPU cores."%ncpu,
                          print " Enabled multithreaded rendering."
            except:
                pass
            cmd.reinitialize("store") # store our adapted state as default

        def launch_gui(self):
            try:
                if sys.platform=='darwin':
                    poll=1
                else:
                    poll=0
                skin = self.invocation.options.skin
                if self.invocation.options.external_gui==1:
                    __import__(self.invocation.options.gui)
                    sys.modules[self.invocation.options.gui].__init__(self,poll,skin)
                elif self.invocation.options.external_gui==3:
                    if not os.environ.has_key('DISPLAY'):
                        os.environ['DISPLAY']=':0.0'
                    os.environ['TCL_LIBRARY']=os.environ['PYMOL_PATH']+"/ext/lib/tcl8.4"
                    os.environ['TK_LIBRARY']=os.environ['PYMOL_PATH']+"/ext/lib/tk8.4"
                    __import__(self.invocation.options.gui)
                    sys.modules[self.invocation.options.gui].__init__(self,poll,skin)

            # -- Greg Landrum's RPC stuff
                if self.invocation.options.rpcServer:
                    import rpc
                    rpc.launch_XMLRPC()
            # --
            except:
                traceback.print_exc()

    def prime_pymol():
        global glutThread
        try:
            glutThread
        except NameError:
            glutThread = thread.get_ident()
        pymol_launch = 0 # never do this again : )
        if (sys.platform=='darwin') and (invocation.options.external_gui==1):
            import os
            xdpyinfo = "/usr/X11R6/bin/xdpyinfo"
            if os.path.exists(xdpyinfo):
                if os.system(xdpyinfo+" >/dev/null 2>&1"):
                    os.system("/usr/bin/open -a X11") # launch X11 (if needed)
            else:
                os.system("/usr/bin/open -a X11") # launch X11 (if needed)
                
    def start_pymol(block_input_hook=0):
        prime_pymol()
        _COb = _cmd._get_global_C_object()        
        _cmd.runpymol(_COb,block_input_hook) # only returns if we are running pretend GLUT
#      from pymol.embed import wxpymol # never returns

    import _cmd
    import cmd

    global _COb

    def thread_launch(pa):
        from pymol import invocation
        invocation.parse_args(pa)
        start_pymol(1)
    
    if pymol_launch == 1: # standard launch (consume main thread)
        cmd._COb = _cmd._get_global_C_object()
        if __name__=='pymol':
            _COb = cmd._COb
        else:
            pymol._COb = cmd._COb
        from pymol import invocation
        invocation.parse_args(pymol_argv)            
        start_pymol(0)

    elif pymol_launch == 2: # threaded launch (create new thread)
        cmd._COb = _cmd._get_global_C_object()
        if __name__=='pymol':
            _COb = cmd._COb
        else:
            pymol._COb = cmd._COb
        # don't do anything else yet...wait for finish_launching() call
        
    elif pymol_launch == 4: # monolithic (embedded) launch
        cmd._COb = _cmd._get_global_C_object()
        if __name__=='pymol':
            _COb = cmd._COb
        else:
            pymol._COb = cmd._COb
        from pymol import invocation
        invocation.parse_args(pymol_argv)            
        prime_pymol()
        # count on host process to actually start PyMOL

    if os.environ.has_key('DISPLAY'): # get X-window support
        from xwin import *

    def finish_launching(args=None):
        if args == None:
            args = pymol_argv+["-K"] # keep PyMOL thread alive
        else:
            args = list(args)
        if pymol_launch == 2: # spawn thread -- 'import pymol'
            global glutThreadObject
            cmd.reaper = threading.currentThread()
            glutThreadObject = threading.Thread(target=thread_launch,
              args=(args,)) 
            glutThreadObject.start()
        _COb = _cmd._get_global_C_object()
        e=threading.Event()
        import pymol # wait for import to complete
        while not _cmd.ready(_COb): # wait for the C library to initialize
            e.wait(0.01)
        while not hasattr(pymol,'xray'): # make sure symmetry module has time to start...
            e.wait(0.01)
            

