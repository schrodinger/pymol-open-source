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

import __main__
if __name__!='__main__':
    import invocation
    
# Global variable "__main__.pymol_launch" tracks how we're launching PyMOL:
#
# 0: old way, now obsolete (e.g. "python launch_pymol.py")
# 1: new way, spawn thread on import: (e.g. from "import pymol")
# 2: new way, consume main thread: (e.g. from "python pymol/__init__.py")
# 3: dry run -- just get PyMOL environment information
# 4: monolithic (embedded) PyMOL.  Prime, but don't start.
# 5: Python embedded launch from within the PyMOL API

if hasattr(__main__,'pymol_launch'):
    pymol_launch = __main__.pymol_launch
else:
    pymol_launch = 2 
    
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

        try_again = 0

        try:
            import pymol
        except ImportError:
            try_again = 1

        if try_again:
            try:
                sys.exc_clear()
                sys.path.insert(0,os.path.dirname(os.path.dirname(__file__)))
            except:
                pass
            import pymol
        
        pymol_launch = 1
    
    elif pymol_launch==-1:

        if hasattr(__main__,"pymol_argv"):
            pymol_argv = __main__.pymol_argv
        else:
            pymol_argv = [ "pymol", "-q" ]

    elif pymol_launch==2:

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
        
        # Create a temporary object "stored" in the PyMOL global namespace
        # for usage with evaluate based-commands such as alter

        class Scratch_Storage:
            pass

        stored = Scratch_Storage()

        # Create a permanent object in the PyMOL global namespace
        # that will be picked and unpickled along with the session

        class Session_Storage:
            pass

        session = Session_Storage()

        # This global will be non-None if logging is active
        # (global variable used for efficiency)

        _log_file = None

        # This global will be non-None if an external gui
        # exists. It mainly exists so that events which occur
        # in the Python thread can be handed off to the
        # external GUI thread through one or more FIFO Queues
        # (global variable used for efficiency)

        _ext_gui = None

        # lists of functions to call when saving and restoring pymol session objects
        # The entry 'None' represents the PyMOL C-API function call

        _session_save_tasks = [ None ]
        _session_restore_tasks = [ None ]

        # standard input reading thread

        _stdin_reader_thread = None
        
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

        # these locks are to be shared by all PyMOL instances within a
        # single Python interpeter
        
        lock_api = threading.RLock() # mutex for API 
        lock_api_c = threading.RLock() # mutex for C management of python threads
        lock_api_status = threading.RLock() # mutex for PyMOL status info
        lock_api_glut = threading.RLock() # mutex for avoiding GLUT
        
        def exec_str(self,s):
            try:
                exec s in self.__dict__, self.__dict__
            except StandardError:
                traceback.print_exc()
            return None

        def exec_deferred():
            if invocation.options.read_stdin:
                try:
                    _stdin_reader_thread = threading.Thread(target=cmd._parser.stdin_reader)
                    _stdin_reader_thread.setDaemon(1)
                    _stdin_reader_thread.start()
                except:
                    trackback.print_exc()
            try:
                if cmd.ready():
                    cmd.config_mouse(quiet=1)
                    for a in invocation.options.deferred:
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
                print "Error: Argument processing aborted due to exception."
            except:
                traceback.print_exc()

        def adapt_to_hardware():
            (vendor,renderer,version) = cmd.get_renderer()
            if vendor[0:6]=='NVIDIA':
                cmd.set('ribbon_smooth',0,quiet=1)
                if renderer[0:7]=='GeForce':
                    if invocation.options.show_splash:
                        print " Adapting to GeForce hardware."
                    cmd.set('line_width','2',quiet=1)
                elif renderer=='NVIDIA GPU OpenGL Engine':
                    if sys.platform=='darwin':
                        if invocation.options.show_splash:
                            print " Adapting to NVIDIA hardware on Mac."
                            cmd.set('line_smooth',0,quiet=1)
                            cmd.set('fog',0.9,quiet=1)
                elif renderer=='NVIDIA GeForce4 GPU OpenGL Engine':
                    if sys.platform=='darwin':
                        if invocation.options.show_splash:
                            cmd.set('stereo_double_pump_mono',1,quiet=1)

                elif renderer[0:6]=='Quadro':
                    if invocation.options.show_splash:
                        print " Adapting to Quadro hardware."
                    cmd.set("stereo_double_pump_mono","1",quiet=1)
                    cmd.set("line_width",1.4,quiet=1)

            elif vendor[0:4]=='Mesa':
                if renderer[0:18]=='Mesa GLX Indirect':
                    cmd.set('ribbon_smooth',0,quiet=1)

            elif vendor[0:3]=='ATI':
                cmd.set('ribbon_smooth',0,quiet=1)
                if renderer[0:17]=='FireGL2 / FireGL3':
                    if invocation.options.show_splash:
                        print " Adapting to FireGL hardware."
                    cmd.set('line_width','2',quiet=1)            
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
                    f=os.popen("egrep -c '^processor[^A-Za-z0-9:]*: [0-9]' /proc/cpuinfo")
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
                          print " Detected %d CPUs."%ncpu,
                          print " Enabled multithreaded rendering."
            except:
                pass

        # NEED SOME CONTRIBUTIONS HERE!

        def launch_gui():
            if sys.platform=='darwin':
                poll=1
       	    else:
                poll=0
            if invocation.options.external_gui==1:
                __import__(invocation.options.gui)
                sys.modules[invocation.options.gui].__init__(sys.modules['pymol'],poll)
            elif invocation.options.external_gui==3:
                os.environ['DISPLAY']=':0.0'
                os.environ['TCL_LIBRARY']=os.environ['PYMOL_PATH']+"/ext/lib/tcl8.4"
                os.environ['TK_LIBRARY']=os.environ['PYMOL_PATH']+"/ext/lib/tk8.4"
                __import__(invocation.options.gui)
                sys.modules[invocation.options.gui].__init__(sys.modules['pymol'],poll)

        # -- Greg Landrum's RPC stuff
            if invocation.options.rpcServer:
                import rpc
                rpc.launch_XMLRPC()
        # --

    def prime_pymol():
        global glutThread
        glutThread = thread.get_ident()
        pymol_launch = 0 # never do this again : )
        if (sys.platform=='darwin') and (invocation.options.external_gui==1):
            import os
            os.system("/usr/bin/open -a X11") # launch X11 if we're going to need it

    def start_pymol(block_input_hook=0):
        prime_pymol()
        _cmd.runpymol(block_input_hook) # only returns if we are running pretend GLUT
#      from pymol.embed import wxpymol # never returns

    import _cmd
    import cmd

    global _COb

    def thread_launch(pa):
        from pymol import invocation
        invocation.parse_args(pa)
        start_pymol(1)

    if pymol_launch==1: # standard launch (absorb main thread)
        cmd._COb = _cmd._get_global_C_object()
        if __name__=='pymol':
            _COb = cmd._COb
        else:
            pymol._COb = cmd._COb
        from pymol import invocation
        invocation.parse_args(pymol_argv)
        start_pymol(0)

    elif pymol_launch==2: # threaded launch (create new thread)
        cmd._COb = _cmd._get_global_C_object()
        if __name__=='pymol':
            _COb = cmd._COb
        else:
            pymol._COb = cmd._COb
        global glutThreadObject
        cmd.reaper = threading.currentThread()
        glutThreadObject = threading.Thread(target=thread_launch,
          args=(pymol_argv+["-K"],)) # keep PyMOL thread alive
                                              # even w/o GUI
        glutThreadObject.start()

    elif pymol_launch==4: # monolithic (embedded) launch
        cmd._COb = _cmd._get_global_C_object()
        if __name__=='pymol':
            _COb = cmd._COb
        else:
            pymol._COb = cmd._COb
        invocation.parse_args(pymol_argv)
        prime_pymol()
        # count on host process to actually start PyMOL

    if os.environ.has_key('DISPLAY'):
        from xwin import *

    def finish_launching():
        e=threading.Event()
        import pymol # wait for import to complete
        while not _cmd.ready(): # wait for the C library to initialize
            e.wait(0.01)
        while not hasattr(pymol,'xray'): # make sure symmetry module has time to start...
            e.wait(0.01)
            

