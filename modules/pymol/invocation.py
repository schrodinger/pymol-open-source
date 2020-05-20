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

# invocation.py
#
# This module unifies argument handling for embedded and modular PyMOL
#

helptext1 = '''Copyright (C) Schrodinger, LLC

Usage: pymol [OPTIONS]... [FILES]... [-- CUSTOM SCRIPT ARGUMENTS]

Options

  --help    display this help and exit
  --version display PyMOL version and exit
  --gldebug use glDebugMessageCallback for GL debugging
  --testing run pymol testing
  --diagnostics dump system diagnostics

  -1        config_mouse one_button
  -2        config_mouse two_button
  -a N      alias for -A
  -A N      application configuration:
    -A1     simple viewer window          (-qxiF -X 68 -Y 100)
    -A3     internal GUI only, no splash  (-qx -X 68 -Y 100)
    -A4     used by PYMOLVIEWER           (-X 68 -Y 100)
    -A5     helper application            (-QxiICUF -X 68 -Y 100)
    -A6     full screen presentation      (-qxieICUPF)
  -b[N]     benchmark wizard
  -B        (DEPRECATED)
  -c        launch in command-line only mode for batch processing
  -C        don't terminate on Ctrl-C
  -d cmd    execute PyMOL command
  -D N      defer_builds_mode=N
  -e        full screen
  -E N      multisampling (GL_MULTISAMPLE_ARB)
  -f N      internal_feedback=N
  -F        internal_feedback=0
  -g file   save image (png) or movie (mpg)
  -G        game mode (DEPRECATED)
  -h        generic helper application (no controls, no feedback)
  -H N      window height in pixels
  -i        internal_gui=0
  -I        auto_reinitialize=1 (Mac only)
  -j        side-by-side stereo (stereo_mode=4)
  -J        cd to user's home directory
  -k        don't load pymolrc or plugins
  -K        keep alive: when running without a GUI, don't quit after the input
            is exhausted
  -l file   run python script in thread (spawn)
  -L file   load file after everything else (only if something was loaded before)
  -m        INTERNAL - do not use (mac external GUI)
  -M        force mono
  -n        INTERNAL - do not use (incentive_product=1)
  -N name   UNSUPPORTED - external gui type (pmg_qt or pmg_tk) (same as -w)
  -o        disable security protections
  -O N      sphere_mode=N
  -p        read commands from STDIN
  -P        handle scenes as if the session were opened in presentation mode
  -q        supress startup message
  -Q        quiet, suppress all text output
  -r file   run python script
  -R        launch RPC Server
  -s file   log to file
  -S        force stereo
  -t N      stereo_mode=N
  -T name   UNSUPPORTED - Tcl/Tk GUI skin
  -u file   resume log file (execute existing content and append new log output)
  -U        UNSUPPORTED reuse the helper application
  -v        use openvr stub instead of a real hardware
  -V N      external GUI window height in pixels
  -w name   UNSUPPORTED - external gui type (pmg_qt or pmg_tk) (same as -N)
  -W N      window width in pixels
  -x        no external gui
  -X N      window x position on screen
  -y        exit on error
  -Y N      window y position on screen
  -z N      window_visible=N
  -Z N      zoom_mode=N

File Extensions

  pdb,sdf,...     molecular structure files
  ccp4,dx,...     map files

  py,pym,pyc      python script
  pml             PyMOL command script

  p5m             implies -A5 (PDB File)
  psw             implies -A6 (PyMOL Show File)
  pwg             PyMOL web GUI

Active "pymolrc" Files
'''

helptext2 = '''
Mail bug reports to https://lists.sourceforge.net/lists/listinfo/pymol-users
'''

if True:

    import copy
    import re
    import os
    import glob
    import sys
    import traceback

    pymolrc_pat1 = '.pymolrc*'
    pymolrc_pat2 = 'pymolrc*'

    ros_pat = 'run_on_startup*'

    class generic:
        pass

    global_options = generic();

    options = global_options

    options.deferred = []
    options.no_gui = 0
    options.internal_gui = 1
    options.internal_feedback = 1
    options.external_gui = 1
    options.force_stereo = -1 if sys.platform == 'darwin' else 0
    options.game_mode = 0
    options.gui = 'pmg_qt'
    options.skin = 'normal'
    options.show_splash = 1
    options.read_stdin = 0
    options.win_x = 640
    options.win_y = 480
    options.win_xy_set = False
    options.win_px = 4
    options.sigint_handler = 1 # terminate on Ctrl-C?
    options.reuse_helper = 0
    options.auto_reinitialize = 0
    options.keep_thread_alive = 0
    options.after_load_script = ""
    options.quiet = 0
    options.multisample = 0
    options.incentive_product = 0
    options.window_visible = 1
    options.presentation = 0
    options.defer_builds_mode = 0
    options.full_screen = 0
    options.sphere_mode = -1
    options.stereo_capable = 0
    options.stereo_mode = 0
    options.zoom_mode = -1
    options.no_quit = 0
    options.plugins = 2
    options.exit_on_error = 0
    options.pymolrc = None
    options.no_spacenav = 0
    options.launch_status = 0
    options.gldebug = 0
    options.testing = 0
    options.openvr_stub = False

    options.win_py = { 'irix':240,
                       'darwin': 214, # hmm...need to set to 192 for Leopard?...
                       'linux2': 220,
                       'win32' : 230}.get(sys.platform,200)

    options.ext_y = 168 # external gui height (eg. for Tcl/Tk top bar)

    options.blue_line = 0

    # Greg Landrum
    options.rpcServer = 0
    # end
    options.security = 1

    script_re = re.compile(r"pymolrc$|\.pml$|\.PML$|\.p1m$|\.P1M$")
    py_re = re.compile(r"\.py$|\.pym$|\.PY$|\.PYM$")

    def get_pwg_options(filename):
        for line in open(filename, 'r'):
            a = line.split()
            if not a or a[0].startswith('#'):
                continue
            if a[0].lower() == 'options':
                return a[1:]
        return []

    def get_personal_folder():
        if sys.platform.startswith('win'):
            try:
                import winreg
                with winreg.OpenKey(winreg.HKEY_CURRENT_USER,
                        r'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders') as key:
                    return winreg.QueryValueEx(key, "Personal")[0]
            except:
                print(' Warning: failed to query "My Documents" from registry')
        return os.path.expanduser('~')

    def get_user_config():
        for d in [os.getcwd(), '$HOME', '$HOMEDRIVE$HOMEPATH', '$PYMOL_PATH']:
            d = os.path.expandvars(d)
            for pat in [pymolrc_pat1, pymolrc_pat2]:
                lst = glob.glob(d + os.sep + pat)
                if lst:
                    break
            if lst:
                break
        # global run_on_startup script (not overridden by pymolrc files, but is disabled by "-k")
        if "PYMOL_PATH" in os.environ:
            first = glob.glob(os.environ['PYMOL_PATH']+"/"+ros_pat)
        else:
            first = []
        second = []
        for a in lst:
            if py_re.search(a):
                first.append(a) # preceeding "_ " cloaks
            elif script_re.search(a):
                second.append(a) # preceeding "_ " cloaks

        first.sort()
        second.sort()
        return first+second

    def parse_args(argv, _pymol=None, options=None, restricted=0):
        if not restricted:
            global _argv
            _argv = copy.deepcopy(argv) # pymol.invocation._argv
            global global_options
            if options is None:
                if _pymol is None:
                    options = global_options
                else:
                    options = _pymol.invocation.options
        av = copy.deepcopy(argv)
        av = av[1:] # throw out the executable path
        av.reverse()
        once_dict = {}
        options.deferred = []
        final_actions = []
        loaded_something = 0
        python_script = None
        # append user settings file as an option
        pymolrc = get_user_config()
        while 1:
            if not len(av):
                break
            a = av.pop()
            a = re.sub(r'''^"|"$|^'|'$''','',a) # strip extra quotes
            if a[0:1]=='-':
                if (a[1:2]=='-'):
                    if a in ('--version', '--help'):
                        import pymol
                        print(pymol.get_version_message())
                        if a == '--help':
                            print(helptext1)
                            if pymolrc:
                                for filename in pymolrc:
                                    print('  ' + filename)
                            else:
                                print('  (no pymolrc file found)')
                            print(helptext2)
                        sys.exit()
                    elif a == "--retina":
                        print("Warning: --retina option has been removed")
                    elif a == "--nospnav":
                        print(' Warning: --nospnav not available in Open-Source PyMOL')
                    elif a == "--gldebug":
                        options.gldebug = 1
                    elif a == "--testing":
                        options.testing = 1
                    elif a == "--diagnostics":
                        # same as: -cd diagnostics
                        options.no_gui=1
                        options.deferred.append("_do_diagnostics")
                    else:
                        # double hypen signals end of PyMOL arguments
                        if python_script is None:
                            python_script = argv[0]
                        rev_av = copy.deepcopy(av)
                        rev_av.reverse()
                        if len(a)>2:
                            sys.argv = [python_script] + [a] + rev_av
                        else:
                            sys.argv = [python_script] + rev_av
                        break
                    continue
                if ("A" in a) or ("a" in a): # application configuration
                    new_args = []
                    # ====== mode 1 - simple viewer window ======
                    if a[2:3] == "1":
                        if 'A1' not in once_dict:
                            once_dict['A1'] = 1
                            new_args = ["-qxiF",
                                "-X","68",
                                "-Y","100",
                                ]
                    # ====== mode 2 - not available -- clashes with -2 =======
                    # ====== mode 3 - internal GUI only no splash ======
                    if a[2:3] == "3":
                        if 'A3' not in once_dict:
                            once_dict['A3'] = 1
                            new_args = ["-qx",
                                "-X","68",
                                "-Y","100",
                                ]
                    # ====== mode 4 - internal GUI only with splash ======
                    if a[2:3] == "4": # used by PYMOLVIEWER
                        if 'A4' not in once_dict:
                            once_dict['A4'] = 1
                            new_args = [
                                "-X","68",
                                "-Y","100",
                                ]
                    # ====== mode 5 - mode 5 helper application ======
                    if a[2:3] == "5":
                        if 'A5' not in once_dict:
                            once_dict['A5'] = 1
                            new_args = ["-QxiICUF",
                                "-X","68",
                                "-Y","100",
                                ]
                    # ====== mode 6 - mode 6 presentation (no GUI) ======
                    if a[2:3] == "6":
                        if 'A6' not in once_dict:
                            once_dict['A6'] = 1
                            new_args = ["-qxieICUPF",
                                ]
                    # ===============================================
                    new_args.reverse()
                    av = av + new_args
                if "1" in a[1:2]:
                    options.deferred.append("_do__ config_mouse one_button")
                if "2" in a[1:2]:
                    options.deferred.append("_do__ config_mouse two_button")
                if "q" in a:
                    options.show_splash = 0
                if "i" in a:
                    options.internal_gui = 0
                if "f" in a:
                    options.internal_feedback = int(av.pop())
                if "F" in a:
                    options.internal_feedback = 0
                if "B" in a:
                    options.blue_line = 1
                if "E" in a:
                    options.multisample = int(av.pop())
                if "P" in a:
                    options.presentation = 1
                if "W" in a:
                    options.win_x = int(av.pop())
                    options.win_xy_set = True
                if "H" in a:
                    options.win_y = int(av.pop())
                    options.win_xy_set = True
                if "X" in a:
                    options.win_px = int(av.pop())
                if "y" in a:
                    options.exit_on_error = 1
                if "Y" in a:
                    options.win_py = int(av.pop())
                if "D" in a:
                    options.defer_builds_mode = int(av.pop())
                if "v" in a:
                    options.openvr_stub = True
                if "V" in a:
                    options.ext_y = int(av.pop())
                if "N" in a: # external gui name...
                    options.gui = av.pop()
                if "x" in a:
                    options.external_gui = 0
                if "n" in a:
                    options.incentive_product = 1
                if "t" in a: # type of stereo to use
                    options.stereo_mode = int(av.pop())
                if "T" in a: # what skin to use?
                    options.skin = str(av.pop())
                if "w" in a: # what gui to use
                    options.gui = str(av.pop())
                if "O" in a:
                    options.sphere_mode = int(av.pop())
                if "z" in a:
                    options.window_visible = 0
                if "Z" in a:
                    options.zoom_mode = int(av.pop())
                    if options.zoom_mode==5:
                        final_actions.append("_do__ zoom")
                if not restricted:
                    if "c" in a:
                        options.no_gui=1
                        options.external_gui=0
                    if "m" in a: # mac external GUI
                        if options.external_gui == 2:
                            options.external_gui = 3
                            if options.win_py == 184: # mac external GUI default
                                options.win_py = 216
                        else:
                            options.external_gui = 2
                            options.win_py = 184

                    if "e" in a:
                        options.full_screen = 1
                    if "G" in a: # Game mode (reqd for Mac stereo)
                        options.game_mode = 1
                        options.win_x = 1024
                        options.win_y = 768
                    if "S" in a: # Force stereo context on stereo-capable hardware
                        options.force_stereo = 1
                        if options.stereo_mode == 0:
                            options.stereo_mode = 1  # quadbuffer
                        if sys.platform=='darwin':
                            options.deferred.append(
                              "_do__ set stereo_double_pump_mono,1,quiet=1")
                    if "M" in a: # Force mono on stereo hardware (all)
                        options.force_stereo = -1
                    if "j" in a: # Geowall: two side-by-side images
                        options.stereo_mode = 4
                        options.deferred.append("_do__ stereo on")
                    if ("d" in a):
                        options.deferred.append(
                            "_do_" + av.pop().replace('%',' '))
                    if ("J" in a): # cd to user's home directory on startup (if possible)
                        path = get_personal_folder()
                        try:
                            # immediatly chdir (was: options.deferred.append(...))
                            os.chdir(path)
                            # clear PYMOL_WD, which may be set by MacPyMOL
                            os.environ.pop('PYMOL_WD', None)
                        except OSError:
                            print(" Error: could not chdir to", repr(path))
                    if ("l" in a):
                        options.deferred.append("_do_spawn %s"%av.pop())
                    if ("r" in a):
                        options.deferred.append("_do_run %s,main"%av.pop())
                    if ("u" in a):
                        options.deferred.append("_do_resume %s"%av.pop())
                    if ("s" in a):
                        options.deferred.append("_do_log_open %s"%av.pop())
                    if ("o" in a):
                        options.security = 0
                    if ("R" in a):
                        options.rpcServer = 1
                    if ("g" in a):
                        filename = av.pop()
                        if '.png' in filename:
                            options.deferred.append("_do__ cmd.png('''%s''')"%filename)
                        elif '.mpg' in filename:
                            options.deferred.append("_do__ movie.produce('''%s''')"%filename)
                    if ("C" in a):
                        options.sigint_handler = 0
                    if ("L" in a):
                        options.after_load_script = av.pop()
                    if ("b" in a): # CPU benchmark
                        options.deferred.append("_do__ feedback disable,all,everything")
                        options.deferred.append("_do__ feedback enable,python,output")
                        options.deferred.append("_do_wizard benchmark")
                        if a[2:]=='':
                            options.deferred.append("_do__ cmd.get_wizard().run_cpu()")
                        if a[2:]=='0':
                            options.deferred.append("_do__ cmd.get_wizard().ray_trace0()")
                        if a[2:]=='1':
                            options.deferred.append("_do__ cmd.get_wizard().ray_trace1()")
                        if a[2:]=='2':
                            options.deferred.append("_do__ cmd.get_wizard().ray_trace2()")

                    if "p" in a:
                        options.read_stdin = 1
                    if "K" in a:
                        options.keep_thread_alive = 1
                if "k" in a: # suppress reading of .pymolrc and related files
                    pymolrc = None
                    options.plugins = 0
                if "U" in a: #
                    options.reuse_helper = 1
                if "Q" in a:
                    options.quiet = 1
                    options.show_splash = 0
                if "I" in a:
                    options.auto_reinitialize = 1
                if "h" in a: # generic helper application
                    options.internal_gui = 0
                    options.external_gui = 0
                    options.internal_feedback = 0
                    options.show_splash = 1
            elif a in ('+1', '+2', '+3', '+4'):
                print('ignoring PyMOLWin.exe argument', a)
            elif not restricted:
                suffix = a[-4:].lower().split('.')[-1]
                if suffix == "p5m":
                    # mode 5 helper application
                    av.append("-A5")
                elif suffix == "psw":
                    # presentation mode
                    av.append("-A6")
                elif suffix in [ 'pym' ,'py', 'pyc' ]:
                    python_script = a
                elif suffix in [ 'pwg' ]:
                    try:
                        pwg_options = get_pwg_options(a)
                        if pwg_options:
                            parse_args(['pymol'] + pwg_options, _pymol, options, 1)
                    except:
                        traceback.print_exc()
                options.deferred.append(a)
                loaded_something = 1
        if pymolrc is not None:
            options.deferred = [('_do__ @' + a) if script_re.search(a) else a
                    for a in pymolrc] + options.deferred
            options.pymolrc = pymolrc
        if options.rpcServer:
            options.deferred.append('_do__ /import pymol.rpc;pymol.rpc.launch_XMLRPC()')
        if options.plugins == 1:
            # Load plugins independent of PMGApp (will not add menu items)
            options.deferred.append('_do__ /import pymol.plugins;pymol.plugins.initialize(-1)')
        if loaded_something and (options.after_load_script!=""):
            options.deferred.append(options.after_load_script)
        options.deferred.extend(final_actions)
        if options.show_splash and not options.no_gui and not restricted:
            options.deferred.insert(0,"_do__ cmd.splash(1)")
        if options.full_screen:
            options.deferred.append("_do__ full_screen on")
