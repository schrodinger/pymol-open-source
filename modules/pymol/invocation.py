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

if __name__=='pymol.invocation':

    import copy
    import re
    import os
    import glob
    import string
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
    options.force_stereo = 0
    options.game_mode = 0
    options.gui = 'pmg_tk'
    options.skin = 'normal'
    options.show_splash = 1
    options.read_stdin = 0
    options.win_x = 640
    options.win_y = 480
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
    pyc_re = re.compile(r"\.pyc$|\.PYC$") # not yet used

    def get_user_config():
        # current working directory
        lst = glob.glob(pymolrc_pat1)
        # users home directory
        if not len(lst): # unix
            if os.environ.has_key("HOME"):
                lst = glob.glob(os.environ['HOME']+"/"+pymolrc_pat1)
        if not len(lst): # unix
            if os.environ.has_key("HOME"):
                lst = glob.glob(os.environ['HOME']+"/"+pymolrc_pat2)
        if not len(lst): # win32
            if os.environ.has_key("HOMEPATH") and os.environ.has_key("HOMEDRIVE"):
                lst = glob.glob(os.environ['HOMEDRIVE']+os.environ['HOMEPATH']+"/"+pymolrc_pat1)
        if not len(lst): # win32
            if os.environ.has_key("HOMEPATH") and os.environ.has_key("HOMEDRIVE"):
                lst = glob.glob(os.environ['HOMEDRIVE']+os.environ['HOMEPATH']+"/"+pymolrc_pat2)
        # installation folder (if known)
        if not len(lst): # all
            if os.environ.has_key("PYMOL_PATH"):
                lst = glob.glob(os.environ['PYMOL_PATH']+"/"+pymolrc_pat1)
        if not len(lst): # all
            if os.environ.has_key("PYMOL_PATH"):
                lst = glob.glob(os.environ['PYMOL_PATH']+"/"+pymolrc_pat2)
        # global run_on_startup script (not overridden by pymolrc files, but is disabled by "-k")
        if os.environ.has_key("PYMOL_PATH"):
            first = glob.glob(os.environ['PYMOL_PATH']+"/"+ros_pat)
        else:
            first = []
        second = []
        for a in lst:
            if py_re.search(a):
                first.append("_do__ run "+a) # preceeding "_ " cloaks
            elif script_re.search(a):
                second.append("_do__ @"+a) # preceeding "_ " cloaks 
    #      elif pyc_re.search(a): # ignore compiled versions for now
    #         first.append("_do__ run "+a) # preceeding "_ " cloaks

        first.sort()
        second.sort()
        return first+second

    def parse_args(argv, _pymol=None, options=None, restricted=0):
        if not restricted:
            global _argv
            _argv = copy.deepcopy(argv) # pymol.invocation._argv
            global global_options
            if options == None:
                if _pymol==None:
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
                    # double hypen signals end of PyMOL arguments
                    if python_script == None:
                        python_script = argv[0]
                    rev_av = copy.deepcopy(av)
                    rev_av.reverse()
                    if len(a)>2:
                        sys.argv = [python_script] + [a] + rev_av
                    else:
                        sys.argv = [python_script] + rev_av
                    break
                if ("A" in a) or ("a" in a): # application configuration
                    new_args = []
                    # ====== mode 1 - simple viewer window ======
                    if a[2:3] == "1": 
                        if not once_dict.has_key('A1'):
                            once_dict['A1'] = 1
                            new_args = ["-qxiF",
                                "-X","68",
                                "-Y","100",
                                ]
                    # ====== mode 2 - not available -- clashes with -2 =======
                    # ====== mode 3 - internal GUI only no splash ======
                    if a[2:3] == "3": 
                        if not once_dict.has_key('A3'):
                            once_dict['A3'] = 1
                            new_args = ["-qx",
                                "-X","68",
                                "-Y","100",
                                ]
                    # ====== mode 4 - internal GUI only with splash ======
                    if a[2:3] == "4": # used by PYMOLVIEWER
                        if not once_dict.has_key('A4'):
                            once_dict['A4'] = 1
                            new_args = [
                                "-X","68",
                                "-Y","100",
                                ]
                    # ====== mode 5 - mode 5 helper application ======
                    if a[2:3] == "5": 
                        if not once_dict.has_key('A5'):
                            once_dict['A5'] = 1
                            new_args = ["-QxiICUF",
                                "-X","68",
                                "-Y","100",
                                ]
                    # ====== mode 6 - mode 6 presentation (no GUI) ======
                    if a[2:3] == "6": 
                        if not once_dict.has_key('A6'):
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
                if "H" in a:
                    options.win_y = int(av.pop())
                if "X" in a:
                    options.win_px = int(av.pop())
                if "Y" in a:
                    options.win_py = int(av.pop())
                if "D" in a:
                    options.defer_builds_mode = int(av.pop())
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
                        options.deferred.append("_do__ full_screen on")
                    if "G" in a: # Game mode (reqd for Mac stereo)
                        options.game_mode = 1
                        options.win_x = 1024
                        options.win_y = 768
                    if "S" in a: # Force stereo context on stereo-capable hardware
                        options.force_stereo = 1
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
                            "_do_%s"%string.replace(av.pop(),'%',' '))
                    if ("J" in a): # cd to user's home directory on startup (if possible)
                        if sys.platform == 'win32':
                            if os.environ.has_key("HOMEDRIVE") and os.environ.has_key("HOMEPATH"):
                                path = os.environ["HOMEDRIVE"] + os.environ["HOMEPATH"]
                                if os.path.isdir(path) and os.path.exists(path):
                                    my_docs = os.path.join(path,"Documents") # for VISTA compatibility
                                    if os.path.isdir(my_docs): # start in Documents (if exists)
                                        path = my_docs
                                    else:
                                        my_docs = os.path.join(path,"My Documents")                                    
                                        if os.path.isdir(my_docs): # start in My Documents (if exists)
                                            path = my_docs
                                    options.deferred.append("_do__ cmd.cd('''%s''',complain=0)"%string.replace(path,"\\","\\\\"))
                        elif os.environ.has_key("HOME"):
                            path = os.environ["HOME"]
                            if os.path.isdir(path):
                                options.deferred.append("_do__ cmd.cd('''%s''',complain=0)"%string.replace(path,"\\","\\\\"))
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
            elif not restricted:
                suffix = string.split(string.lower(a[-4:]),'.')[-1]
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
                        lines = open(a,'rb').readlines()
                        pseudo_argv = ["pymol"]
                        for line in lines:
                            line = line.strip()
                            if len(line) and line[0:1] != '#':
                                input = line.split(None,1)
                                if len(input) and input[0]!='#':
                                    keyword = input[0].lower()
                                    if keyword == 'options':
                                        if len(input)>1:
                                            pseudo_argv = ['pymol'] + input[1].split()
                                            parse_args(pseudo_argv, _pymol=_pymol,
                                                       options=options, restricted=1)
                    except:
                        traceback.print_exc()
                options.deferred.append(a)
                loaded_something = 1
        if pymolrc != None:
            options.deferred = pymolrc + options.deferred
        if loaded_something and (options.after_load_script!=""):
            options.deferred.append(options.after_load_script)
        options.deferred.extend(final_actions)
        if options.show_splash and not options.no_gui and not restricted:
            options.deferred.insert(0,"_do__ cmd.splash(1)")
        
