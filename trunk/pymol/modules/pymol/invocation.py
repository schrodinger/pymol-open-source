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
    
    pattern1 = '.pymolrc*'
    pattern2 = 'pymolrc*'

    class generic:
        pass

    options = generic();

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
    options.passive_stereo = 0
    options.zoom_mode = -1
    
    if sys.platform[0:4] == 'irix':
        options.win_py = 240
    elif sys.platform == 'darwin':
        options.win_py = 214
    elif sys.platform != 'win32':
        options.win_py = 200
    else:
        options.win_py = 230

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
        lst = glob.glob(pattern1)
        if not len(lst): # unix
            if os.environ.has_key("HOME"):
                lst = glob.glob(os.environ['HOME']+"/"+pattern1)
        if not len(lst): # unix
            if os.environ.has_key("HOME"):
                lst = glob.glob(os.environ['HOME']+"/"+pattern2)
        if not len(lst): # win32
            if os.environ.has_key("HOMEPATH") and os.environ.has_key("HOMEDRIVE"):
                lst = glob.glob(os.environ['HOMEDRIVE']+os.environ['HOMEPATH']+"/"+pattern1)
        if not len(lst): # win32
            if os.environ.has_key("HOMEPATH") and os.environ.has_key("HOMEDRIVE"):
                lst = glob.glob(os.environ['HOMEDRIVE']+os.environ['HOMEPATH']+"/"+pattern2)
        if not len(lst): # all
            if os.environ.has_key("PYMOL_PATH"):
                lst = glob.glob(os.environ['PYMOL_PATH']+"/"+pattern1)
        if not len(lst): # all
            if os.environ.has_key("PYMOL_PATH"):
                lst = glob.glob(os.environ['PYMOL_PATH']+"/"+pattern2)

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

    def parse_args(argv):
        av = copy.deepcopy(argv)
        av = av[1:] # throw out the executable path
        av.reverse()
        global options
        once_dict = {}
        options.deferred = []
        final_actions = []
        loaded_something = 0
        # append user settings file as an option
        options.deferred.extend(get_user_config())
        while 1:
            if not len(av):
                break
            a = av.pop()
            a = re.sub(r'''^"|"$|^'|'$''','',a) # strip extra quotes
            if a[0:1]=='-':
                if a[1:2]=='-':
                    break # double hypen signals end of PyMOL arguments
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
                if "c" in a:
                    options.no_gui=1
                    options.external_gui=0
                if "2" in a:
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
                if "m" in a: # mac external GUI
                    if options.external_gui == 2:
                        options.external_gui = 3
                        if options.win_py == 184: # mac external GUI default
                            options.win_py = 216
                    else:
                        options.external_gui = 2
                        options.win_py = 184 
                if "n" in a:
                    options.incentive_product = 1
                if "t" in a:
                    options.gui = 'pmg_tk'
                if "T" in a: # what skin to use? 
                    options.skin = str(av.pop())
                if "w" in a:
                    options.gui = 'pmg_wx'
                if "O" in a:
                    options.sphere_mode = int(av.pop())
                if "z" in a:
                    options.window_visible = 0
                if "Z" in a:
                    options.zoom_mode = int(av.pop())
                    if options.zoom_mode==5:
                        final_actions.append("_do__ zoom")
                if "d" in a:
                    options.deferred.append(
                        "_do_%s"%string.replace(av.pop(),'%',' '))
                if "e" in a:
                    options.full_screen = 1
                    options.deferred.append("_do__ full_screen on")
                if "G" in a: # Game mode (reqd for Mac stereo)
                    options.game_mode = 1
                    options.win_x = 1024
                    options.win_y = 768
                if "S" in a: # Force stereo on stereo hardware (OSX only)
                    options.force_stereo = 1
                    options.deferred.append("_do__ stereo on")
                    if sys.platform=='darwin': 
                        options.deferred.append(
                          "_do__ set stereo_double_pump_mono,1,quiet=1")
                if "M" in a: # Force mono on stereo hardware (all)
                    options.force_stereo = -1
                if "j" in a: # Geowall: two side-by-side images
                    options.passive_stereo = 1
                    options.deferred.append("_do__ stereo on")
                if "l" in a:
                    options.deferred.append("_do_spawn %s"%av.pop())
                if "r" in a:
                    options.deferred.append("_do_run %s,main"%av.pop())
                if "u" in a:
                    options.deferred.append("_do_resume %s"%av.pop())
                if "s" in a:
                    options.deferred.append("_do_log_open %s"%av.pop())
                if "p" in a:
                    options.read_stdin = 1
                if "o" in a:
                    options.security = 0
                if "R" in a:
                    options.rpcServer = 1
                if "K" in a:
                    options.keep_thread_alive = 1
                if "g" in a:
                    options.deferred.append("_do_png %s"%av.pop())
                if "C" in a:
                    options.sigint_handler = 0
                if "U" in a:
                    options.reuse_helper = 1
                if "Q" in a:
                    options.quiet = 1
                    options.show_splash = 0
                if "I" in a:
                    options.auto_reinitialize = 1
                if "L" in a:
                    options.after_load_script = av.pop()
                if "h" in a: # generic helper application
                    options.internal_gui = 0
                    options.external_gui = 0
                    options.internal_feedback = 0
                    options.show_splash = 1
                if "b" in a: # CPU benchmark
                    options.deferred.append("_do__ feedback disable,all,everything")
                    options.deferred.append("_do__ feedback enable,python,output")
                    options.deferred.append("_do_ wizard benchmark")
                    if a[2:]=='':
                        options.deferred.append("_do_ cmd.get_wizard().run_cpu()")
                    if a[2:]=='0':
                        options.deferred.append("_do_ cmd.get_wizard().ray_trace0()")
                    if a[2:]=='1':
                        options.deferred.append("_do_ cmd.get_wizard().ray_trace1()")
                    if a[2:]=='2':
                        options.deferred.append("_do_ cmd.get_wizard().ray_trace2()")
                        
            else: 
                if a[-4:] in (".p5m",".P5M"):
                    # mode 5 helper application 
                    av.append("-A5")
                if a[-4:] in (".psw",".PSW"):
                    # presentation mode
                    av.append("-A6")               
                options.deferred.append(a)
                loaded_something = 1
        if loaded_something and (options.after_load_script!=""):
            options.deferred.append(options.after_load_script)
        options.deferred.extend(final_actions)
        if options.show_splash and not options.no_gui:
            options.deferred.insert(0,"_do__ cmd.splash(1)")
        
