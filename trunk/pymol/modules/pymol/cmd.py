
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

# cmd.py 
# Python interface module for PyMol
#
# **This is the only module which should be/need be imported by 
# ** PyMol API Based Programs

# NEW CALL RETURN CONVENTIONS for _cmd.so C-layer
#

# (1) Calls into C (_cmd) should return results/status and print
#     errors and feedback (according to mask) BUT NEVER RAISE EXCEPTIONS
#     from within the C code itself.

# (2) Effective with version 0.99, standard Python return conventions
# apply, but haven't yet been fully implemented.  In summary:

#     Unless explicitly specified in the function:

#     ==> Success with no information should return None

#     ==> Failure should return a negative number as follows:
#        -1 = a general, unspecified failure
          
#     Upon an error, exceptions will be raised by the Python wrapper
#     layer if the "raise_exceptions" setting is on.

#     ==> Boolean queries should return 1 for true/yes and 0 for false/no.

#     ==> Count queries should return 0 or a positive number

# (3) If _cmd produces a specific return result, be sure to include an
#     error result as one of the possibilities outside the range of the
#     expected return value.  For example, a negative distance
#
# (4) cmd.py API wrappers can then raise exceptions and return values.
#
# NOTE: Output tweaking via the "quiet" parameter of API functions.
#
# Many PyMOL API functions have a "quiet" parameter which is used to
# customize output depending on when and where the routine was called.
#
# As defined, quiet should be 1.  Called from an external module, output
# to the console should be minimal.
#
# However, when a command is run through the parser (often manually) the
# user expects a little more feedback.  The parser will automatically
# set "quiet" to zero
#
# In rare cases, certain nonserious error or warning output should
# also be suppressed.  Set "quiet" to 2 for this behavior.

if __name__=='pymol.cmd':

    import traceback
    import sys

    try:
        
        import re
        import _cmd
        import string
        import thread
        import threading
        import types
        import pymol
        import os
        import parsing
        import __main__
        import time
        import urllib
        
        from shortcut import Shortcut

        from chempy import io

        #######################################################################
        # symbols for early export
        #######################################################################

        file_ext_re = re.compile(string.join([
            "\.pdb$|\.pdb1$|\.ent$|\.mol$|\.p5m$|",
            r"\.mmod$|\.mmd$|\.dat$|\.out$|\.mol2$|",
            r"\.xplor$|\.pkl$|\.sdf$|\.pqr|", 
            r"\.r3d$|\.xyz$|\.xyz_[0-9]*$|", 
            r"\.cc1$|\.cc2$|", # ChemDraw 3D
            r"\.cube$|", # Gaussian Cube
            r"\.dx$|", # DX files (APBS)
            r"\.pse$|", # PyMOL session (pickled dictionary)
            r"\.pmo$|", # Experimental molecular object format
            r"\.moe$|", # MOE (proprietary)
            r"\.mae$|", # MAE (proprietary)            
            r"\.ccp4$|", # CCP4
            r"\.top$|", # AMBER Topology
            r"\.trj$|", # AMBER Trajectory
            r"\.crd$|", # AMBER coordinate file
            r"\.rst$|", # AMBER restart
            r"\.cex$|", # CEX format (used by metaphorics)
            r"\.phi$|", # PHI format (delphi)
            r"\.fld$|", # FLD format (AVS)
            r"\.trj$|\.trr$|\.xtc$|\.gro$|\.g96$|\.dcd$|", # Trajectories
            r"\.o$|\.omap$|\.dsn6$|\.brix$|", # BRIX/O format
            r"\.grd$", # InsightII Grid format
            ],''), re.I)

        reaper = None
        safe_oname_re = re.compile(r"\ |\+|\(|\)|\||\&|\!|\,")  # quash reserved characters
        sanitize_list_re = re.compile(r"[^0-9\.\-\[\]\,]+")
        sanitize_alpha_list_re = re.compile(r"[^a-zA-Z0-9_\'\"\.\-\[\]\,]+")
        nt_hidden_path_re = re.compile(r"\$[\/\\]")
        quote_alpha_list_re = re.compile(
            r'''([\[\,]\s*)([a-zA-Z_][a-zA-Z0-9_\ ]*[a-zA-Z0-9_]*)(\s*[\,\]])''')
        def safe_list_eval(st):
            return eval(sanitize_list_re.sub('',st))

        def safe_alpha_list_eval(st):
            st = sanitize_alpha_list_re.sub('',st)
            st = quote_alpha_list_re.sub(r'\1"\2"\3',st) # need to do this twice
            st = quote_alpha_list_re.sub(r'\1"\2"\3',st)
            return eval(sanitize_alpha_list_re.sub('',st))

        QuietException = parsing.QuietException
        
        DEFAULT_ERROR = -1
        DEFAULT_SUCCESS = None
        
        #--------------------------------------------------------------------
        # shortcuts...

        toggle_dict = {'on':1,'off':0,'1':1,'0':0,'toggle':-1, '-1':-1}
        toggle_sc = Shortcut(toggle_dict.keys())

        stereo_dict = {'on':1,'off':0,'1':1,'0':0,'swap':-1,
                       'crosseye':2,'quadbuffer':3,
                       'walleye':4,'geowall':5,'sidebyside':6}
        
        stereo_sc = Shortcut(stereo_dict.keys())

        space_sc = Shortcut(['cmyk','rgb','pymol'])

        window_dict = { 'show' : 1, 'hide' : 0, 'position' : 2, 'size' : 3,
                        'box' : 4, 'maximize' : 5, 'fit' : 6, 'focus' : 7,
                        'defocus' : 8 }
        window_sc = Shortcut(window_dict.keys())

        repres = {
            'everything'    : -1,
            'sticks'        : 0,
            'spheres'       : 1,
            'surface'       : 2,
            'labels'        : 3,
            'nb_spheres'    : 4,
            'cartoon'       : 5,
            'ribbon'        : 6,
            'lines'         : 7,
            'mesh'          : 8,
            'dots'          : 9,
            'dashes'        :10,
            'nonbonded'     :11,
            'cell'          :12,
            'cgo'           :13,
            'callback'      :14,
            'extent'        :15,
            'slice'         :16,
            'angles'        :17,
            'dihedrals'     :18,         
        }
        repres_sc = Shortcut(repres.keys())

        boolean_dict = {
            'yes'           : 1,
            'no'            : 0,
            '1'             : 1,
            '0'             : 0,
            'on'            : 1,
            'off'           : 0
            }

        boolean_sc = Shortcut(boolean_dict.keys())

        palette_dict = {
            'rainbow_cycle'           : ('o',3,0  ,999), # perceptive rainbow
            'rainbow_cycle_rev'       : ('o',3,999,  0),      
            'rainbow'                 : ('o',3,107,893),
            'rainbow_rev'             : ('o',3,893,107),
            'rainbow2'                : ('s',3,167, 833), # cartesian rainbow 
            'rainbow2_rev'            : ('s',3,833,167),

            'gcbmry' : ('r',3,166,999),
            'yrmbcg' : ('r',3,999,166),

            'cbmr'   : ('r',3,166,833),
            'rmbc'   : ('r',3,833,166),      

            'green_yellow_red'        : ('s',3,500,833),
            'red_yellow_green'        : ('s',3,833,500),      

            'yellow_white_blue'       : ('w',3,  0, 83),
            'blue_white_yellow'       : ('w',3, 83,  0),      

            'blue_white_red'          : ('w',3, 83,167),
            'red_white_blue'          : ('w',3,167, 83),

            'red_white_green'         : ('w',3,167,250),
            'green_white_red'         : ('w',3,250,167),

            'green_white_magenta'     : ('w',3,250,333),
            'magenta_white_green'     : ('w',3,333,250),      

            'magenta_white_cyan'      : ('w',3,333,417),
            'cyan_white_magenta'      : ('w',3,417,333),

            'cyan_white_yellow'       : ('w',3,417,500),
            'yellow_cyan_white'       : ('w',3,500,417),

            'yellow_white_green'      : ('w',3,500,583),
            'green_white_yellow'      : ('w',3,583,500),

            'green_white_blue'        : ('w',3,583,667),
            'blue_white_green'        : ('w',3,667,583),      

            'blue_white_magenta'      : ('w',3,667,750),
            'magenta_white_blue'      : ('w',3,750,667),

            'magenta_white_yellow'    : ('w',3,750,833),
            'yellow_white_magenta'    : ('w',3,833,750),

            'yellow_white_red'        : ('w',3,833,917),
            'red_white_yellow'        : ('w',3,817,833),

            'red_white_cyan'          : ('w',3,916,999),
            'cyan_white_red'          : ('w',3,999,916),      

            'yellow_blue'       : ('c',3,  0, 83),
            'blue_yellow'       : ('c',3, 83,  0),      

            'blue_red'          : ('c',3, 83,167),
            'red_blue'          : ('c',3,167, 83),

            'red_green'         : ('c',3,167,250),
            'green_red'         : ('c',3,250,167),

            'green_magenta'     : ('c',3,250,333),
            'magenta_green'     : ('c',3,333,250),      

            'magenta_cyan'      : ('c',3,333,417),
            'cyan_magenta'      : ('c',3,417,333),

            'cyan_yellow'       : ('c',3,417,500),
            'yellow_cyan'       : ('c',3,500,417),

            'yellow_green'      : ('c',3,500,583),
            'green_yellow'      : ('c',3,583,500),

            'green_blue'        : ('c',3,583,667),
            'blue_green'        : ('c',3,667,583),      

            'blue_magenta'      : ('c',3,667,750),
            'magenta_blue'      : ('c',3,750,667),

            'magenta_yellow'    : ('c',3,750,833),
            'yellow_magenta'    : ('c',3,833,750),

            'yellow_red'        : ('c',3,833,917),
            'red_yellow'        : ('c',3,817,833),

            'red_cyan'          : ('c',3,916,999),
            'cyan_red'          : ('c',3,999,916),      
            }

        palette_sc = Shortcut(palette_dict.keys())

        def null_function():
            pass

        #--------------------------------------------------------------------
        # convenient type-checking

        def is_string(obj):
            return isinstance(obj,types.StringType)

        def is_list(obj):
            return isinstance(obj,types.ListType)

        def is_tuple(obj):
            return isinstance(obj,types.TupleType)

        def is_sequence(obj):
            return isinstance(obj,types.ListType) or isinstance(obj,types.TupleType)

        #-------------------------------------------------------------------
        # path expansion, including our fixes for Win32

        def _nt_expandvars(path): # allow for //share/folder$/file
            path = nt_hidden_path_re.sub(r"$$\\",path)
            return os.path.expandvars(path)
        
        if "nt" in sys.builtin_module_names:
            _expandvars = _nt_expandvars
        else:
            _expandvars = os.path.expandvars
        
        def exp_path(path):
            return _expandvars(os.path.expanduser(path))
        
        #--------------------------------------------------------------------
        # locks and threading

        # the following lock is used by both C and Python to insure that no more than
        # one active thread enters PyMOL at a given time. 

        lock_api = pymol.lock_api
        lock_api_c = pymol.lock_api_c
        lock_api_status = pymol.lock_api_status
        lock_api_glut = pymol.lock_api_glut
        
        # WARNING: internal routines, subject to change      
        def lock_c(): 
            lock_api_c.acquire(1)

        def unlock_c():
            lock_api_c.release()

        def lock_status_attempt():
            return lock_api_status.acquire(0)

        def lock_status(): 
            lock_api_status.acquire(1)

        def unlock_status():
            lock_api_status.release()

        def lock_glut(): 
            lock_api_glut.acquire(1)

        def unlock_glut():
            lock_api_glut.release()

        def lock_without_glut():
            try:
                lock_glut()
                lock()
            finally:
                unlock_glut()

        def lock(): # INTERNAL -- API lock
    #      print " lock: acquiring as 0x%x"%thread.get_ident(),(thread.get_ident() == pymol.glutThread)
            if not lock_api.acquire(0):
                w = 0.001
                while 1:
    #            print " lock: ... as 0x%x"%thread.get_ident(),(thread.get_ident() == pymol.glutThread)
                    e = threading.Event() 
                    e.wait(w)  
                    del e
                    if lock_api.acquire(0):
                        break
                    if w<0.1:
                        w = w * 2 # wait twice as long each time until flushed
    #      print "lock: acquired by 0x%x"%thread.get_ident()

        def lock_attempt(): # INTERNAL
            return lock_api.acquire(blocking=0)

        def unlock(result=None): # INTERNAL
            if (thread.get_ident() == pymol.glutThread):
                global reaper
                if reaper:
                    try:
                        if not reaper.isAlive():
                            if pymol.invocation.options.no_gui:
                                _cmd.quit()
                            else:
                                reaper = None
                    except:
                        pass
                lock_api.release()
    #         print "lock: released by 0x%x (glut)"%thread.get_ident()
                if result==None: # don't flush if we have an incipient error (negative input)
                    _cmd.flush_now()
                elif is_ok(result):
                    _cmd.flush_now()
            else:
    #         print "lock: released by 0x%x (not glut), waiting queue"%thread.get_ident()
                lock_api.release()
                if _cmd.wait_queue(): # commands waiting to be executed?
                    e = threading.Event() # abdicate control for a 100 usec for quick tasks
                    e.wait(0.0001)
                    del e
                    # then give PyMOL increasingly longer intervals to get its work done...
                    w = 0.0005  # NOTE: affects API perf. for "do" and delayed-exec
                    while _cmd.wait_queue(): 
                        e = threading.Event() # abdicate control for a 100 usec for quick tasks
                        e.wait(w)
                        del e
                        if w > 0.1: # wait up 0.2 sec max for PyMOL to flush queue
                            if _feedback(fb_module.cmd,fb_mask.debugging):
                                fb_debug.write("Debug: avoiding possible dead-lock?\n")
    #                  print "dead locked as 0x%x"%thread.get_ident()
                            break
                        w = w * 2 # wait twice as long each time until flushed

        def is_glut_thread(): # internal
            if thread.get_ident() == pymol.glutThread:
                return 1
            else:
                return 0

        def get_progress(reset=0):
            r = -1.0
            try:
                lock_status()
                r = _cmd.get_progress(int(reset))
            finally:
                unlock_status()
            return r

        def check_redundant_open(file):
            found = 0
            for a in pymol.invocation.options.deferred:
                if a == file:
                    found = 1
                    break
            for a in sys.argv:
                if a == file:
                    found = 1
                    break;
            return found

        #--------------------------------------------------------------------
        # Feedback

        class fb_action:
            set = 0
            enable = 1
            disable = 2
            push = 3
            pop = 4

        fb_action_sc = Shortcut(fb_action.__dict__.keys())

        class fb_module:

        # This first set represents internal C systems

            all                       =0
            isomesh                   =1
            map                       =2
            matrix                    =3
            mypng                     =4
            triangle                  =5
            match                     =6
            raw                       =7
            isosurface                =8
            opengl                    =9

            color                     =10
            cgo                       =11
            feedback                  =12
            scene                     =13
            threads                   =14  
            symmetry                  =15
            ray                       =16
            setting                   =17
            object                    =18
            ortho                     =19
            movie                     =20
            python                    =21
            extrude                   =22
            rep                       =23
            shaker                    =24

            coordset                  =25
            distset                   =26
            gadgetset                 =27

            objectmolecule            =30
            objectmap                 =31
            objectmesh                =32
            objectdist                =33 
            objectcgo                 =34
            objectcallback            =35
            objectsurface             =36
            objectgadget              =37
            objectslice               =38

            repangle                  =43
            repdihederal              =44
            repwirebond               =45
            repcylbond                =46
            replabel                  =47
            repsphere                 =49
            repsurface                =50
            repmesh                   =51
            repdot                    =52
            repnonbonded              =53
            repnonbondedsphere        =54
            repdistdash               =55
            repdistlabel              =56
            repribbon                 =57
            repcartoon                =58
            sculpt                    =59
            vfont                     =60

            executive                 =70
            selector                  =71
            editor                    =72

            export                    =75
            ccmd                      =76
            api                       =77   

            main                      =80  

        # This second set, with negative indices
        # represent "python-only" subsystems

            parser                    =-1
            cmd                       =-2

        fb_module_sc = Shortcut(fb_module.__dict__.keys())

        class fb_mask:
            output =              0x01 # Python/text output
            results =             0x02
            errors =              0x04
            actions =             0x08
            warnings =            0x10
            details =             0x20
            blather =             0x40
            debugging =           0x80
            everything =          0xFF

        fb_mask_sc = Shortcut(fb_mask.__dict__.keys())

        fb_dict ={}

        for a in fb_module.__dict__.keys():
            vl = getattr(fb_module,a)
            if vl<0:
                fb_dict[vl] = 0x1F # default mask

        fb_debug = sys.stderr # can redirect python debugging output elsewhere if desred...

        def feedback(action="?",module="?",mask="?"):
            '''
DESCRIPTION

    "feedback" allows you to change the amount of information output by pymol.

USAGE

    feedback action,module,mask

    action is one of ['set','enable','disable']
    module is a space-separated list of strings or simply "all"
    mask is a space-separated list of strings or simply "everything"

NOTES:

    "feedback" alone will print a list of the available module choices

PYMOL API

    cmd.feedback(string action,string module,string mask)

EXAMPLES

    feedback enable, all , debugging
    feedback disable, selector, warnings actions
    feedback enable, main, blather

DEVELOPMENT TO DO

    Add a way of querying the current feedback settings.
    Check C source code to make source correct modules are being used.
    Check C source code to insure that all output is properly
    Update Python API and C source code to use "quiet" parameter as well.
        '''
            r = None

            # validate action

            if action=="?":
                print " feedback: possible actions: \nset, enable, disable"
                act_int = 0
            else:
                act_kee = fb_action_sc.interpret(action)
                if act_kee == None:
                    print "Error: invalid feedback action '%s'."%action
                    if _raising():
                        raise QuietException
                    else:
                        return None
                elif not is_string(act_kee):
                    print "Error: ambiguous feedback action '%s'."%action
                    print action_amb
                    if _raising():
                        raise QuietException
                    else:
                        return None
                act_int = int(getattr(fb_action,act_kee))

            if (act_int<3) and ("?" in [action,module,mask]):
                if module=="?":
                    print " feedback: Please specify module names:"
                    lst = fb_module.__dict__.keys()
                    lst.sort()
                    for a in lst:
                        if a[0]!='_':
                            print "   ",a
                if mask=="?":
                    print " feedback: Please specify masks:"
                    lst = fb_mask.__dict__.keys()
                    lst.sort()
                    for a in lst:
                        if a[0]!='_':
                            print "   ",a
            else:
                if (act_int>=3):
                    module='all'
                    mask='everything'

                # validate and combine masks

                mask_int = 0
                mask_lst = string.split(mask)
                for mask in mask_lst:
                    mask_kee = fb_mask_sc.interpret(mask)
                    if mask_kee == None:
                        print "Error: invalid feedback mask '%s'."%mask
                        if _raising(): raise QuietException
                        else: return None
                    elif not is_string(mask_kee):
                        print "Error: ambiguous feedback mask '%s'."%mask
                        if _raising(): raise QuietException
                        else: return None
                    mask_int = int(getattr(fb_mask,mask_kee))

                # validate and iterate modules

                mod_lst = string.split(module)
                for module in mod_lst:
                    mod_kee = fb_module_sc.interpret(module)
                    if mod_kee == None:
                        print "Error: invalid feedback module '%s'."%module
                        if _raising(): raise QuietException
                        else: return None
                    elif not is_string(mod_kee):
                        print "Error: ambiguous feedback module '%s'."%module
                        if _raising(): raise QuietException
                        else: return None
                    mod_int = int(getattr(fb_module,mod_kee))
                    if mod_int>=0:
                        try:
                            lock()
                            r = _cmd.set_feedback(act_int,mod_int,mask_int)
                        finally:
                            unlock()
                    if mod_int<=0:
                        if mod_int:
                            if act_int==0:
                                fb_dict[mod_int] = mask_int
                            elif act_int==1:
                                fb_dict[mod_int] = fb_dict[mod_int] | mask_int
                            elif act_int==2:
                                fb_dict[mod_int] = fb_dict[mod_int] & ( 0xFF - mask_int )
                        else:
                            for mod_int in fb_dict.keys():
                                if act_int==0:
                                    fb_dict[mod_int] = mask_int
                                elif act_int==1:
                                    fb_dict[mod_int] = fb_dict[mod_int] | mask_int
                                elif act_int==2:
                                    fb_dict[mod_int] = fb_dict[mod_int] & ( 0xFF - mask_int )
                        if _feedback(fb_module.feedback,fb_mask.debugging):
                             sys.stderr.write(" feedback: mode %d on %d mask %d\n"%(
                                 act_int,mod_int,mask_int))
            return r

        #--------------------------------------------------------------------
        # internal API routines

        def _ray_anti_spawn(thread_info):
            # WARNING: internal routine, subject to change      
            # internal routine to support multithreaded raytracing
            thread_list = []
            for a in thread_info[1:]:
                t = threading.Thread(target=_cmd.ray_anti_thread,
                                            args=(a,))
                t.setDaemon(1)
                thread_list.append(t)
            for t in thread_list:
                t.start()
            _cmd.ray_anti_thread(thread_info[0])
            for t in thread_list:
                t.join()

        def _ray_hash_spawn(thread_info):
            # WARNING: internal routine, subject to change      
            # internal routine to support multithreaded raytracing
            thread_list = []
            for a in thread_info[1:]:
                if a != None:
                    t = threading.Thread(target=_cmd.ray_hash_thread,
                                         args=(a,))
                    t.setDaemon(1)
                    t.start()
                    thread_list.append(t)
            _cmd.ray_hash_thread(thread_info[0])
            for t in thread_list:
                t.join()

        def _ray_spawn(thread_info):
            # WARNING: internal routine, subject to change      
            # internal routine to support multithreaded raytracing
            thread_list = []
            for a in thread_info[1:]:
                t = threading.Thread(target=_cmd.ray_trace_thread,
                                            args=(a,))
                t.setDaemon(1)
                thread_list.append(t)
            for t in thread_list:
                t.start()
            _cmd.ray_trace_thread(thread_info[0])
            for t in thread_list:
                t.join()

        def _coordset_update_thread(list_lock,thread_info):
            # WARNING: internal routine, subject to change
            while 1:
                list_lock.acquire()
                if not len(thread_info):
                    list_lock.release()
                    break
                else:
                    info = thread_info.pop(0)
                    list_lock.release()
                _cmd.coordset_update_thread(info)

        def _coordset_update_spawn(thread_info,n_thread):
            # WARNING: internal routine, subject to change
            if len(thread_info):
                list_lock = threading.Lock() # mutex for list
                thread_list = []
                for a in range(1,n_thread):
                    t = threading.Thread(target=_coordset_update_thread,
                                                args=(list_lock,thread_info))
                    t.setDaemon(1)
                    thread_list.append(t)
                for t in thread_list:
                    t.start()
                _coordset_update_thread(list_lock,thread_info)
                for t in thread_list:
                    t.join()

        def _object_update_thread(list_lock,thread_info):
            # WARNING: internal routine, subject to change
            while 1:
                list_lock.acquire()
                if not len(thread_info):
                    list_lock.release()
                    break
                else:
                    info = thread_info.pop(0)
                    list_lock.release()
                _cmd.object_update_thread(info)
        
        def _object_update_spawn(thread_info,n_thread):
            # WARNING: internal routine, subject to change
            if len(thread_info):
                list_lock = threading.Lock() # mutex for list
                thread_list = []
                for a in range(1,n_thread):
                    t = threading.Thread(target=_object_update_thread,
                                                args=(list_lock,thread_info))
                    t.setDaemon(1)
                    thread_list.append(t)
                for t in thread_list:
                    t.start()
                _object_update_thread(list_lock,thread_info)
                for t in thread_list:
                    t.join()

        # status reporting

        def _feedback(module,mask): # feedback query routine
            # WARNING: internal routine, subject to change      
            r = 0
            module = int(module)
            mask = int(mask)
            if module>0:
                try:
                    lock()
                    r = _cmd.feedback(module,mask)
                finally:
                    unlock(-1)
            else:
                if fb_dict.has_key(module):
                    r = fb_dict[module]&mask
            return r

        # do command (while API already locked)

        def _do(cmmd,log=0,echo=1):
            return _cmd.do(cmmd,log,echo)

        # movie rendering

        def _mpng(*arg): # INTERNAL
            # WARNING: internal routine, subject to change
            try:
                lock()   
                fname = arg[0]
                if re.search("\.png$",fname):
                    fname = re.sub("\.png$","",fname)
                fname = exp_path(fname)
                if len(arg)>3:
                    preserve = arg[3]
                else:
                    preserve = 0
                r = _cmd.mpng_(str(fname),int(arg[1]),int(arg[2]),int(preserve))
            finally:
                unlock(-1)
            return r

        # loading

        
        def _load(oname,finfo,state,ftype,finish,discrete,
                     quiet=1,multiplex=0,zoom=-1):
            # WARNING: internal routine, subject to change
            # caller must already hold API lock
            # NOTE: state index assumes 1-based state
            r = DEFAULT_ERROR
            size = 0
            if ftype not in (loadable.model,loadable.brick):
                if ftype == loadable.r3d:
                    import cgo
                    obj = cgo.from_r3d(finfo)
                    if is_ok(obj):
                        r = _cmd.load_object(str(oname),obj,int(state)-1,loadable.cgo,
                                              int(finish),int(discrete),int(quiet),
                                              int(zoom))
                    else:
                        print "Load-Error: Unable to open file '%s'."%finfo
                elif ftype == loadable.cc1: # ChemDraw 3D
                    obj = io.cc1.fromFile(finfo)
                    if obj:
                        r = _cmd.load_object(str(oname),obj,int(state)-1,loadable.model,
                                              int(finish),int(discrete),
                                              int(quiet),int(zoom))
                elif ftype == loadable.moe:
                    try:
                        # BEGIN PROPRIETARY CODE SEGMENT
                        from epymol import moe

                        if (string.find(finfo,":")>1):
                            moe_file = urllib.urlopen(finfo)
                        else:
                            moe_file = open(finfo)
                        moe_str = moe_file.read()
                        moe_file.close()
                        r = moe.read_moestr(moe_str,str(oname),int(state),
                                        int(finish),int(discrete),int(quiet),int(zoom))
                        
                        # END PROPRIETARY CODE SEGMENT
                    except ImportError:
                        print "Error: .MOE format not supported by this PyMOL build."
                        if _raising(-1): raise pymol.CmdException
                        
                elif ftype == loadable.mae:
                    try:
                        # BEGIN PROPRIETARY CODE SEGMENT
                        from epymol import schrodinger

                        if (string.find(finfo,":")>1):
                            mae_file = urllib.urlopen(finfo)
                        else:
                            mae_file = open(finfo)
                        mae_str = mae_file.read()
                        mae_file.close()
                        r = schrodinger.read_maestr(mae_str,str(oname),
                                                    int(state),
                                                    int(finish),int(discrete),
                                                    int(quiet),int(zoom))
                        
                        # END PROPRIETARY CODE SEGMENT
                    except ImportError:
                        print "Error: .MAE format not supported by this PyMOL build."
                        if _raising(-1): raise pymol.CmdException
                        
                else:
                    if ftype in _load2str.keys() and (string.find(finfo,":")>1):
                        try:
                            tmp_file = urllib.urlopen(finfo)
                        except:
                            print "Error: unable to open URL '%s'"%finfo
                            traceback.print_exc()
                            return 0
                        finfo = tmp_file.read(tmp_file) # WARNING: will block and hang -- thread instead?
                        ftype = _load2str[ftype]
                        tmp_file.close()
                        
                    r = _cmd.load(str(oname),finfo,int(state)-1,int(ftype),
                                  int(finish),int(discrete),int(quiet),
                                  int(multiplex),int(zoom))
            else:
                try:
                    x = io.pkl.fromFile(finfo)
                    if isinstance(x,types.ListType) or isinstance(x,types.TupleType):
                        for a in x:
                            r = _cmd.load_object(str(oname),a,int(state)-1,
                                                        int(ftype),0,int(discrete),int(quiet),
                                                        int(zoom))
                            if(state>0):
                                state = state + 1
                        _cmd.finish_object(str(oname))
                    else:
                        r = _cmd.load_object(str(oname),x,
                                                    int(state)-1,int(ftype),
                                                    int(finish),int(discrete),
                                                    int(quiet),int(zoom))
                except:
                    print "Load-Error: Unable to load file '%s'." % finfo
            return r

        # function keys and other specials

        def _special(k,x,y,m=0): # INTERNAL (invoked when special key is pressed)
            # WARNING: internal routine, subject to change
            k=int(k)
            m=int(m)
            my_special = special
            if(m>0) and (m<5):
                my_special = (special, shft_special, ctrl_special, ctsh_special, alt_special)[m]
            if my_special.has_key(k):
                if my_special[k][1]:
                    apply(my_special[k][1],my_special[k][2],my_special[k][3])
                else:
                    key = my_special[k][0]
                    if(m>0) and (m<5):
                        key = ('','SHFT-','CTRL-','CTSH-','ALT-')[m] + key
                    if viewing.scene_dict.has_key(key):
                        scene(key)
                    elif viewing.view_dict.has_key(key):
                        view(key)
            return None

        # control keys

        def _ctrl(k):
            # WARNING: internal routine, subject to change
            if ctrl.has_key(k):
                ck = ctrl[k]
                if ck[0]!=None:
                    apply(ck[0],ck[1],ck[2])
            return None

        # alt keys

        def _alt(k):
            # WARNING: internal routine, subject to change
            if alt.has_key(k):
                ak=alt[k]
                if ak[0]!=None:
                    apply(ak[0],ak[1],ak[2])
            return None

        # writing PNG files (thread-unsafe)

        def _png(a,width=0,height=0,dpi=-1.0,ray=0,quiet=1):
            # INTERNAL - can only be safely called by GLUT thread 
            # WARNING: internal routine, subject to change
            try:
                lock()   
                fname = a
                if not re.search("\.png$",fname):
                    fname = fname +".png"
                fname = exp_path(fname)
                r = _cmd.png(str(fname),int(width),int(height),float(dpi),int(ray),int(quiet))
            finally:
                unlock(-1)
            return r

        # quitting (thread-specific)

        def _quit():
            # WARNING: internal routine, subject to change
            try:
                lock()
                try: # flush and close log if possible to avoid threading exception
                    if pymol._log_file!=None:
                        try:
                            pymol._log_file.flush()
                        except:
                            pass
                        pymol._log_file.close()
                        del pymol._log_file
                except:
                    pass
                if reaper!=None:
                    try:
                        reaper.join()
                    except:
                        pass
                r = _cmd.quit()
            finally:
                unlock(-1)
            return r

        # screen redraws (thread-specific)

        def _refresh(swap_buffers=1):  # Only call with GLUT thread!
            # WARNING: internal routine, subject to change
            r = None
            try:
                lock()
                if thread.get_ident() == pymol.glutThread:
                    if swap_buffers:
                        r = _cmd.refresh_now()
                    else:
                        r = _cmd.refresh()
                else:
                    print "Error: Ignoring an unsafe call to cmd._refresh"
            finally:
                unlock(-1)
            return r

        # stereo (platform dependent )

        def _sgi_stereo(flag): # SGI-SPECIFIC - bad bad bad
            # WARNING: internal routine, subject to change
            if sys.platform[0:4]=='irix':
                if os.path.exists("/usr/gfx/setmon"):
                    if flag:
                        mode = os.environ.get('PYMOL_SGI_STEREO','1024x768_96s')
                        os.system("/usr/gfx/setmon -n "+mode)
                    else:
                        mode = os.environ.get('PYMOL_SGI_MONO','72hz')
                        os.system("/usr/gfx/setmon -n "+mode)

        # color alias interpretation

        def _interpret_color(_self,color):
            # WARNING: internal routine, subject to change
            _validate_color_sc(_self)
            new_color = _self.color_sc.interpret(color)
            if new_color:
                if is_string(new_color):
                    return new_color
                else:
                    _self.color_sc.auto_err(color,'color')
            else:
                return color

        def _validate_color_sc(_self):
            # WARNING: internal routine, subject to change
            if _self.color_sc == None: # update color shortcuts if needed
                lst = _self.get_color_indices()
                lst.extend([('default',-1),('auto',-2),('current',-3),('atomic',-4)])
                _self.color_sc = Shortcut(map(lambda x:x[0],lst))
                color_dict = {}
                for a in lst: color_dict[a[0]]=a[1]

        def _invalidate_color_sc(_self):
            # WARNING: internal routine, subject to change
            _self.color_sc = None

        def _get_color_sc(_self):
            # WARNING: internal routine, subject to change
            _validate_color_sc(_self)
            return _self.color_sc

        def _get_feedback(): # INTERNAL
            # WARNING: internal routine, subject to change
            l = []
            if lock_attempt():
                try:
                    r = _cmd.get_feedback()
                    while r:
                        l.append(r)
                        r = _cmd.get_feedback()
                finally:
                    unlock(-1)
            return l
        get_feedback = _get_feedback # for legacy compatibility

        def _fake_drag(): # internal
            lock()
            try:
                _cmd.fake_drag()
            finally:
                unlock(-1)
            return 1

        def interrupt(): # asynch -- no locking!
            _cmd.interrupt(1)
            return None

        # testing tools

        # for comparing floating point numbers calculated using
        # different FPUs and which may show some wobble...

        def _dump_floats(lst,format="%7.3f",cnt=9):
            # WARNING: internal routine, subject to change
            c = cnt
            for a in lst:
                print format%a,
                c = c -1
                if c<=0:
                    print
                    c=cnt
            if c!=cnt:
                print

        def _dump_ufloats(lst,format="%7.3f",cnt=9):
            # WARNING: internal routine, subject to change
            c = cnt
            for a in lst:
                print format%abs(a),
                c = c -1
                if c<=0:
                    print
                    c=cnt
            if c!=cnt:
                print

        # HUH?
        def _adjust_coord(a,i,x):
            a.coord[i]=a.coord[i]+x
            return None

        _intType = types.IntType
        
        def is_error(result): # errors are always negative numbers
            if isinstance(result,_intType):
                return (result<0)
            return 0

        def is_ok(result): # something other than a negative number
            if isinstance(result,_intType):
                return (result>=0)
            return 1

        def _raising(code=-1):
            # WARNING: internal routine, subject to change
            if isinstance(code, _intType):
                if code<0:
                    return get_setting_legacy("raise_exceptions")
            return 0

        #######################################################################
        # now import modules which depend on the above
        #######################################################################

        import editor

        #######################################################################
        # cmd module functions...
        #######################################################################

        def ready(): # INTERNAL
            # WARNING: internal routine, subject to change
            return _cmd.ready()

        def setup_global_locks(): # INTERNAL, OBSOLETE?
            # WARNING: internal routine, subject to change
            pass

        def null():
            pass

        # for extending the language

        def extend(name,function):
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
            keyword[name] = [function, 0,0,',',parsing.STRICT]
            kwhash.append(name)
            help_sc.append(name)

        # for aliasing compound commands to a single keyword

        def alias(name, command):
            '''
DESCRIPTION

    "alias" allows you to bind routinely used command-line input to a
    new PyMOL command keyword.

USAGE

    alias name, literal-command-input

ARGUMENTS

    literal-command-input may contain multiple commands separated by semicolons.
    
EXAMPLE

    alias my_scene, hide; show ribbon, polymer; show sticks, organic; show nonbonded, solvent
    my_scene

NOTES

    For security reasons, aliased commands are not saved or restored
    in sessions.  

SEE ALSO

    extend, api
            '''
            keyword[name] = [eval("lambda :do('''%s ''')"%command), 0,0,',',parsing.STRICT]
            kwhash.append(name)

        def write_html_ref(file):
            lst = globals()
            f=open(file,'w')
            head = 'H2'
            f.write("<HTML><BODY><H1>Reference</H1>")
            kees = lst.keys()
            kees.sort()
            for a in kees:
                if hasattr(lst[a],'__doc__'):
                    if (a[0:1]!='_' and
                        (a not in ['string','thread',
                                   'setup_global_locks',
                                   'real_system', 'sys','imp','glob','vl','time',
                                   'threading', 'repres','re','python_help','os',
                                   'fb_debug', 'fb_dict','ctrl','auto_arg','alt','a',
                                   'help_only', 'special','stereo_dict','toggle_dict',
                                   'palette_dict', 'types' ])):
                        doc = lst[a].__doc__
                        if is_string(doc):
                            if len(doc):
                                doc = string.strip(doc)
                                doc = string.replace(doc,"<","&lt;")
                                f.write("<HR SIZE=1><%s>"%head+a+"</%s>\n"%head)
                                f.write("<PRE>"+string.strip(doc)+"\n\n</PRE>")
            f.write("</BODY></HTML>")
            f.close()

        def dummy(*arg):
            '''
DESCRIPTION

    This is a dummy function which returns None.
            '''
            #'
            return None

        def python_help(*arg):
            r'''
DESCRIPTION

    You have asked for help on a Python keyword which is available
    from within the PyMOL command language.  Please consult the
    official Python documentation at http://www.python.org for
    detailed information on Python keywords.

    You may include Python blocks in your PyMOL command scripts, but do
    note that multi-line blocks of Python in PyMOL command files will
    require explicit continuation syntax in order to execute properly
    (see below).

    Generally, if you want to write Python block which span multiple
    lines, you will want to use ".py" file, and then use "extend" in
    order to expose your new code to the PyMOL command language.  This
    will give you better error checking and more predictable results.

EXAMPLES

    a=1
    while a<10: \
        print a \
        a=a+1

SEE ALSO

    extend, run, @
            '''
            return None

        #####################################################################
        # Here is where the PyMOL Command Language and API are built.
        #####################################################################

        # first we need to import a set of symbols into the local namespace

        from api import *

        # now we create the command langauge
        
        import keywords
        keyword = keywords.get_command_keywords()
        kw_list = keyword.keys()

        # remove legacy commands from the shortcut list
        
        kw_list.remove('matrix_transfer')
        kw_list.remove('util.mroll')
        kw_list.remove('util.mrock')
        
        kwhash = Shortcut(kw_list)

        # Prepare for Python 2.6 (not hashed)

        keyword['show_as'] = keyword['as']
        
        # Aliases for Mother England (not hashed)

        keyword['colour'] = keyword['color']
        keyword['set_colour'] = keyword['set_color']
        keyword['recolour'] = keyword['recolor']
        keyword['bg_colour'] = keyword['bg_color']
    
        # informational or API-only functions which don't exist in the
        # PyMOL command language namespace

        help_only = keywords.get_help_only_keywords()
        help_sc = Shortcut(keyword.keys()+help_only.keys())

        def auto_measure():
            lst = get_names("selections")
            if "pk1" in lst:
                if "pk2" in lst:
                    if "pk3" in lst:
                        if "pk4" in lst:
                            dihedral()
                        else:
                            angle()
                    else:
                        distance()
            unpick()   

        # keyboard configuration
        
        import keyboard
        
        special = keyboard.get_special()

        shft_special = keyboard.get_shft_special()        
        alt_special = keyboard.get_alt_special()        
        ctrl_special = keyboard.get_ctrl_special()
        ctsh_special = keyboard.get_ctsh_special()

        ctrl = keyboard.get_ctrl()        
        alt = keyboard.get_alt()

        selection_sc = lambda sc=Shortcut,gn=get_names:sc(gn('public')+['all'])
        object_sc = lambda sc=Shortcut,gn=get_names:sc(gn('objects'))
        map_sc = lambda sc=Shortcut,gnot=get_names_of_type:sc(gnot('object:map'))
        contour_sc =  lambda sc=Shortcut,gnot=get_names_of_type:sc(
            gnot('object:mesh')+gnot('object:surface'))
        
        # Table for argument autocompletion

        import completing
        
        auto_arg = completing.get_auto_arg_list()

        color_sc = None

        _load2str = { loadable.pdb : loadable.pdbstr,
                      loadable.mol : loadable.molstr,
                      loadable.xplor : loadable.xplorstr,
                      loadable.mol2 : loadable.mol2str,
                      loadable.mmod : loadable.mmodstr,
                      loadable.ccp4 : loadable.ccp4str,
                      loadable.sdf2 : loadable.sdf2str}

    except:
        print "Error: unable to initalize the pymol.cmd module"
        traceback.print_exc()
        sys.exit(0)
        
else:
    from pymol.cmd import *
    

