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

        file_ext_re= re.compile(string.join([
            "\.pdb$|\.pdb1$|\.ent$|\.mol$|\.p5m$|",
            r"\.PDB$|\.pdb1$|\.ENT$|\.MOL$|\.P5M$|",
            r"\.mmod$|\.mmd$|\.dat$|\.out$|\.mol2$|",
            r"\.MMOD$|\.MMD$|\.DAT$|\.OUT$|\.MOL2$|",
            r"\.xplor$|\.pkl$|\.sdf$|\.pqr|", 
            r"\.XPLOR$|\.PKL$|\.SDF$|\.PQR|",                        
            r"\.r3d$|\.xyz$|\.xyz_[0-9]*$|", 
            r"\.R3D$|\.XYZ$|\.XYZ_[0-9]*$|",
            r"\.cc1$|\.cc2$|", # ChemDraw 3D
            r"\.CC1$|\.CC2$|",
            r"\.DX$|\.dx$|", # DX files (APBS)
            r"\.pse$|\.PSE$|", # PyMOL session (pickled dictionary)
            r"\.pmo$|", # Experimental molecular object format
            r"\.PMO$|",
            r"\.ccp4$|\.CCP4$|", # CCP4
            r"\.top$|\.TOP$|", # AMBER Topology
            r"\.trj$|\.TRJ$|", # AMBER Trajectory
            r"\.crd$|\.CRD$|", # AMBER coordinate file
            r"\.rst$|\.RST$|", # AMBER restart
            r"\.cex$|\.CEX$|", # CEX format (used by metaphorics)
            r"\.phi$|\.PHI$|", # PHI format (delphi)
            r"\.fld$|\.FLD$|", # FLD format (AVS)
            r"\.o$|\.O$|\.omap$|\.OMAP$|\.dsn6$|\.DSN6$|\.brix$|\.BRIX$", # BRIX/O format
            r"|\.grd$|\.GRD$", # InsightII Grid format
            ],''))

        reaper = None
        safe_oname_re = re.compile(r"\ |\+|\(|\)|\||\&|\!|\,")  # quash reserved characters
        sanitize_list_re = re.compile(r"[^0-9\.\-\[\]\,]+")
        sanitize_alpha_list_re = re.compile(r"[^a-zA-Z0-9\'\"\.\-\[\]\,]+")
        nt_hidden_path_re = re.compile(r"\$[\/\\]")
        
        def safe_list_eval(st):
            return eval(sanitize_list_re.sub('',st))

        def safe_alpha_list_eval(st):
            return eval(sanitize_alpha_list_re.sub('',st))

        QuietException = parsing.QuietException
        
        DEFAULT_ERROR = -1
        DEFAULT_SUCCESS = None
        
        #--------------------------------------------------------------------
        # shortcuts...

        toggle_dict = {'on':1,'off':0,'1':1,'0':0}
        toggle_sc = Shortcut(toggle_dict.keys())

        stereo_dict = {'on':1,'off':0,'1':1,'0':0,'swap':-1,
                            'crosseye':2,'quadbuffer':3,'walleye':4}
        stereo_sc = Shortcut(stereo_dict.keys())

        space_sc = Shortcut(['cmyk','rgb','pymol'])

        window_dict = { 'show' : 1, 'hide' : 0, 'position' : 2, 'size' : 3,
                             'box' : 4, 'maximize' : 5, 'fit' : 6 }
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
            result = lock_api.acquire(blocking=0)
            return result

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
                if result==None: # don't flush if we have an incipient error...
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
            for a in thread_info[1:2]:
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
                    unlock()
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
                r = _cmd.mpng_(str(fname),int(arg[1]),int(arg[2]))
            finally:
                unlock()
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

        def _png(a,width=0,height=0,quiet=1): # INTERNAL - can only be safely called by GLUT thread 
            # WARNING: internal routine, subject to change
            try:
                lock()   
                fname = a
                if not re.search("\.png$",fname):
                    fname = fname +".png"
                fname = exp_path(fname)
                r = _cmd.png(str(fname),int(width),int(height),int(quiet))
            finally:
                unlock()
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
                unlock()
            return r

        # screen redraws (thread-specific)

        def _refresh(swap_buffers=1):  # Only call with GLUT thread!
            # WARNING: internal routine, subject to change
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
                unlock()
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

        def _interpret_color(color):
            # WARNING: internal routine, subject to change
            _validate_color_sc()
            new_color = color_sc.interpret(color)
            if new_color:
                if is_string(new_color):
                    return new_color
                else:
                    color_sc.auto_err(color,'color')
            else:
                return color

        def _validate_color_sc():
            # WARNING: internal routine, subject to change
            global color_sc
            if color_sc == None: # update color shortcuts if needed
                lst = get_color_indices()
                lst.extend([('default',-1),('auto',-2),('current',-3),('atomic',-4)])
                color_sc = Shortcut(map(lambda x:x[0],lst))
                color_dict = {}
                for a in lst: color_dict[a[0]]=a[1]

        def _invalidate_color_sc():
            # WARNING: internal routine, subject to change
            global color_sc
            color_sc = None

        def _get_color_sc():
            # WARNING: internal routine, subject to change
            _validate_color_sc()
            return color_sc

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
                    unlock()
            return l
        get_feedback = _get_feedback # for legacy compatibility

        def _fake_drag(): # internal
            lock()
            try:
                _cmd.fake_drag()
            finally:
                unlock()
            return 1
        # testing tools

        # for comparing floating point numbers calculated using
        # different FPUs

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

        def alias(name,command):
            '''
DESCRIPTION

    "alias" allows you to bind a commonly used command to a single PyMOL keyword.

USAGE

    alias name, command-sequence

PYMOL API

    cmd.alias(string name,string command)

EXAMPLES

    alias go,load $TUT/1hpv.pdb; zoom 200/; show sticks, 200/ around 8
    go

NOTES

    For security reasons, new PyMOL commands created using "extend" are
    not saved or restored in sessions.

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
                    if a[0:1]!='_' and (a not in ['string','thread',
                                                            'setup_global_locks',
                                                            'real_system','sys','imp','glob','vl','time',
                                                            'threading','repres','re','python_help','os',
                                                            'fb_debug','fb_dict','ctrl','auto_arg','alt','a',
                                                            'help_only','special','stereo_dict','toggle_dict',
                                                            'palette_dict','types'
                                                            ]):
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

    You've asked for help on a Python keyword which is available from
    within the PyMOL command language.  Please consult the official
    Python documentation at http://www.python.org for detailed
    information on Python keywords.

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


        # first we need to import a set of symbols into this module's local
        # namespace

        #--------------------------------------------------------------------
        from importing import \
              finish_object,      \
              load,               \
              load_brick,         \
              load_callback,      \
              load_cgo,           \
              load_embedded,      \
              load_map,           \
              load_model,         \
              load_object,        \
              load_traj,          \
              loadable,           \
              read_mmodstr,       \
              read_molstr,        \
              read_pdbstr,        \
              read_xplorstr,      \
              set_session,        \
              space              

        #--------------------------------------------------------------------
        import creating
        from creating import \
              copy,               \
              create,             \
              fragment,           \
              isodot,             \
              isolevel,           \
              isomesh,            \
              isosurface,         \
              slice_new,          \
              symexp,             \
              map_new,            \
              ramp_new

        #--------------------------------------------------------------------
        from commanding import \
              cls,                \
              delete,             \
              do,                 \
              log,                \
              log_close,          \
              log_open,           \
              quit,               \
              resume,             \
              splash,             \
              reinitialize,       \
              sync

        #--------------------------------------------------------------------
        import controlling
        from controlling import \
              button,             \
              config_mouse,       \
              mouse,              \
              mask,               \
              order,              \
              set_key,            \
              unmask,             \
              edit_mode

        #--------------------------------------------------------------------
        from querying import \
              angle,              \
              count_atoms,        \
              count_frames,       \
              count_states,       \
              dist,               \
              dihedral,           \
              distance,           \
              export_dots,        \
              find_pairs,         \
              get_angle,          \
              get_area,           \
              get_chains,         \
              get_color_index,    \
              get_color_indices,  \
              get_object_color_index, \
              get_object_list,    \
              get_color_tuple,    \
              get_atom_coords,    \
              get_dihedral,       \
              get_distance,       \
              get_extent,         \
              get_model,          \
              get_movie_locked,   \
              get_names,          \
              get_names_of_type,  \
              get_object_matrix, \
              get_phipsi,         \
              get_position,       \
              get_povray,         \
              get_renderer,       \
              get_symmetry,       \
              get_title,          \
              get_type,           \
              id_atom,            \
              identify,           \
              index,              \
              overlap,            \
              phi_psi

        #--------------------------------------------------------------------
        from selecting import \
              deselect,           \
              indicate,           \
              select,             \
              select_list,        \
              pop

        #--------------------------------------------------------------------
        from exporting import \
              png,                \
              export_coords,      \
              get_pdbstr,         \
              get_session,        \
              multisave,          \
              save               

        #--------------------------------------------------------------------
        import editing
        from editing import \
              alter,              \
              alter_list,         \
              alter_state,        \
              attach,             \
              bond,               \
              cycle_valence,      \
              deprotect,          \
              dss,                \
              edit,               \
              fix_chemistry,      \
              flag,               \
              fuse,               \
              h_add,              \
              h_fill,             \
              invert,             \
              iterate,            \
              iterate_state,      \
              map_set_border,     \
              map_double,         \
              matrix_transfer,    \
              matrix_reset,       \
              protect,            \
              push_undo,          \
              redo,               \
              remove,             \
              remove_picked,      \
              rename,             \
              replace,            \
              rotate,             \
              sculpt_purge,       \
              sculpt_deactivate,  \
              sculpt_activate,    \
              sculpt_iterate,     \
              set_dihedral,       \
              set_name,           \
              set_geometry,       \
              set_object_ttt,     \
              set_symmetry,       \
              set_title,          \
              smooth,             \
              sort,               \
              split_states,       \
              torsion,            \
              transform_object,   \
              transform_selection,\
              translate,          \
              translate_atom,     \
              unbond,             \
              undo,               \
              unpick,             \
              update

        #--------------------------------------------------------------------

        from externing import \
              cd,                 \
              ls,                 \
              paste,              \
              pwd,                \
              system

        #--------------------------------------------------------------------
        from wizarding import \
              get_wizard,         \
              get_wizard_stack,   \
              refresh_wizard,     \
              replace_wizard,     \
              set_wizard,         \
              set_wizard_stack,   \
              dirty_wizard,       \
              wizard

        #--------------------------------------------------------------------
        from fitting import \
              align,             \
              fit,               \
              rms,               \
              rms_cur,           \
              intra_fit,         \
              intra_rms,         \
              intra_rms_cur,     \
              pair_fit          

        #--------------------------------------------------------------------
        from preset import \
              simple,            \
              technical,         \
              pretty,         \
              publication

        #--------------------------------------------------------------------
        import moving
        from moving import \
              madd,              \
              mset,              \
              mclear,            \
              mdo,               \
              mappend,           \
              mmatrix,           \
              mdump,             \
              accept,            \
              decline,           \
              mpng,              \
              mview,             \
              forward,           \
              backward,          \
              rewind,            \
              middle,            \
              ending,            \
              mplay,             \
              mtoggle,           \
              mstop,             \
              mpng,              \
              mray,              \
              frame,             \
              get_movie_playing, \
              get_state,         \
              get_frame         

        #--------------------------------------------------------------------
        import viewing
        from viewing import \
              as,                 \
              bg_color,           \
              bg_colour,          \
              cartoon,            \
              clip,               \
              color,              \
              colour,             \
              del_colorection,    \
              dirty,              \
              disable,            \
              draw,               \
              enable,             \
              full_screen,        \
              get_colorection,    \
              get_view,           \
              get_vis,            \
              get_scene_dict,     \
              hide,               \
              label,              \
              load_png,           \
              meter_reset,        \
              move,               \
              orient,             \
              origin,             \
              center,             \
              ray,                \
              rebuild,            \
              recolor,            \
              recolour,           \
              refresh,            \
              reset,              \
              rock,               \
              scene,              \
              set_color,          \
              set_colour,         \
              set_colorection,    \
              set_colorection_name,\
              set_vis,            \
              set_view,           \
              show,               \
              spectrum,           \
              stereo,             \
              toggle,             \
              turn,               \
              view,               \
              viewport,           \
              window,             \
              zoom
    #        rgbfunction,        \
    #        slice_lock,         \
    #        slice_unlock,       \
    #        slice_heightmap,    \

        #--------------------------------------------------------------------
        import setting
        from setting import \
              set,                 \
              get,                 \
              unset,               \
              get_setting_legacy,  \
              get_setting_tuple,   \
              get_setting_updates, \
              get_setting_text

        #--------------------------------------------------------------------
        import helping
        from helping import \
              show_help,           \
              help,                \
              commands

        #--------------------------------------------------------------------
        from experimenting import \
              check,              \
              dump,               \
              expfit,             \
              get_bond_print,     \
              fast_minimize,      \
              import_coords,      \
              load_coords,        \
              mem,                \
              minimize,           \
              spheroid,           \
              test

        #--------------------------------------------------------------------
        #from m4x import \
        #     metaphorics

        #--------------------------------------------------------------------
        # Modules which contain programs used explicity as "module.xxx"

        import util
        import movie

        # This is the main dictionary

        keyword = {

            # keyword : [ command, # min_arg, max_arg, separator, mode ]

            # NOTE: min_arg, max_arg, and separator, are hold-overs from the
            #       original PyMOL parser which will eventually be removed.
            #       all new commands should use NO_CHECK or STRICT modes
            #       which make much better use of built-in python features.
            'abort'         : [ dummy             , 0 , 0 , ''  , parsing.ABORT ],
            'accept'        : [ accept            , 0 , 0 , ''  , parsing.STRICT ],
            'alias'         : [ alias             , 0 , 0 , ''  , parsing.LITERAL1 ], # insecure
            'align'         : [ align             , 0 , 0 , ''  , parsing.STRICT ],
            'alter'         : [ alter             , 0 , 0 , ''  , parsing.LITERAL1 ], # insecure
            'alter_state'   : [ alter_state       , 0 , 0 , ''  , parsing.LITERAL2 ], # insecure
            'angle'         : [ angle             , 0 , 0 , ''  , parsing.STRICT ],          
            'as'            : [ as                , 0 , 0 , ''  , parsing.STRICT ],          
            'assert'        : [ python_help       , 0 , 0 , ''  , parsing.PYTHON ], 
            'attach'        : [ attach            , 0 , 0 , ''  , parsing.STRICT ],
            'backward'      : [ backward          , 0 , 0 , ''  , parsing.STRICT ],
            'bg_color'      : [ bg_color          , 0 , 0 , ''  , parsing.STRICT ],
            'bond'          : [ bond              , 0 , 0 , ''  , parsing.STRICT ],
            'break'         : [ python_help       , 0 , 0 , ''  , parsing.PYTHON ],   
            'button'        : [ button            , 0 , 0 , ''  , parsing.STRICT ],
            'cartoon'       : [ cartoon           , 0 , 0 , ''  , parsing.STRICT ],
            'cd'            : [ cd                , 0 , 0 , ''  , parsing.STRICT ],
            'center'        : [ center            , 0 , 0 , ''  , parsing.STRICT ],     
            'check'         : [ check             , 0 , 0 , ''  , parsing.STRICT ],
            'class'         : [ python_help       , 0 , 0 , ''  , parsing.PYTHON ], 
            'clip'          : [ clip              , 0 , 0 , ''  , parsing.STRICT ],
            'cls'           : [ cls               , 0 , 0 , ''  , parsing.STRICT ],
            'color'         : [ color             , 0 , 0 , ''  , parsing.STRICT ],
            'commands'      : [ helping.commands  , 0 , 0 , ''  , parsing.STRICT ],
            'config_mouse'  : [ config_mouse      , 0 , 0 , ''  , parsing.STRICT ],
            'continue'      : [ python_help       , 0 , 0 , ''  , parsing.PYTHON ],
            'copy'          : [ copy              , 0 , 0 , ''  , parsing.LEGACY ],
            'count_atoms'   : [ count_atoms       , 0 , 0 , ''  , parsing.STRICT ],
            'count_frames'  : [ count_frames      , 0 , 0 , ''  , parsing.STRICT ],   
            'count_states'  : [ count_states      , 0 , 0 , ''  , parsing.STRICT ],
            'cycle_valence' : [ cycle_valence     , 0 , 0 , ''  , parsing.STRICT ],
            'create'        : [ create            , 0 , 0 , ''  , parsing.LEGACY ],
            'decline'       : [ decline           , 0 , 0 , ''  , parsing.STRICT ],      
            'delete'        : [ delete            , 0 , 0 , ''  , parsing.STRICT ],
            'def'           : [ python_help       , 0 , 0 , ''  , parsing.PYTHON ],   
            'del'           : [ python_help       , 0 , 0 , ''  , parsing.PYTHON ],
            'deprotect'     : [ deprotect         , 0 , 0 , ''  , parsing.STRICT ],
            'deselect'      : [ deselect          , 0 , 0 , ''  , parsing.STRICT ],
            'dihedral'      : [ dihedral          , 0 , 0 , ''  , parsing.STRICT ],
            'dir'           : [ ls                , 0 , 0 , ''  , parsing.STRICT ],
            'disable'       : [ disable           , 0 , 0 , ''  , parsing.STRICT ],
            'distance'      : [ distance          , 0 , 0 , ''  , parsing.LEGACY ],
            'draw'          : [ draw              , 0 , 0 , ''  , parsing.STRICT ],
            'dss'           : [ dss               , 0 , 0 , ''  , parsing.STRICT ],
            'dump'          : [ dump              , 0 , 0 , ''  , parsing.STRICT ],
            'dummy'         : [ dummy             , 0 , 0 , ''  , parsing.STRICT ],   
            'edit'          : [ edit              , 0 , 0 , ''  , parsing.STRICT ],
            'edit_mode'     : [ edit_mode         , 0 , 0 , ''  , parsing.STRICT ],
            'elif'          : [ python_help       , 0 , 0 , ''  , parsing.PYTHON ],
            'else'          : [ python_help       , 0 , 0 , ''  , parsing.PYTHON ],
            'embed'         : [ dummy             , 0 , 3 , ',' , parsing.EMBED  ],
            'enable'        : [ enable            , 0 , 0 , ''  , parsing.STRICT ],
            'ending'        : [ ending            , 0 , 0 , ''  , parsing.STRICT ],
            'except'        : [ python_help       , 0 , 0 , ''  , parsing.PYTHON ],      
            'exec'          : [ python_help       , 0 , 0 , ''  , parsing.PYTHON ],   
            'export_dots'   : [ export_dots       , 0 , 0 , ''  , parsing.STRICT ],
            'extend'        : [ extend            , 0 , 0 , ''  , parsing.STRICT ],
            'fast_minimize' : [ fast_minimize     , 1,  4 , ',' , parsing.SIMPLE ], # TO REMOVE
            'feedback'      : [ feedback          , 0,  0 , ''  , parsing.STRICT ],
            'fit'           : [ fit               , 0 , 0 , ''  , parsing.STRICT ],
            'finally'       : [ python_help       , 0 , 0 , ''  , parsing.PYTHON ],
            'fix_chemistry' : [ fix_chemistry     , 0 , 0 , ''  , parsing.STRICT ],
            'flag'          : [ flag              , 0 , 0 , ''  , parsing.LEGACY ],
            'for'           : [ python_help       , 0 , 0 , ''  , parsing.PYTHON ],
            'fork'          : [ helping.spawn     , 1 , 2 , ',' , parsing.SPAWN  ],
            'forward'       : [ forward           , 0 , 0 , ''  , parsing.STRICT ],
            'fragment'      : [ fragment          , 0 , 0 , ''  , parsing.STRICT ],
            'full_screen'   : [ full_screen       , 0 , 0 , ''  , parsing.STRICT ],
            'fuse'          : [ fuse              , 0 , 0 , ''  , parsing.STRICT ],
            'frame'         : [ frame             , 0 , 0 , ''  , parsing.STRICT ],
            'get'           : [ get               , 0 , 0 , ''  , parsing.STRICT ],
            'get_angle'     : [ get_angle         , 0 , 0 , ''  , parsing.STRICT ],      
            'get_area'      : [ get_area          , 0 , 0 , ''  , parsing.STRICT ],
            'get_chains'    : [ get_chains        , 0 , 0 , ''  , parsing.STRICT ],
            'get_dihedral'  : [ get_dihedral      , 0 , 0 , ''  , parsing.STRICT ],
            'get_distance'  : [ get_distance      , 0 , 0 , ''  , parsing.STRICT ],
            'get_extent'    : [ get_extent        , 0 , 0 , ''  , parsing.STRICT ],
            'get_position'  : [ get_position      , 0 , 0 , ''  , parsing.STRICT ],
            'get_symmetry'  : [ get_symmetry      , 0 , 0 , ''  , parsing.STRICT ],
            'get_title'     : [ get_title         , 0 , 0 , ''  , parsing.STRICT ],   
            'get_type'      : [ get_type          , 0 , 0 , ''  , parsing.STRICT ],
            'get_view'      : [ get_view          , 0 , 0 , ''  , parsing.STRICT ],
            'global'        : [ python_help       , 0 , 0 , ''  , parsing.PYTHON ],   
            'h_add'         : [ h_add             , 0 , 0 , ''  , parsing.STRICT ],
            'h_fill'        : [ h_fill            , 0 , 0 , ''  , parsing.STRICT ],
            'help'          : [ help              , 0 , 0 , ''  , parsing.STRICT ],
            'hide'          : [ hide              , 0 , 0 , ''  , parsing.STRICT ],
            'id_atom'       : [ id_atom           , 0 , 0 , ''  , parsing.STRICT ],
            'identify'      : [ identify          , 0 , 0 , ''  , parsing.STRICT ],
            'if'            : [ python_help       , 0 , 0 , ''  , parsing.PYTHON ],
            'import'        : [ python_help       , 0 , 0 , ''  , parsing.PYTHON ],   
            'index'         : [ index             , 0 , 0 , ''  , parsing.STRICT ],
            'indicate'      : [ indicate          , 0 , 0 , ''  , parsing.STRICT ],   
            'intra_fit'     : [ intra_fit         , 0 , 0 , ''  , parsing.STRICT ],
            'intra_rms'     : [ intra_rms         , 0 , 0 , ''  , parsing.STRICT ],
            'intra_rms_cur' : [ intra_rms_cur     , 0 , 0 , ''  , parsing.STRICT ],
            'invert'        : [ invert            , 0 , 0 , ''  , parsing.STRICT ],
            'isodot'        : [ isodot            , 0 , 0 , ''  , parsing.LEGACY ],
            'isolevel'      : [ isolevel           , 0 , 0 , '' , parsing.STRICT ],      
            'isomesh'       : [ isomesh           , 0 , 0 , ''  , parsing.LEGACY ],
            'isosurface'    : [ isosurface        , 0 , 0 , ''  , parsing.LEGACY ],   
            'iterate'       : [ iterate           , 0 , 0 , ''  , parsing.LITERAL1 ], # insecure
            'iterate_state' : [ iterate_state     , 0 , 0 , ''  , parsing.LITERAL2 ], # insecure
            'label'         : [ label             , 0 , 0 , ''  , parsing.LITERAL1 ], # insecure
            'load'          : [ load              , 0 , 0 , ''  , parsing.STRICT ],
            'space'         : [ space             , 0 , 0 , ''  , parsing.STRICT ],
            'load_embedded' : [ load_embedded     , 0 , 0 , ''  , parsing.STRICT ],
            'load_png'      : [ load_png          , 0 , 0 , ''  , parsing.STRICT ],
            'load_traj'     : [ load_traj         , 0 , 0 , ''  , parsing.STRICT ],
            'log'           : [ log               , 0 , 0 , ''  , parsing.STRICT ],
            'log_close'     : [ log_close         , 0 , 0 , ''  , parsing.STRICT ],
            'log_open'      : [ log_open          , 0 , 0 , ''  , parsing.STRICT ],
            'ls'            : [ ls                , 0 , 0 , ''  , parsing.STRICT ],
            'madd'          : [ madd              , 0 , 0 , ''  , parsing.STRICT ],
            'mask'          : [ mask              , 0 , 0 , ''  , parsing.STRICT ],
            'map_set_border': [ map_set_border    , 0 , 0 , ''  , parsing.STRICT ],
            'map_double'    : [ map_double        , 0 , 0 , ''  , parsing.STRICT ],      
            'map_new'       : [ map_new           , 0 , 0 , ''  , parsing.STRICT ],    
            'mappend'       : [ mappend           , 2 , 2 , ':' , parsing.MOVIE  ],
            'matrix_reset'  : [ matrix_reset      , 0 , 0 , ''  , parsing.STRICT ],         
            'matrix_transfer': [ matrix_transfer  , 0 , 0 , ''  , parsing.STRICT ],
            'mem'           : [ mem               , 0 , 0 , ''  , parsing.STRICT ],
            'meter_reset'   : [ meter_reset       , 0 , 0 , ''  , parsing.STRICT ],
            'move'          : [ move              , 0 , 0 , ''  , parsing.STRICT ],
            'mset'          : [ mset              , 0 , 0 , ''  , parsing.STRICT ],
            'mdo'           : [ mdo               , 2 , 2 , ':' , parsing.MOVIE  ],
            'mdump'         : [ mdump             , 0 , 0 , ''  , parsing.STRICT ],      
            'mpng'          : [ mpng              , 0 , 0 , ''  , parsing.STRICT ],
            'mplay'         : [ mplay             , 0 , 0 , ''  , parsing.STRICT ],
            'mtoggle'       : [ mtoggle           , 0 , 0 , ''  , parsing.STRICT ],         
            'mray'          : [ mray              , 0 , 0 , ''  , parsing.STRICT ],
            'mstop'         : [ mstop             , 0 , 0 , ''  , parsing.STRICT ],
            'mclear'        : [ mclear            , 0 , 0 , ''  , parsing.STRICT ],
            'middle'        : [ middle            , 0 , 0 , ''  , parsing.STRICT ],
            'minimize'      : [ minimize          , 0 , 4 , ',' , parsing.SIMPLE ], # TO REMOVE
            'mmatrix'       : [ mmatrix           , 0 , 0 , ''  , parsing.STRICT ],
            'mouse'         : [ mouse             , 0 , 0 , ''  , parsing.STRICT ],
            'multisave'     : [ multisave         , 0 , 0 , ''  , parsing.STRICT ],
            'mview'         : [ mview             , 0 , 0 , ''  , parsing.STRICT ],
            'origin'        : [ origin            , 0 , 0 , ''  , parsing.STRICT ],
            'orient'        : [ orient            , 0 , 0 , ''  , parsing.STRICT ],
            'overlap'       : [ overlap           , 0 , 0 , ''  , parsing.STRICT ],
            'pair_fit'      : [ pair_fit          , 2 ,98 , ',' , parsing.SIMPLE ],
            'pass'          : [ python_help       , 0 , 0 , ''  , parsing.PYTHON ],
            'phi_psi'       : [ phi_psi           , 0 , 0 , ''  , parsing.STRICT ],
            'pop'           : [ pop               , 0 , 0 , ''  , parsing.STRICT ],
            'protect'       : [ protect           , 0 , 0 , ''  , parsing.STRICT ],
            'push_undo'     : [ push_undo         , 0 , 0 , ''  , parsing.STRICT ],   
            'pwd'           : [ pwd               , 0 , 0 , ''  , parsing.STRICT ],
            'raise'         : [ python_help       , 0 , 0 , ''  , parsing.PYTHON ],
            'ramp_new'      : [ ramp_new          , 0 , 0 , ''  , parsing.STRICT ],      
            'ray'           : [ ray               , 0 , 0 , ''  , parsing.STRICT ],
            'rebuild'       : [ rebuild           , 0 , 0 , ''  , parsing.STRICT ],
            'recolor'       : [ recolor           , 0 , 0 , ''  , parsing.STRICT ],   
            'redo'          : [ redo              , 0 , 0 , ''  , parsing.STRICT ],
            'reinitialize'  : [ reinitialize      , 0 , 0 , ''  , parsing.STRICT ],      
            'refresh'       : [ refresh           , 0 , 0 , ''  , parsing.STRICT ],
            'refresh_wizard': [ refresh_wizard    , 0 , 0 , ''  , parsing.STRICT ],
            'remove'        : [ remove            , 0 , 0 , ''  , parsing.STRICT ],
            'remove_picked' : [ remove_picked     , 0 , 0 , ''  , parsing.STRICT ],
            'rename'        : [ rename            , 0 , 0 , ''  , parsing.STRICT ],
            'order'         : [ order             , 0 , 0 , ''  , parsing.STRICT ],
            'replace'       : [ replace           , 0 , 0 , ''  , parsing.STRICT ],
            'replace_wizard': [ replace_wizard    , 0 , 0 , ''  , parsing.STRICT ],
            'reset'         : [ reset             , 0 , 0 , ''  , parsing.STRICT ],
            'resume'        : [ resume            , 0 , 0 , ''  , parsing.STRICT ],
            'return'        : [ python_help       , 0 , 0 , ''  , parsing.PYTHON ],   
            'rewind'        : [ rewind            , 0 , 0 , ''  , parsing.STRICT ],
    #      'rgbfunction'   : [ rgbfunction       , 0 , 0 , ''  , parsing.LEGACY ],         
            'rock'          : [ rock              , 0 , 0 , ''  , parsing.STRICT ],
            'rotate'        : [ rotate            , 0 , 0 , ''  , parsing.STRICT ],
            'run'           : [ helping.run       , 1 , 2 , ',' , parsing.RUN    ], # insecure
            'rms'           : [ rms               , 0 , 0 , ''  , parsing.STRICT ],
            'rms_cur'       : [ rms_cur           , 0 , 0 , ''  , parsing.STRICT ],
            'save'          : [ save              , 0 , 0 , ''  , parsing.STRICT ],
            'scene'         : [ scene             , 0 , 0 , ''  , parsing.STRICT ],
            'sculpt_purge'  : [ sculpt_purge      , 0 , 0 , ''  , parsing.STRICT ],   
            'sculpt_deactivate': [ sculpt_deactivate,0, 0 , ''  , parsing.STRICT ],
            'sculpt_activate': [ sculpt_activate  , 0 , 0 , ''  , parsing.STRICT ],
            'sculpt_iterate': [ sculpt_iterate    , 0 , 0 , ''  , parsing.STRICT ],
            'spectrum'      : [ spectrum          , 0 , 0 , ''  , parsing.STRICT ],
            'select'        : [ select            , 0 , 0 , ''  , parsing.LEGACY ],
            'set'           : [ set               , 0 , 0 , ''  , parsing.LEGACY ],
            'set_color'     : [ set_color         , 0 , 0 , ''  , parsing.LEGACY ],
            'set_dihedral'  : [ set_dihedral      , 0 , 0 , ''  , parsing.STRICT ],
            'set_name'      : [ set_name          , 0 , 0 , ''  , parsing.STRICT ],
            'set_geometry'  : [ set_geometry      , 0 , 0 , ''  , parsing.STRICT ],
            'set_symmetry'  : [ set_symmetry      , 0 , 0 , ''  , parsing.STRICT ],         
            'set_title'     : [ set_title         , 0 , 0 , ''  , parsing.STRICT ],   
            'set_key'       : [ set_key           , 0 , 0 , ''  , parsing.STRICT ], # API only
            'set_view'      : [ set_view          , 0 , 0 , ''  , parsing.STRICT ],   
            'show'          : [ show              , 0 , 0 , ''  , parsing.STRICT ],
            'slice_new'     : [ slice_new         , 0 , 0 , ''  , parsing.STRICT ],
    #      'slice_lock'    : [ slice_lock        , 0 , 0 , ''  , parsing.LEGACY ],
    #      'slice_unlock'  : [ slice_unlock      , 0 , 0 , ''  , parsing.LEGACY ],
            'smooth'        : [ smooth            , 0 , 0 , ''  , parsing.STRICT ],
            'sort'          : [ sort              , 0 , 0 , ''  , parsing.STRICT ],
            'spawn'         : [ helping.spawn     , 1 , 2 , ',' , parsing.SPAWN  ], # insecure
            'spheroid'      : [ spheroid          , 0 , 0 , ''  , parsing.STRICT ],
            'splash'        : [ splash            , 0 , 0 , ''  , parsing.STRICT ],
            'split_states'  : [ split_states      , 0 , 0 , ''  , parsing.STRICT ],
            '_special'      : [ _special          , 0 , 0 , ''  , parsing.STRICT ],
            'stereo'        : [ stereo            , 0 , 0 , ''  , parsing.STRICT ],
            'symexp'        : [ symexp            , 0 , 0 , ''  , parsing.LEGACY ],
            'system'        : [ system            , 0 , 0 , ''  , parsing.LITERAL ],
            'test'          : [ test              , 0 , 0 , ''  , parsing.STRICT ],
            'toggle'        : [ toggle            , 0 , 0 , ''  , parsing.STRICT ],      
            'torsion'       : [ torsion           , 0 , 0 , ''  , parsing.STRICT ], # vs toggle_object
            'translate'     : [ translate         , 0 , 0 , ''  , parsing.STRICT ],
            'try'           : [ python_help       , 0 , 0 , ''  , parsing.PYTHON ],
            'turn'          : [ turn              , 0 , 0 , ''  , parsing.STRICT ],
            'quit'          : [ quit              , 0 , 0 , ''  , parsing.STRICT ],
            '_quit'         : [ _quit             , 0 , 0 , ''  , parsing.STRICT ],
            'png'           : [ png               , 0 , 0 , ''  , parsing.STRICT ],
            'unbond'        : [ unbond            , 0 , 0 , ''  , parsing.STRICT ],
            'unpick'        : [ unpick            , 0 , 0 , ''  , parsing.STRICT ],
            'undo'          : [ undo              , 0 , 0 , ''  , parsing.STRICT ],
            'unmask'        : [ unmask            , 0 , 0 , ''  , parsing.STRICT ],
            'unprotect'     : [ deprotect         , 0 , 0 , ''  , parsing.STRICT ],
            'unset'         : [ unset             , 0 , 0 , ''  , parsing.STRICT ],   
            'update'        : [ update            , 0 , 0 , ''  , parsing.STRICT ],
            'view'          : [ view              , 0 , 0 , ''  , parsing.STRICT ],   
            'viewport'      : [ viewport          , 0 , 0 , ''  , parsing.STRICT ],
            'window'        : [ window            , 0 , 0 , ''  , parsing.STRICT ],         
            'while'         : [ python_help       , 0 , 0 , ''  , parsing.PYTHON ],   
            'wizard'        : [ wizard            , 0 , 0 , ''  , parsing.STRICT ],
            'zoom'          : [ zoom              , 0 , 0 , ''  , parsing.STRICT ],
        # utility programs
            'util.cbag'     : [ util.cbag         , 0 , 0 , ''  , parsing.STRICT ],
            'util.cbac'     : [ util.cbac         , 0 , 0 , ''  , parsing.STRICT ],
            'util.cbay'     : [ util.cbay         , 0 , 0 , ''  , parsing.STRICT ],
            'util.cbas'     : [ util.cbas         , 0 , 0 , ''  , parsing.STRICT ],
            'util.cbap'     : [ util.cbap         , 0 , 0 , ''  , parsing.STRICT ],
            'util.cbak'     : [ util.cbak         , 0 , 0 , ''  , parsing.STRICT ],
            'util.cbaw'     : [ util.cbaw         , 0 , 0 , ''  , parsing.STRICT ],
            'util.cbab'     : [ util.cbab         , 0 , 0 , ''  , parsing.STRICT ],
            'util.cbao'     : [ util.cbao         , 0 , 0 , ''  , parsing.STRICT ],
            'util.cbam'     : [ util.cbam         , 0 , 0 , ''  , parsing.STRICT ],
            'util.cbc'      : [ util.cbc          , 0 , 0 , ''  , parsing.STRICT ],
            'util.chainbow' : [ util.chainbow     , 0 , 0 , ''  , parsing.STRICT ],
            'util.colors'   : [ util.colors       , 0 , 0 , ''  , parsing.STRICT ],
            'util.mrock'    : [ util.mrock        , 0 , 0 , ''  , parsing.STRICT ], # LEGACY
            'util.mroll'    : [ util.mroll        , 0 , 0 , ''  , parsing.STRICT ], # LEGACY
            'util.ss'       : [ util.ss           , 0 , 0 , ''  , parsing.STRICT ],# secondary structure
            'util.rainbow'  : [ util.rainbow      , 0 , 0 , ''  , parsing.STRICT ],
        # movie programs
            'movie.rock'    : [ movie.rock        , 0 , 0 , ''  , parsing.STRICT ],
            'movie.roll'    : [ movie.roll        , 0 , 0 , ''  , parsing.STRICT ],
            'movie.load'    : [ movie.load        , 0 , 0 , ''  , parsing.STRICT ],
            'movie.zoom'    : [ movie.zoom        , 0 , 0 , ''  , parsing.STRICT ],
            'movie.screw'   : [ movie.screw       , 0 , 0 , ''  , parsing.STRICT ],
            'movie.sweep'   : [ movie.sweep       , 0 , 0 , ''  , parsing.STRICT ],
            'movie.pause'   : [ movie.pause       , 0 , 0 , ''  , parsing.STRICT ],               
            'movie.nutate'  : [ movie.nutate      , 0 , 0 , ''  , parsing.STRICT ],
            'movie.tdroll'  : [ movie.tdroll      , 0 , 0 , ''  , parsing.STRICT ],
        # activate metaphorics extensions
        #   'metaphorics'   : [ metaphorics       , 0 , 0 , ''  , parsing.STRICT ],
            }

        kwhash = Shortcut(keyword.keys())

        # Aliases for Mother England (NOTE: not hashed)

        keyword['colour'] = keyword['color']
        keyword['set_colour'] = keyword['set_color']
        keyword['recolour'] = keyword['recolor']
        keyword['bg_colour'] = keyword['bg_color']

        # informational, or API-only functions which don't exist in the
        # PyMOL command language namespace

        help_only = {  
            'api'           : [ helping.api          , 0 , 0 , '' , 0 ],
            'editing'       : [ helping.editing      , 0 , 0 , '' , 0 ],  
            'edit_keys'     : [ helping.edit_keys    , 0 , 0 , '' , 0 ],
            'examples'      : [ helping.examples     , 0 , 0 , '' , 0 ],
            'faster'        : [ helping.faster       , 0 , 0 , '' , 0 ],
            'get_pdbstr'    : [ get_pdbstr        , 0 , 0 , ''  , 0 ],
            'get_names'     : [ get_names             , 0 , 0 , '' , 0 ],
            'get_type'      : [ get_type              , 0 , 0 , '' , 0 ],
            'keyboard'      : [ helping.keyboard     , 0 , 0 , '' , 0 ],
            'launching'     : [ helping.launching    , 0 , 0 , '' , 0 ],
            'load_model'    : [ load_model            , 0 , 0 , '' , 0 ],
            'mouse'         : [ helping.mouse        , 0 , 0 , '' , 0 ],
            'movies'        : [ helping.movies       , 0 , 0 , '' , 0 ],
            'python_help'   : [ python_help          , 0 , 0 , '' , 0 ],        
            'povray'        : [ helping.povray       , 0 , 0 , '' , 0 ],
            'read_molstr'   : [ read_molstr           , 0 , 0 , '' , 0 ],
            'read_pdbstr'   : [ read_pdbstr           , 0 , 0 , '' , 0 ],      
            'release'       : [ helping.release      , 0 , 0 , '' , 0 ],   
            'selections'    : [ helping.selections   , 0 , 0 , '' , 0 ],
            'sync'          : [ sync                 , 0 , 0 , '' , 0],
            'transparency'  : [ helping.transparency , 0 , 0 , '' , 0 ],
            '@'             : [ helping.at_sign      , 0 , 0 , '' , 0 ],  
        }

        help_sc = Shortcut(keyword.keys()+help_only.keys())


        special = {
            1        : [ 'F1'        , None                   , () , {} ],
            2        : [ 'F2'        , None                   , () , {} ],
            3        : [ 'F3'        , None                   , () , {} ],
            4        : [ 'F4'        , None                   , () , {} ],
            5        : [ 'F5'        , None                   , () , {} ],
            6        : [ 'F6'        , None                   , () , {} ],
            7        : [ 'F7'        , None                   , () , {} ],
            8        : [ 'F8'        , None                   , () , {} ],
            9        : [ 'F9'        , None                   , () , {} ],
            10       : [ 'F10'       , None                   , () , {} ],
            11       : [ 'F11'       , None                   , () , {} ],
            12       : [ 'F12'       , None                   , () , {} ],
            100      : [ 'left'      , backward               , () , {} ],
            101      : [ 'up'        , None                   , () , {} ],
            102      : [ 'right'     , forward                , () , {} ],
            103      : [ 'down'      , None                   , () , {} ],
            104      : [ 'pgup'      , scene                  , ('','previous') , {} ],
            105      : [ 'pgdn'      , scene                  , ('','next') , {} ],
            106      : [ 'home'      , zoom                   , () ,  {'animate':-1} ],
            107      : [ 'end'       , mtoggle                , () , {} ],
            108      : [ 'insert'    , rock                   , () , {} ]   
        }

        shft_special = {
            1        : [ 'F1'        , None                   , () , {} ],
            2        : [ 'F2'        , None                   , () , {} ],
            3        : [ 'F3'        , None                   , () , {} ],
            4        : [ 'F4'        , None                   , () , {} ],
            5        : [ 'F5'        , None                   , () , {} ],
            6        : [ 'F6'        , None                   , () , {} ],
            7        : [ 'F7'        , None                   , () , {} ],
            8        : [ 'F8'        , None                   , () , {} ],
            9        : [ 'F9'        , None                   , () , {} ],
            10       : [ 'F10'       , None                   , () , {} ],
            11       : [ 'F11'       , None                   , () , {} ],
            12       : [ 'F12'       , None                   , () , {} ],
            100      : [ 'left'      , backward               , () , {} ],
            101      : [ 'up'        , None                   , () , {} ],
            102      : [ 'right'     , forward                , () , {} ],
            103      : [ 'down'      , None                   , () , {} ],
            104      : [ 'pgup'      , scene                  , ('','previous') , {} ],
            105      : [ 'pgdn'      , scene                  , ('','next') , {} ],
            106      : [ 'home'      , rewind                 , () ,  {'animate':-1} ],
            107      : [ 'end'       , ending                 , () , {} ],
            108      : [ 'insert'    , rock                   , () , {} ]   
        }

        alt_special = { # NOTE: some OSes/Windowing systems intercept ALT-Fn keys.
            1        : [ 'F1'        , None                   , () , {} ],
            2        : [ 'F2'        , None                   , () , {} ],
            3        : [ 'F3'        , None                   , () , {} ],
            4        : [ 'F4'        , None                   , () , {} ],
            5        : [ 'F5'        , None                   , () , {} ],
            6        : [ 'F6'        , None                   , () , {} ],
            7        : [ 'F7'        , None                   , () , {} ],
            8        : [ 'F8'        , None                   , () , {} ],
            9        : [ 'F9'        , None                   , () , {} ],
            10       : [ 'F10'       , None                   , () , {} ],
            11       : [ 'F11'       , None                   , () , {} ],
            12       : [ 'F12'       , None                   , () , {} ],
            100      : [ 'left'      , backward               , () , {} ],
            101      : [ 'up'        , None                   , () , {} ],
            102      : [ 'right'     , forward                , () , {} ],
            103      : [ 'down'      , None                   , () , {} ],
            104      : [ 'pgup'      , rewind                 , () , {} ],
            105      : [ 'pgdn'      , ending                 , () , {} ],
            106      : [ 'home'      , zoom                   , () ,  {'animate':-1} ],
            107      : [ 'end'       , ending                 , () , {} ],
            108      : [ 'insert'    , rock                   , () , {} ]   
        }

        ctrl_special = { # NOTE: some OSes/Windowing systems intercept CTRL-Fn keys.
            1        : [ 'F1'        , scene  , ('F1','store') , {} ],
            2        : [ 'F2'        , scene,('F2','store')    , {} ],
            3        : [ 'F3'        , scene,('F3','store')    , {} ],
            4        : [ 'F4'        , scene,('F4','store')    , {} ],
            5        : [ 'F5'        , scene,('F5','store')    , {} ],
            6        : [ 'F6'        , scene,('F6','store')    , {} ],
            7        : [ 'F7'        , scene,('F7','store')    , {} ],
            8        : [ 'F8'        , scene,('F8','store')    , {} ],
            9        : [ 'F9'        , scene,('F9','store')    , {} ],
            10       : [ 'F10'       , scene,('F10','store')   , {} ],
            11       : [ 'F11'       , scene,('F11','store')   , {} ],
            12       : [ 'F12'       , scene,('F12','store')   , {} ],
            100      : [ 'left'      , backward               , () , {} ],
            101      : [ 'up'        , None                   , () , {} ],
            102      : [ 'right'     , forward                , () , {} ],
            103      : [ 'down'      , None                   , () , {} ],
            104      : [ 'pgup'      , scene                  , ('','insert_before') , {} ],
            105      : [ 'pgdn'      , scene                  , ('','insert_after') , {} ],
            106      : [ 'home'      , zoom                   , () , {'animate':-1} ],
            107      : [ 'end'       , scene                  , ('new','store') , {} ],
            108      : [ 'insert'    , scene                  , ('auto','store') , {} ]   
        }

        ctsh_special = { # NOTE: some OSes/Windowing systems intercept CTRL-Fn keys.
            1        : [ 'F1'        , scene,('SHFT-F1','store') , {} ],
            2        : [ 'F2'        , scene,('SHFT-F2','store')    , {} ],
            3        : [ 'F3'        , scene,('SHFT-F3','store')    , {} ],
            4        : [ 'F4'        , scene,('SHFT-F4','store')    , {} ],
            5        : [ 'F5'        , scene,('SHFT-F5','store')    , {} ],
            6        : [ 'F6'        , scene,('SHFT-F6','store')    , {} ],
            7        : [ 'F7'        , scene,('SHFT-F7','store')    , {} ],
            8        : [ 'F8'        , scene,('SHFT-F8','store')    , {} ],
            9        : [ 'F9'        , scene,('SHFT-F9','store')    , {} ],
            10       : [ 'F10'       , scene,('SHFT-F10','store')   , {} ],
            11       : [ 'F11'       , scene,('SHFT-F11','store')   , {} ],
            12       : [ 'F12'       , scene,('SHFT-F12','store')   , {} ],
            100      : [ 'left'      , backward               , () , {} ],
            101      : [ 'up'        , None                   , () , {} ],
            102      : [ 'right'     , forward                , () , {} ],
            103      : [ 'down'      , None                   , () , {} ],
            104      : [ 'pgup'      , scene                  , ('','insert_before') , {} ],
            105      : [ 'pgdn'      , ending                 , ('','insert_after') , {} ],
            106      : [ 'home'      , zoom                 , () ,  {'animate':-1} ],
            107      : [ 'end'       , mtoggle                , () , {} ],
            108      : [ 'insert'    , rock                   , () , {} ]   
        }

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
            

        ctrl = {
            'A' : [ redo                   , () , {}],
            'B' : [ replace                , ('Br',1,1), {} ],
            'C' : [ replace                , ('C',4,4), {} ],
            'D' : [ remove_picked          , () , {'quiet':0} ],
            'E' : [ invert                 , () , {'quiet':0} ],      
            'F' : [ replace                , ('F',1,1), {} ],   
            'G' : [ replace                , ('H',1,1), {} ],
            'I' : [ replace                , ('I',1,1), {} ],
            'J' : [ alter                  , ('pk1','formal_charge=-1.0'), {} ],
            'K' : [ alter                  , ('pk1','formal_charge =1.0'), {} ],
            'L' : [ replace                , ('Cl',1,1) , {}],
            'N' : [ replace                , ('N',4,3) , {}],
            'O' : [ replace                , ('O',4,2) , {}],   
            'P' : [ replace                , ('P',4,1) , {}],
            'Q' : [ None                   , () , {}],   
            'R' : [ h_fill                 , () , {} ],   
            'S' : [ replace                , ('S',4,2) , {}],
            'T' : [ lambda a,b,c=bond,u=unpick: (c(a,b),u()) , ('pk1','pk2') , {} ],   
            'U' : [ alter                  , ('pk1','formal_charge =0.0') , {}],
            'W' : [ cycle_valence          , () , {}],
#         'X' : [ lambda a,b,c,d=dist,u=unpick:(d(a,b,c),u()), (None,'pk1','pk2') , {} ],
            'X' : [ auto_measure           , () , {} ],   
            'Y' : [ attach                 , ('H',1,1) , {} ],
            'Z' : [ undo                   , () , {} ],   
            }

        alt = {
            '1' : [ editor.attach_fragment  , ("pk1","formamide",5,0), {}],
            '2' : [ editor.attach_fragment  , ("pk1","formamide",4,0), {}],
            '3' : [ editor.attach_fragment  , ("pk1","sulfone",3,1), {}],
            '4' : [ editor.attach_fragment  , ("pk1","cyclobutane",4,0), {}],
            '5' : [ editor.attach_fragment  , ("pk1","cyclopentane",5,0), {}],
            '6' : [ editor.attach_fragment  , ("pk1","cyclohexane",7,0), {}],
            '7' : [ editor.attach_fragment  , ("pk1","cycloheptane",8,0), {}],
            '8' : [ editor.attach_fragment  , ("pk1","cyclopentadiene",5,0), {}],
            '9' : [ editor.attach_fragment  , ("pk1","benzene",6,0), {}],
            '0' : [ editor.attach_fragment  , ("pk1","formaldehyde",2,0), {}],
            'a' : [ editor.attach_amino_acid, ("pk1","ala"), {}],
            'b' : [ editor.attach_amino_acid, ("pk1","ace"), {}],                                 
            'c' : [ editor.attach_amino_acid, ("pk1","cys"), {}],
            'd' : [ editor.attach_amino_acid, ("pk1","asp"), {}],
            'e' : [ editor.attach_amino_acid, ("pk1","glu"), {}],
            'f' : [ editor.attach_amino_acid, ("pk1","phe"), {}],

            'g' : [ editor.attach_amino_acid, ("pk1","gly"), {}],
            'h' : [ editor.attach_amino_acid, ("pk1","his"), {}],
            'i' : [ editor.attach_amino_acid, ("pk1","ile"), {}],

            'j' : [ editor.attach_fragment,   ("pk1","acetylene",2,0), {}],
            'k' : [ editor.attach_amino_acid, ("pk1","lys"), {}],
            'l' : [ editor.attach_amino_acid, ("pk1","leu"), {}],

            'm' : [ editor.attach_amino_acid, ("pk1","met"), {}],
            'n' : [ editor.attach_amino_acid, ("pk1","asn"), {}],
            'p' : [ editor.attach_amino_acid, ("pk1","pro"), {}],
            'q' : [ editor.attach_amino_acid, ("pk1","gln"), {}],
            'r' : [ editor.attach_amino_acid, ("pk1","arg"), {}],

            's' : [ editor.attach_amino_acid, ("pk1","ser"), {}],
            't' : [ editor.attach_amino_acid, ("pk1","thr"), {}],
            'v' : [ editor.attach_amino_acid, ("pk1","val"), {}],
            'w' : [ editor.attach_amino_acid, ("pk1","trp"), {}],
            'y' : [ editor.attach_amino_acid, ("pk1","tyr"), {}],
            'z' : [ editor.attach_amino_acid, ("pk1","nme"), {}],
            }


        selection_sc = lambda sc=Shortcut,gn=get_names:sc(gn('public')+['all'])
        object_sc = lambda sc=Shortcut,gn=get_names:sc(gn('objects'))
        map_sc = lambda sc=Shortcut,gnot=get_names_of_type:sc(gnot('object:map'))

        # Table for argument autocompletion

        auto_arg =[
    # 1st
            {
            'align'          : [ selection_sc           , 'selection'       , ','  ],
            'alter'          : [ selection_sc           , 'selection'       , ''   ],
            'as'             : [ repres_sc              , 'representation'  , ', ' ],
            'bg'             : [ _get_color_sc          , 'color'           , ''   ],      
            'button'         : [ controlling.button_sc  , 'button'          , ', ' ],
            'color'          : [ _get_color_sc          , 'color'           , ', ' ],
            'cartoon'        : [ viewing.cartoon_sc     , 'cartoon'         , ', ' ],
            'center'         : [ selection_sc           , 'selection'       , ''   ],   
            'clip'           : [ viewing.clip_action_sc , 'clipping action' , ', ' ],
            'count_atoms'    : [ selection_sc           , 'selection'       , ''   ],
            'delete'         : [ selection_sc           , 'selection'       , ''   ],
            'deprotect'      : [ selection_sc           , 'selection'       , ''   ],   
            'full_screen'    : [ toggle_sc              , 'option'          , ''   ],
            'feedback'       : [ fb_action_sc           , 'action'          , ', ' ],
            'flag'           : [ editing.flag_sc        , 'flag'            , ', ' ],
            'get'            : [ setting.setting_sc     , 'setting'         , ','  ],      
            'help'           : [ help_sc                , 'selection'       , ''   ],
            'hide'           : [ repres_sc              , 'representation'  , ', ' ],
            'iterate'        : [ selection_sc           , 'selection'       , ''   ],
            'iterate_state'  : [ selection_sc           , 'selection'       , ''   ],
            'indicate'       : [ selection_sc           , 'selection'       , ''   ],   
            'mask'           : [ selection_sc           , 'selection'       , ''   ],
            'mview'          : [ moving.mview_action_sc , 'action'          , ''   ],
            'map_double'     : [ map_sc                 , 'map object'      , ', ' ],
            'order'          : [ selection_sc           , 'name'            , ''   ],
            'origin'         : [ selection_sc           , 'selection'       , ''   ],
            'protect'        : [ selection_sc           , 'selection'       , ''   ],
            'ramp_new'       : [ object_sc              , 'ramp'            , ''   ],
            'scene'          : [ viewing.scene_dict_sc  , 'scene'           , ''   ],
            'set'            : [ setting.setting_sc     , 'setting'         , ','  ],
            'show'           : [ repres_sc              , 'representation'  , ', ' ],
            'smooth'         : [ selection_sc           , 'selection'       , ''   ],
            'space'          : [ space_sc               , 'space'           , ''   ],      
            'split_states'   : [ object_sc              , 'object'          , ','  ],
            'stereo'         : [ stereo_sc              , 'option'          , ''   ],      
            'view'           : [ viewing.view_dict_sc   , 'view'            , ''   ],                              
            'unmask'         : [ selection_sc           , 'selection'       , ''   ],
            'unset'          : [ setting.setting_sc     , 'setting'         , ','  ],
            'window'         : [ window_sc              , 'action'          , ','  ],      
            'zoom'           : [ selection_sc           , 'selection'       , ''   ],
            },
    # 2nd
            {
            'align'          : [ selection_sc           , 'selection'       , ''   ],
            'as'             : [ selection_sc           , 'selection'       , ''   ],
            'feedback'       : [ fb_module_sc           , 'module'          , ', ' ],
            'button'         : [ controlling.but_mod_sc , 'modifier'        , ', ' ],
            'show'           : [ selection_sc           , 'selection'       , ''   ],
            'hide'           : [ selection_sc           , 'selection'       , ''   ],
            'color'          : [ selection_sc           , 'selection'       , ''   ],
            'select'         : [ selection_sc           , 'selection'       , ''   ],
            'save'           : [ selection_sc           , 'selection'       , ', ' ],
            'flag'           : [ selection_sc           , 'selection'       , ', ' ],   
            'load'           : [ selection_sc           , 'selection'       , ', ' ],
            'load_traj'      : [ selection_sc           , 'selection'       , ', ' ],
            'create'         : [ selection_sc           , 'selection'       , ', ' ],
            'map_new'        : [ creating.map_type_sc   , 'map type'        , ', ' ],
            'spectrum'       : [ palette_sc             , 'palette'         , ''   ],      
            'order'          : [ boolean_sc             , 'sort'            , ','  ],
            'symexp'         : [ object_sc              , 'object'          , ', ' ],   
            'isomesh'        : [ map_sc                 , 'map object'      , ', ' ],
            'isosurface'     : [ map_sc                 , 'map object'      , ', ' ],
            'slice_new'      : [ map_sc                 , 'map object'      , ', ' ],
            'view'           : [ viewing.view_sc        , 'view action'     , ''   ],
            'scene'          : [ viewing.scene_action_sc, 'scene action'    , ','   ],                  
            'unset'          : [ selection_sc           , 'selection'        , ','  ],
            'ramp_new'       : [ map_sc                 , 'map object'       , ''   ],      
            },
    #3rd
            {
            'spectrum'       : [ selection_sc           , 'selection'       , ''   ],
            'feedback'       : [ fb_mask_sc             , 'mask'            , ''   ],
            'order'          : [ controlling.location_sc, 'location'        , ','  ],
            'button'         : [ controlling.but_act_sc , 'button action'   , ''   ],
            'flag'           : [ editing.flag_action_sc , 'flag action'     , ''   ],
            'set'            : [ selection_sc            , 'selection'         , ','  ],
            },
    #4th
            {
            'ramp_new'       : [ creating.ramp_spectrum_sc , 'ramp color spectrum'        , ', ' ],      
            'map_new'        : [ selection_sc           , 'selection'       , ', ' ],
            'isosurface'     : [ selection_sc           , 'selection'       , ', ' ],
            'isomesh'        : [ selection_sc           , 'selection'       , ', ' ],      
            }
            ]

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
    

