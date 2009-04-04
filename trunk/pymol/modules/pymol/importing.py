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

if __name__=='pymol.importing':
    
    import re
    import string
    import os
    import cmd
    from cmd import _cmd,lock,unlock,Shortcut, \
          _feedback,fb_module,fb_mask, \
          file_ext_re,gz_ext_re,safe_oname_re, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error, \
          _load, is_list, space_sc, safe_list_eval, is_string, loadable
    import setting
    
    import selector
    try:
        from pymol import m4x
    except ImportError:
        m4x = None

    from pymol import parser
    from chempy.sdf import SDF,SDFRec
    from chempy.cif import CIF,CIFRec
    from chempy import io,PseudoFile
    import pymol
    import copy
    import traceback
    
    loadable_sc = Shortcut(loadable.__dict__.keys()) 

    def set_session(session,partial=0,quiet=1,cache=1,steal=-1,_self=cmd):
        r = DEFAULT_SUCCESS
        if is_string(session): # string implies compressed session data 
            import zlib
            session = io.pkl.fromString(zlib.decompress(session))
            if steal<0:
                steal = 1
        elif steal<0:
            steal = 0
        # use the pymol instance to store state, not the code module
        _pymol = _self._pymol
        for a in _pymol._session_restore_tasks:
            if a==None:
                try:
                    _self.lock(_self)
                    r = _cmd.set_session(_self._COb,session,int(partial),int(quiet))
                finally:
                    _self.unlock(r,_self)
                try:
                    if session.has_key('session'):
                        if steal:
                            _pymol.session = session['session']
                            del session['session']
                        else:
                            _pymol.session = copy.deepcopy(session['session'])
                    else:
                        _pymol.session = pymol.Session_Storage()
                    if cache:
                        if session.has_key('cache'):
                            cache = session['cache']
                            if len(cache):
                                if steal:
                                    _pymol._cache = session['cache']
                                    del session['cache']
                                else:
                                    _pymol._cache = copy.deepcopy(session['cache'])
                except:
                    traceback.print_exc()
            else:
                if not apply(a,(session,),{'_self':_self}): # don't stop on errors...try to complete anyway
                    r = DEFAULT_ERROR
        if _self.get_movie_locked()>0: # if the movie contains commands...activate security
            _self.wizard("security")
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def load_object(type,object,name,state=0,finish=1,discrete=0,
                         quiet=1,zoom=-1,_self=cmd):
        '''
DESCRIPTION

    "load_object" is a general developer function for loading Python objects
    into PyMOL.

PYMOL API

    cmd.load_object(type,object,name,state=0,finish=1,discrete=0,quiet=1)

    type = one one of the numeric cmd.loadable types
    object = 
    name = object name (string)
    finish = perform (1) or defer (0) post-processing of structure after load
    discrete = treat each state as an independent, unrelated set of atoms
    quiet = suppress chatter (default is yes)
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.load_object(_self._COb,str(name),object,int(state)-1,
                                        int(type),int(finish),int(discrete),
                                        int(quiet),int(zoom))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def load_brick(*arg,**kw):
        _self = kw.get('_self',cmd)
        '''
    Temporary routine for GAMESS-UK project.
    '''
        lst = [loadable.brick]
        lst.extend(list(arg))
        return apply(_self.load_object,lst,kw)

    def load_map(*arg,**kw):
        _self = kw.get('_self',cmd)
        '''
    Temporary routine for the Phenix project.
    '''

        lst = [loadable.map]
        lst.extend(list(arg))
        return apply(_self.load_object,lst,kw)

    def space(space="", gamma=1.0, quiet=0, _self=cmd):
        '''
DESCRIPTION

    "space" selects a color palette (or color space).
    
USAGE

    space space [, gamma]

ARGUMENTS

    space = rgb, cmyk, or pymol: {default: rgb}

    gamma = floating point gamma transformation
    
EXAMPLES

    space rgb
    space cmyk
    space pymol

NOTES

    Whereas computer displays use the RGB color space, computer
    printers typically use the CMYK color space.  The two spaces are
    non-equivalent, meaning that certain RGB colors cannot be
    expressed in the CMYK space and vice-versa.  And as a result,
    molecular graphics images prepared using RGB often turn out poorly
    when converted to CMYK, with purplish blues or yellowish greens.
    
    "space cmyk" forces PyMOL to restrict its use of the RGB color
    space to subset that can be reliably converted to CMYK using
    common tools such as Adobe Photoshop.  Thus, what you see on the
    screen is much closer to what you will get in print.

    Analog video systems as well as digital video compression codecs
    based on the YUV color space also have incompatibilities with RGB.
    Oversaturated colors usually cause the most problems.

    Although PyMOL lacks "space yuv", "space pymol" will help PyMOL
    avoid oversaturated colors can cause problems when exporting
    animations to video.

PYMOL API

    cmd.space(string space, float gamma)
    
SEE ALSO

    color
    
    '''
        r = DEFAULT_ERROR
        
        tables = { 'cmyk' : "$PYMOL_PATH/data/pymol/cmyk.png",
                   'pymol' : 'pymol',
                   'rgb' : 'rgb' }
        
        space_auto = space_sc.interpret(space)
        if (space_auto != None) and not is_list(space_auto):
            space = space_auto

        if space=="": 
            filename = ""
        else:         
            filename = tables.get(string.lower(space),"")
            if filename == "":
                print "Error: unknown color space '%s'."%space
                filename = None
        if filename!=None:
            try:
                if filename!="":
                    filename = _self.exp_path(filename)
                _self.lock(_self)
                r = _cmd.load_color_table(_self._COb,str(filename),float(gamma),int(quiet))
            finally:
                _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def load_callback(*arg,**kw):
        '''
DESCRIPTION

    "load_callback" is used to load a generic Python callback object.
    These objects are called every time the screen is updated and can be used
    to trigger OpenGL rendering calls (such as with PyOpenGL).

PYMOL API

    cmd.load_callback(object,name,state,finish,discrete)

    '''
        _self = kw.get('_self',cmd)
        lst = [loadable.callback]
        lst.extend(list(arg))
        return apply(_self.load_object,lst)

    def load_cgo(*arg,**kw):
        '''
DESCRIPTION

    "load_cgo" is used to load a compiled graphics object, which is
    actually a list of floating point numbers built using the constants
    in the $PYMOL_PATH/modules/pymol/cgo.py file.

PYMOL API

    cmd.load_cgo(object,name,state,finish,discrete)

    '''
        _self = kw.get('_self',cmd)
        lst = [loadable.cgo]
        lst.extend(list(arg))
        if not is_list(lst[1]): 
           lst[1] = list(lst[1]) 
        return apply(_self.load_object,lst,kw)

    def load_model(*arg,**kw):
        '''
DESCRIPTION

    "load_model" reads a ChemPy model into an object

PYMOL API

    cmd.load_model(model, object [,state [,finish [,discrete ]]])
        '''
        _self = kw.get('_self',cmd)
        lst = [loadable.model]
        lst.extend(list(arg))
        return apply(_self.load_object,lst,kw)

    def load_traj(filename,object='',state=0,format='',interval=1,
                      average=1,start=1,stop=-1,max=-1,selection='all',image=1,
                      shift="[0.0,0.0,0.0]",plugin="",_self=cmd):
        '''
DESCRIPTION

    "load_traj" reads trajectory files (currently just AMBER files).
    The file extension is used to determine the format.

    AMBER files must end in ".trj" 

USAGE

    load_traj filename [,object [,state [,format [,interval [,average ]
                             [,start [,stop [,max [,selection [,image [,shift
                             ]]]]]]]]]

PYMOL API

    cmd.load_traj(filename,object='',state=0,format='',interval=1,
                  average=1,start=1,stop=-1,max=-1,selection='all',image=1,
                  shift="[0.0,0.0,0.0]")

NOTES

    You must first load a corresponding topology file before attempting
    to load a trajectory file.

    PyMOL does not know how to wrap the truncated octahedron used by Amber
    You will need to use the "ptraj" program first to do this.

    The average option is not a running average.  To perform this type of
    average, use the "smooth" command after loading the trajectory file.

SEE ALSO

    load
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            type = format
            ftype = -1
            state = int(state)
            interval = int(interval)
            average = int(average)
            start = int(start)
            stop = int(stop)
            max = int(max)
            image = int(image)
            shift = safe_list_eval(str(shift)) # dangerous
            if is_list(shift):
                shift = [float(shift[0]),float(shift[1]),float(shift[2])]
            else:
                shift = [float(shift),float(shift),float(shift)]

            # preprocess selection
            selection = selector.process(selection)
            #   

            fname = _self.exp_path(filename)
            
            if not len(str(type)):
                # determine file type if possible
                if re.search("\.trj$",filename,re.I):
                    ftype = loadable.trj
                    try: # autodetect gromacs TRJ
                        magic = map(ord,open(fname,'r').read(4))
                        if (201 in magic) and (7 in magic):
                            ftype = loadable.trj2
                            if plugin=="": plugin = "trj"
                    except:
                        traceback.print_exc()
                elif re.search("\.xtc$",filename,re.I):
                    ftype = loadable.xtc
                    if plugin=="": plugin = "xtc" 
                elif re.search("\.trr$",filename,re.I):
                    ftype = loadable.trr
                    if plugin=="": plugin = "trr"
                elif re.search("\.gro$",filename,re.I):
                    ftype = loadable.gro
                    if plugin=="": plugin = "gro"
                elif re.search("\.g96$",filename,re.I):
                    ftype = loadable.g96
                    if plugin=="": plugin = "g96"
                elif re.search("\.dcd$",filename,re.I):
                    ftype = loadable.dcd
                    if plugin=="": plugin = "dcd"
                else:
                    raise pymol.CmdException
            elif _self.is_string(type):
                try:
                    ftype = int(type)
                except:
                    type = loadable_sc.auto_err(type,'file type')
                    if hasattr(loadable,type):
                        ftype = getattr(loadable,type)
                    else:
                        print "Error: unknown type '%s'",type
                        raise pymol.CmdException
            else:
                ftype = int(type)

    # get object name
            if len(str(object))==0:
                oname = re.sub(r".*\/|.*\\","",filename) # strip path
                oname = gz_ext_re.sub("",oname) # strip gz
                oname = file_ext_re.sub("",oname) # strip extension
                oname = safe_oname_re.sub("_",oname)
                if not len(oname): # safety
                    oname = 'obj01'
            else:
                oname = string.strip(object)

            if ftype>=0:
                r = _cmd.load_traj(_self._COb,str(oname),fname,int(state)-1,int(ftype),
                                         int(interval),int(average),int(start),
                                         int(stop),int(max),str(selection),
                                         int(image),
                                         float(shift[0]),float(shift[1]),
                                         float(shift[2]),str(plugin))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r
    
    def _processCIF(cif,oname,state,quiet,discrete,_self=cmd):
        while 1:
            rec = cif.read()
            if not rec: break
            if len(rec.model.atom):
                r = cmd.load_model(rec.model,oname,state,quiet=quiet,
                                   _self=_self,discrete=discrete)
        del cif
#        _cmd.finish_object(_self._COb,str(oname))
#        if _cmd.get_setting(_self._COb,"auto_zoom")==1.0:
#            _self._do("zoom (%s)"%oname)

    def _processSDF(sdf,oname,state,quiet,_self=cmd):
        while 1:
            rec = sdf.read()


            if not rec: break
            r = _load(oname,string.join(rec.get('MOL'),''),state,
                      loadable.molstr,0,1,quiet,_self=_self)
        del sdf
        _cmd.finish_object(_self._COb,str(oname))
        if _cmd.get_setting(_self._COb,"auto_zoom")==1.0:
            _self._do("zoom (%s)"%oname)

    def _processALN(fname,quiet=1,_self=cmd):
        legal_dict = {}
        seq_dict = {}
        seq_order = []
        header_seen = 0
        for line in open(fname).readlines():
            if not header_seen:
                if line[0:7]=='CLUSTAL':
                    header_seen = 1
            else:
                key = line[0:16].strip()
                if key!='':
                    legal_key = _self.get_legal_name(key)
                    if not legal_dict.has_key(key):
                        seq_order.append(legal_key)
                    legal_dict[key] = legal_key
                    key = legal_key
                    seq = line[16:].strip()
                    if seq != '':
                        seq_dict[key] = seq_dict.get(key,'') + seq
        for key in seq_order:
            raw_seq = seq_dict[key].replace('-','')
            _self.fab(raw_seq, key, quiet=quiet)

    def _processFASTA(fname,quiet=1,_self=cmd):
        legal_dict = {}
        seq_dict = {}
        seq_order = []
        for line in open(fname).readlines():
            line = line.strip()
            if len(line):
                if line[0:1] == '>':
                    key = line[1:].strip()
                    legal_key = _self.get_legal_name(key)
                    if not legal_dict.has_key(key):
                        seq_order.append(legal_key)
                    legal_dict[key] = legal_key
                    key = legal_key
                elif key:
                    seq = line
                    seq_dict[key] = seq_dict.get(key,'') + seq
        for key in seq_order:
            raw_seq = seq_dict[key].replace('-','')
            _self.fab(raw_seq, key, quiet=quiet)
        
    def _processPWG(fname,_self=cmd):
        r = DEFAULT_ERROR
        try:
            from consortia.web.d0903.pymolhttpd import PymolHttpd
            browser_flag = 0
            launch_flag = 0
            report_url = None
            logging = 1
            root = None
            port = 0
            if (string.find(fname,":")>1):
                import urllib
                lines = urllib.urlopen(fname).readlines()
            else:
                lines = open(fname).readlines()
            for line in lines:
                line = line.strip()
                if len(line) and line[0:1] != '#':
                    input = line.split(None,1)
                    if len(input) and input[0]!='#':
                        keyword = input[0].lower()
                        if keyword == 'port': # will be assigned dynamically if not specified
                            if len(input)>1:
                                port = int(input[1].strip())
                                launch_flag = 1
                        elif keyword == 'logging': 
                            if len(input)>1:
                                logging = int(input[1].strip())
                        elif keyword == 'root': # must encode a valid filesystem path to local content
                            if len(input)>1:
                                root = input[1].strip()
                                root = _self.exp_path(root) # allow for env var substitution
                                if os.path.exists(root):
                                    launch_flag = 1
                                else:
                                    print "Error: requested path '%s' does not exist."%root
                            else:
                                print "Error: missing path to root content"
                        elif keyword == 'browser':
                            # could perhaps interpret a browser name here
                            browser_flag = 1
                        elif keyword == 'launch': # launch the module named in the file (must exist!)
                            if len(input)>1:
                                mod_name = input[1]
                                try:
                                    __builtin__.__import__(mod_name)
                                    mod = sys.modules[mod_name]
                                    if hasattr(mod,'__launch__'):
                                        mod.__launch__(_self)
                                        r = DEFAULT_SUCCESS
                                except:
                                    traceback.print_exc()
                                    print "Error: unable to launch web application'%s'."%mode_name
                        elif keyword == 'report':
                            if len(input)>1:
                                report_url = input[1]
                        elif keyword == 'delete':
                            os.unlink(fname)
                        else:
                            print "Error: unrecognized input '%s'"%input
            if launch_flag:
                server = PymolHttpd(port,root,logging)
                if port == 0:
                    port = server.port # get the dynamically assigned port number
                server.start()
                if browser_flag: # fire up a local browser
                    import webbrowser
                    webbrowser.open("http://localhost:%d"%port,new=1)
                    r = DEFAULT_SUCCESS
                else:
                    r = DEFAULT_SUCCESS
                if report_url != None: # report port back to server url (is this secure?)
                    import urllib
                    try:
                        report_url = report_url + str(port)
                        print " Reporting back pymol port via: '%s'"%report_url
                        urllib.urlretrieve(report_url)
                    except:
                        print " Report attempt may have failed."
        except ImportError:
            traceback.print_exc()

        if is_error(r):
            print "Error: unable to handle PWG file"
        return r
    
    def load(filename, object='', state=0, format='', finish=1,
             discrete=-1, quiet=1, multiplex=None, zoom=-1, partial=0,
             _self=cmd):
        '''
DESCRIPTION

    "load" can by used to read molecules, crystallographic maps and
    other volumetric data, PyMOL sessions, and some other types of
    content.

USAGE

    load filename [,object [,state [,format ]]

ARGUMENTS

    filename = string: file path or URL

    object = string: name of the object {default: filename prefix}

    state = integer: number of the state into which
    the content should be loaded, or 0 for append {default:0}

    format = pdb, ccp4, etc. {default: use file extension}): format of
    data file    
    
EXAMPLES

    load 1dn2.pdb

    load file001.pdb, ligand

    load http://delsci.com/sample.pdb
    
NOTES

    The file extension is used to determine the format unless the
    format is provided explicitly.

    If an object name is specified, then the file is loaded into that
    object.  Otherwise, an object is created with the same name as the
    file prefix.

    If a state value is not specified, then the content is appended
    after the last existing state (if any).

    Supported molecular file formats include: pdb, mol, mol2, sdf,
    xyz, and others.

    Supported map formats include: xplor, ccp4, phi, and others.

    All supported file formats are covered in the reference
    documentation under File Formats.

PYMOL API

    cmd.load(string filename, string object-name, integer state,
             string format, int finish, int discrete, int quiet,
             int multiplex, int zoom, int partial)
    
SEE ALSO

    save
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            type = format
            ftype = 0
            state = int(state)
            finish = int(finish)
            zoom = int(zoom)
            if discrete == -1:
                discrete_default = 1
                discrete = 0
            else:
                discrete_default = 0
                discrete = int(discrete)
            if multiplex==None:
                multiplex=-2
            fname = _self.exp_path(filename)
            go_to_first_scene = 0
            if not len(str(type)):
                fname_no_gz = gz_ext_re.sub("",filename) # strip gz
                # determine file type if possible
                if re.search("\.pdb$|\.pdb1$|\.ent$|\.p5m",fname_no_gz,re.I):
                    ftype = loadable.pdb
                elif re.search("\.mol$",fname_no_gz,re.I):
                    ftype = loadable.mol
                elif re.search("\.mmod$|\.mmd$|\.dat$|\.out$",fname_no_gz,re.I):
                    ftype = loadable.mmod
                elif re.search("\.xplor$",fname_no_gz,re.I):
                    ftype = loadable.xplor
                elif re.search("\.ccp4$",fname_no_gz,re.I):
                    ftype = loadable.ccp4
                elif re.search("\.pkl$",fname_no_gz,re.I):
                    ftype = loadable.model
                elif re.search("\.r3d$",fname_no_gz,re.I):
                    ftype = loadable.r3d
                elif re.search("\.xyz$",fname_no_gz,re.I):
                    ftype = loadable.xyz
                elif re.search("\.cc1$|\.cc2$",fname_no_gz,re.I):
                    ftype = loadable.cc1
                elif re.search("\.xyz_[0-9]*$",fname_no_gz,re.I):
                    ftype = loadable.xyz
                elif re.search("\.sdf$|\.sd$",fname_no_gz,re.I): 
                    ftype = loadable.sdf2 # now using the C-based SDF reader by default...
                elif re.search("\.cex$",fname_no_gz,re.I):
                    ftype = loadable.cex
                elif re.search("\.pmo$",fname_no_gz,re.I):
                    ftype = loadable.pmo
                elif re.search("\.top$",fname_no_gz,re.I):
                    ftype = loadable.top
                elif re.search("\.trj$",fname_no_gz,re.I):
                    ftype = loadable.trj
                elif re.search("\.trr$",fname_no_gz,re.I):
                    ftype = loadable.trr
                elif re.search("\.xtc$",fname_no_gz,re.I):
                    ftype = loadable.xtc
                elif re.search("\.gro$",fname_no_gz,re.I):
                    ftype = loadable.gro
                elif re.search("\.g96$",fname_no_gz,re.I):
                    ftype = loadable.g96
                elif re.search("\.dcd$",fname_no_gz,re.I):
                    ftype = loadable.dcd
                elif re.search("\.crd$",fname_no_gz,re.I):
                    ftype = loadable.crd
                elif re.search("\.rst$",fname_no_gz,re.I):
                    ftype = loadable.crd
                elif re.search("\.pse$",fname_no_gz,re.I):
                    ftype = loadable.pse
                elif re.search("\.psw$",fname_no_gz,re.I):
                    ftype = loadable.psw
                elif re.search("\.phi$",fname_no_gz,re.I):
                    ftype = loadable.phi
                elif re.search("\.mol2$",fname_no_gz,re.I):
                    ftype = loadable.mol2
                elif re.search("\.dx$",fname_no_gz,re.I):
                    ftype = loadable.dx
                elif re.search("\.fld$",fname_no_gz,re.I):
                    ftype = loadable.fld
                elif re.search("\.pqr$",fname_no_gz,re.I):
                    ftype = loadable.pqr
                elif re.search("\.o$|\.dsn6$|\.brix$|\.omap$",fname_no_gz,re.I):
                    ftype = loadable.brix
                elif re.search("\.grd$",fname_no_gz,re.I):
                    ftype = loadable.grd
                elif re.search("\.acnt$",fname_no_gz,re.I):
                    ftype = loadable.acnt
                elif re.search("\.p1m$",fname_no_gz,re.I):
                    ftype = loadable.p1m
                elif re.search("\.png$",fname_no_gz,re.I):
                    ftype = loadable.png
                elif re.search("\.moe$",fname_no_gz,re.I):
                    ftype = loadable.moe
                elif re.search("\.mae$",fname_no_gz,re.I):
                    ftype = loadable.mae
                elif re.search("\.cube$",fname_no_gz,re.I):
                    ftype = loadable.cube
                elif re.search("\.cif$",fname_no_gz,re.I):
                    ftype = loadable.cif1
                elif re.search("\.pim$",fname_no_gz,re.I):
                    ftype = loadable.pim
                elif re.search("\.pwg$",fname_no_gz,re.I):
                    ftype = loadable.pwg
                elif re.search("\.aln$",fname_no_gz,re.I):
                    ftype = loadable.aln
                elif re.search("\.fasta$",fname_no_gz,re.I):
                    ftype = loadable.fasta
                elif re.search("\.map$",fname_no_gz,re.I):
                    r = DEFAULT_ERROR
                    print 'Error: .map is ambiguous.  Please add format or use another extension:'
                    print 'Error: For example, "load fofc.map, format=ccp4" or "load 2fofc.xplor".'
                    
                    if _self._raising(r,_self):
                        raise pymol.CmdException
                    else:
                        return r
                else:
                    ftype = loadable.pdb # default is PDB
            elif _self.is_string(type):
                if hasattr(loadable,type):
                    ftype = getattr(loadable,type)
                else:
                    try:
                        ftype = int(type) # for some reason, these exceptions aren't always caught...
                    except:
                        type = loadable_sc.auto_err(type,'file type')
                        if hasattr(loadable,type):
                            ftype = getattr(loadable,type)
                        else:
                            print "Error: unknown type '%s'",type
                            raise pymol.CmdException
            else:
                ftype = int(type)

    # special handling for PSW files 
            if ftype == loadable.psw:
                go_to_first_scene = 1                
                ftype = loadable.pse
            elif ftype == loadable.pse:
                if int(_self.get_setting_legacy("presentation"))!=0:
                    go_to_first_scene = 1
                    
    # get object name
            if len(str(object))==0:
                oname = re.sub(r".*\/|.*\\","",filename) # strip path
                oname = gz_ext_re.sub("",oname) # strip gz                
                oname = file_ext_re.sub("",oname) # strip extension
                oname = safe_oname_re.sub("_",oname)
                if not len(oname): # safety
                    oname = 'obj01'
            else:
                oname = string.strip(object)

    # loadable.sdf1 is for the old Python-based SDF file reader
            if ftype == loadable.sdf1:
                sdf = SDF(fname)
                _processSDF(sdf,oname,state,quiet,_self)
                r = DEFAULT_SUCCESS
                ftype = -1

     # loadable.sdf1 is the Python-based CIF file reader
            if ftype == loadable.cif1:
                try:
                    cif = CIF(fname)
                    _processCIF(cif,oname,state,quiet,discrete,_self)
                except:
                    traceback.print_exc()
                r = DEFAULT_SUCCESS
                ftype = -1

    # png images 
            if ftype == loadable.png:
                r = _self.load_png(str(fname),quiet=quiet)
                ftype = -1

    # p1m embedded data script files (more secure)
    
            if ftype == loadable.p1m:
                _self._do("_ @"+fname)
                ftype = -1
                r = DEFAULT_SUCCESS

    # pim import files (unrestricted scripting -- insecure)
    
            if ftype == loadable.pim:
                _self._do("_ @"+fname)
                ftype = -1
                r = DEFAULT_SUCCESS

    # pwg launch (PyMOL Web GUI / http server launch)
    
            if ftype == loadable.pwg:
                ftype = -1
                r = _processPWG(fname)

    # aln CLUSTAL
    
            if ftype == loadable.aln:
                ftype = -1
                r = _processALN(fname,quiet=quiet)

    # fasta
    
            if ftype == loadable.fasta:
                ftype = -1
                r = _processFASTA(fname,quiet=quiet)

    # special handling for trj failes (autodetect AMBER versus GROMACS)
            if ftype == loadable.trj:
                try: # autodetect gromacs TRJ
                    magic = map(ord,open(fname,'r').read(4))
                    if (201 in magic) and (7 in magic):
                        ftype = loadable.trj2
                except:
                    traceback.print_exc()

    # special handling of cex files
            if ftype == loadable.cex:
                ftype = -1
                if m4x!=None:
                    r = m4x.readcex(fname,str(oname)) # state, format, discrete?
                else:
                    print " Error: CEX format not currently supported"
                    raise pymol.CmdException

    # special handling of pse files

            if ftype == loadable.pse:
                ftype = -1
                r = _self.set_session(io.pkl.fromFile(fname),quiet=quiet,
                                      partial=partial,steal=1)
                if not partial:
                    _self.set("session_file",fname,quiet=1)
                
    # special handling for multi-model files 

            if ftype in ( loadable.mol2, loadable.sdf1, loadable.sdf2, loadable.mae ):
                if discrete_default: # make such files discrete by default
                    discrete = -1

    # standard file handling
            if ftype>=0:
                r = _load(oname,fname,state,ftype,finish,
                          discrete,quiet,multiplex,zoom,_self=_self)
        finally:
            _self.unlock(r,_self)
        if go_to_first_scene:
            if int(_self.get_setting_legacy("presentation_auto_start"))!=0:
                _self.scene("auto","start",animate=0)

        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def load_embedded(key=None, name=None, state=0, finish=1, discrete=1,
                      quiet=1, zoom=-1, multiplex=-2, _self=cmd):
        '''
DESCRIPTION

    "load_embedded" loads content previously defined in the current
    PyMOL command script using the "embed" command.

USAGE

    load_embedded [ key [, name [, state [, finish [, discrete [, quiet ]]]]]]        

EXAMPLE

    embed wats, pdb
    HETATM    1  O   WAT     1       2.573  -1.034  -1.721
    HETATM    2  H1  WAT     1       2.493  -1.949  -1.992
    HETATM    3  H2  WAT     1       2.160  -0.537  -2.427
    HETATM    4  O   WAT     2       0.705   0.744   0.160
    HETATM    5  H1  WAT     2      -0.071   0.264   0.450
    HETATM    6  H2  WAT     2       1.356   0.064  -0.014
    embed end

    load_embedded wats

NOTES

    This approach only works with text data files.
    
    '''
        r = DEFAULT_ERROR
        list = _self._parser.get_embedded(key)
        if list == None:
            print "Error: embedded data '%s' not found."%key
        else:
            if name == None:
                if key != None:
                    name = key
                else:
                    name = _self.cmd._parser.get_default_key()
            type = list[0]
            data = list[1]
            try:
                ftype = int(type)
            except:
                type = loadable_sc.auto_err(type,'file type')
                if hasattr(loadable,type):
                    ftype = getattr(loadable,type)
                else:
                    print "Error: unknown type '%s'",type
                    raise pymol.CmdException
            if ftype==loadable.pdb:
                r = read_pdbstr(string.join(data,''),name,state,finish,
                                discrete,quiet,zoom,multiplex)
            elif ftype==loadable.mol:
                r = read_molstr(string.join(data,''),name,state,finish,
                                discrete,quiet,zoom)
            elif ftype==loadable.mol2:
                r = read_mol2str(string.join(data,''),name,state,finish,
                                 discrete,quiet,zoom,multiplex)
            elif ftype==loadable.xplor:
                r = read_xplorstr(string.join(data,''),name,state,finish,discrete,quiet)
            elif ftype==loadable.mae:
                try:
                    # BEGIN PROPRIETARY CODE SEGMENT
                    from epymol import mae
                    if discrete_default:
                        discrete = -1
                    r = mae.read_maestr(string.join(data,''),
                                        name,state,finish,discrete,
                                        quiet,zoom,multiplex)
                    # END PROPRIETARY CODE SEGMENT
                except ImportError:
                    print "Error: .MAE format not supported by this PyMOL build."
                    if raising(-1): raise pymol.CmdException
            elif ftype==loadable.sdf1: # Python-based SDF reader
                sdf = SDF(PseudoFile(data),'pf')
                r = _processSDF(sdf,name,state,quiet,_self)
            elif ftype==loadable.sdf2: # C-based SDF reader (much faster)
                r = read_sdfstr(string.join(data,''),name,state,finish,
                                discrete,quiet,zoom,multiplex)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    _raw_dict = {
        loadable.pdb  : loadable.pdbstr,
        loadable.mol  : loadable.molstr,
        loadable.sdf  : loadable.sdf2str,
        loadable.ccp4 : loadable.ccp4str,
        loadable.xplor: loadable.xplorstr
        }

    def load_raw(content,  format='', object='', state=0, finish=1,
                 discrete=-1, quiet=1, multiplex=None, zoom=-1,_self=cmd):
        r = DEFAULT_ERROR
        if multiplex==None:
            multiplex=-2
        ftype = None
        type = loadable_sc.auto_err(format,'data format')
        if hasattr(loadable,type):
            ftype = getattr(loadable,type)
        else:
            print "Error: unknown format '%s'",format
            if _self._raising(r,_self): raise pymol.CmdException            
        if ftype!=None:
            if ftype in (loadable.pdb,
                         loadable.mol,
                         loadable.sdf,
                         loadable.ccp4,
                         loadable.xplor):
                try:
                    _self.lock(_self)
                    r = _cmd.load(_self._COb,str(object),str(content),int(state)-1,
                                  _raw_dict[ftype],int(finish),int(discrete),
                                  int(quiet),int(multiplex),int(zoom))
                finally:
                    _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r
        
    def read_sdfstr(sdfstr,name,state=0,finish=1,discrete=1,quiet=1,
                    zoom=-1,multiplex=-2,_self=cmd):
        '''
DESCRIPTION

    "read_sdfstr" reads an MDL MOL format file as a string

PYMOL API ONLY

    cmd.read_sdfstr( string molstr, string name, int state=0,
        int finish=1, int discrete=1 )

NOTES

    "state" is a 1-based state index for the object, or 0 to append.

    "finish" is a flag (0 or 1) which can be set to zero to improve
    performance when loading large numbers of objects, but you must
    call "finish_object" when you are done.

    "discrete" is a flag (0 or 1) which tells PyMOL that there will be
    no overlapping atoms in the file being loaded.  "discrete"
    objects save memory but can not be edited.
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.load(_self._COb,str(name),str(sdfstr),int(state)-1,
                              loadable.sdf2str,int(finish),int(discrete),
                              int(quiet),int(multiplex),int(zoom))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def read_molstr(molstr,name,state=0,finish=1,discrete=1,quiet=1,
                         zoom=-1,_self=cmd):
        '''
DESCRIPTION

    "read_molstr" reads an MDL MOL format file as a string

PYMOL API ONLY

    cmd.read_molstr( string molstr, string name, int state=0,
        int finish=1, int discrete=1 )

NOTES

    "state" is a 1-based state index for the object, or 0 to append.

    "finish" is a flag (0 or 1) which can be set to zero to improve
    performance when loading large numbers of objects, but you must
    call "finish_object" when you are done.

    "discrete" is a flag (0 or 1) which tells PyMOL that there will be
    no overlapping atoms in the file being loaded.  "discrete"
    objects save memory but can not be edited.
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.load(_self._COb,str(name),str(molstr),int(state)-1,
                          loadable.molstr,int(finish),int(discrete),
                          int(quiet),0,int(zoom))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def read_mmodstr(*arg,**kw):
        '''
DESCRIPTION

    "read_mmodstr" reads a macromodel format structure from a Python
    string.

    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            ftype = loadable.mmodstr
            if kw.has_key('quiet'):
                quiet = int(kw['quiet'])
            else:
                quiet = 1
            if kw.has_key('zoom'):
                zoom = int(kw['zoom'])
            else:
                zoom = -1
            if len(arg)==2:
                oname = string.strip(arg[1])
                r = _cmd.load(_self._COb,str(oname),arg[0],-1,int(ftype),1,1,quiet,0,zoom)
            elif len(arg)==3:
                oname = string.strip(arg[1])
                r = _cmd.load(_self._COb,str(oname),arg[0],int(arg[2])-1,
                                  int(ftype),1,1,quiet,0,zoom)
            else:
                print "argument error."
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def read_pdbstr(pdb,name,state=0,finish=1,discrete=0,quiet=1,
                         zoom=-1,multiplex=-2,_self=cmd):
        '''
DESCRIPTION

    "read_pdbstr" in an API-only function which reads a pdb file from a
    Python string.  This feature can be used to load or update
    structures into PyMOL without involving any temporary files.

PYMOL API ONLY

    cmd.read_pdbstr( string pdb-content, string object name 
        [ ,int state [ ,int finish [ ,int discrete ] ] ] )

NOTES

    "state" is a 1-based state index for the object.

    "finish" is a flag (0 or 1) which can be set to zero to improve
    performance when loading large numbers of objects, but you must
    call "finish_object" when you are done.

    "discrete" is a flag (0 or 1) which tells PyMOL that there will be
    no overlapping atoms in the PDB files being loaded.  "discrete"
    objects save memory but can not be edited.
    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            ftype = loadable.pdbstr
            oname = string.strip(str(name))
            r = _cmd.load(_self._COb,str(oname),pdb,int(state)-1,int(ftype),
                              int(finish),int(discrete),int(quiet),
                              int(multiplex),int(zoom))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def read_mol2str(mol2,name,state=0,finish=1,discrete=0,
                          quiet=1,zoom=-1,multiplex=-2,_self=cmd):
        '''
DESCRIPTION

    "read_mol2str" in an API-only function which reads a mol2 file from a
    Python string.  This feature can be used to load or update
    structures into PyMOL without involving any temporary files.

PYMOL API ONLY

    cmd.read_mol2str( string mol2-content, string object name 
        [ ,int state [ ,int finish [ ,int discrete ] ] ] )

NOTES

    "state" is a 1-based state index for the object.

    "finish" is a flag (0 or 1) which can be set to zero to improve
    performance when loading large numbers of objects, but you must
    call "finish_object" when you are done.

    "discrete" is a flag (0 or 1) which tells PyMOL that there will be
    no overlapping atoms in the MOL2 files being loaded.  "discrete"
    objects save memory but can not be edited.
    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            ftype = loadable.mol2str
            oname = string.strip(str(name))
            r = _cmd.load(_self._COb,str(oname),mol2,int(state)-1,int(ftype),
                              int(finish),int(discrete),int(quiet),
                              int(multiplex),int(zoom))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def read_xplorstr(xplor,name,state=0,finish=1,discrete=0,
                            quiet=1,zoom=-1,_self=cmd):
        '''
DESCRIPTION

    "read_xplorstr" in an API-only function which reads an XPLOR map
    from a Python string.  This feature can be used to bypass
    temporary files.

PYMOL API ONLY

    cmd.read_xplorstr( string xplor-content, string object name 
        [ ,int state ] )

NOTES

    "state" is a 1-based state index for the object.

    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            ftype = loadable.xplorstr
            oname = string.strip(str(name))
            r = _cmd.load(_self._COb,str(oname),xplor,int(state)-1,int(ftype),
                              int(finish),int(discrete),int(quiet),
                              0,int(zoom))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def finish_object(name,_self=cmd):
        '''
DESCRIPTION

    "finish_object" is used in cases where many individual states are
    being loaded and it is advantageos to avoid processing them until
    all states have been loaded into RAM.  This function should always
    be called after loading an object with the finish flag set to zero.

PYMOL API

    cmd.finish(string name)

    "name" should be the name of the object
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.finish_object(_self._COb,name)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def _fetch(code,name,state,finish,discrete,multiplex,zoom,type,path,file,quiet,_self=cmd):
        import urllib
        import gzip
        import os
        import string
        import time
        
        fobj = None
        fname = None
        auto_close_file = 1
        if path and not file:
            file = 1
        if (file==1) or (file=='1') or (file=='auto'):
            if path:
                fname = os.path.join(path,string.lower(code)+"."+type)
            else:
                fname = string.lower(code)+"."+type
        elif is_string(file):
            fname = file
        elif file:
            fobj = file
            auto_close_file = 0
        if fname and not fobj:
            if os.path.exists(fname):
                return _self.load(fname,name,state,'pdb',finish,discrete,quiet,
                                multiplex,zoom)
        tries = 0
        r = DEFAULT_ERROR
        done = 0
        while (done == 0) and (tries<3): # try loading URL up to 3 times
            tries = tries + 1
            if (type=='pdb') or (type=='pdb1'):
                remoteCode = string.upper(code)
                try:
                    url = None
                    if type=='pdb':
#                        filename = urllib.urlretrieve(
#                            'http://www.rcsb.org/pdb/cgi/export.cgi/' +
#                            remoteCode + '.pdb.gz?format=PDB&pdbId=' +
#                            remoteCode + '&compression=gz')[0]
#                        url = ('http://www.rcsb.org/pdb/download/downloadFile.do?' +
#                               'fileFormat=PDB&structureId=' +
#                               remoteCode + '&compression=GZIP')
                        url = ('http://www.rcsb.org/pdb/files/'+
                               remoteCode + '.pdb.gz')
                    elif type=='pdb1':
                        url = ('http://www.rcsb.org/pdb/files/' +
                               remoteCode + '.pdb1.gz')
                    if url!=None:
                        filename = urllib.urlretrieve(url)[0]
                except:
#                    print traceback.print_exc()
                    pass
                else:
                    if os.path.exists(filename):
                        if (os.path.getsize(filename) > 0): # If 0, then pdb code was invalid
                            try:
                                abort = 0
                                pdb_str = gzip.open(filename).read()
                                if fname and not fobj:
                                    fobj = open(fname,'wb')
                                if fobj:
                                    fobj.write(pdb_str)
                                    fobj.flush()
                                    if auto_close_file:
                                        fobj.close()
                                r = _self.read_pdbstr(pdb_str,name,state,finish,discrete,quiet,
                                                      zoom,multiplex)
                                done = 1
                            except IOError:
#                                print traceback.print_exc()
                                pass
                        else:
#                            print traceback.print_exc()
                            pass
                        try:
                            os.remove(filename)
                        except:
                            pass
            if done == 0:
                time.sleep(0.1)
        if done == 0:
            r = DEFAULT_ERROR
            print "Error-fetch: unable to load '%s'."%code
        return r
    
    def _multifetch(code,name,state,finish,discrete,multiplex,zoom,type,path,file,quiet,_self):
        import string
        r = DEFAULT_SUCCESS
        code_list = string.split(str(code))
        name = string.strip(str(name))
        if (name!='') and (len(code_list)>1) and (discrete<0):
            discrete = 1 # by default, select discrete  when loading
            # multiple PDB entries into a single object
        for obj_code in code_list:
            obj_code = string.strip(obj_code)
            if len(obj_code):
                if name=='':
                    obj_name = obj_code
                else:
                    obj_name = name
                r = _fetch(obj_code,obj_name,state,finish,
                           discrete,multiplex,zoom,type,path,file,quiet,_self)
        return r
    
    def fetch(code, name='', state=0, finish=1, discrete=-1,
              multiplex=-2, zoom=-1, type='pdb', async=-1, path=None,
              file=None, quiet=1, _self=cmd):
        
        '''
DESCRIPTION

    "fetch" downloads a file from the internet (if possible)

USAGE

    fetch code [,name [,state]]

ARGUMENTS

    code = a single PDB identifier or a list of identifiers.

    name = the object name into which the file should be loaded.

    state = the state number into which the file should loaded.

PYMOL API

    cmd.fetch(string code, string name, int state, init finish,
              int discrete, int multiplex, int zoom, string type,
              int async, string path, string file, int quiet)
              
NOTES

    When running in interactive mode, the fetch command loads
    structures asyncronously by default, meaning that the next command
    may get executed before the structures have been loaded.  If you
    need synchronous behavior in order to insure that all structures
    are loaded before the next command is executed, please provide the
    optional argument "async=0".

    Fetch requires a direct connection to the internet and thus may
    not work behind certain types of network firewalls.
    
        '''
        import threading
        r = DEFAULT_SUCCESS
        if path==None:
            path = setting.get('fetch_path',_self=_self)
        if async<0: # by default, run asynch when interactive, sync when not
            async = not quiet
        if not int(async):
            try:
                _self.block_flush(_self)
                r = _multifetch(code,name,state,finish,
                                discrete,multiplex,zoom,type,path,file,quiet,_self)
            except:
                _self.unblock_flush(_self)
        else:
            t = threading.Thread(target=_multifetch,
                                 args=(code,name,state,finish,
                                       discrete,multiplex,zoom,type,path,file,quiet,_self))
            t.setDaemon(1)
            t.start()
        return r
        
