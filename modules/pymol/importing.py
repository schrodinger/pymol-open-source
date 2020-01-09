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

from __future__ import print_function, absolute_import

if True:

    import re
    import os
    import sys
    import copy
    import traceback
    import pymol
    cmd = sys.modules["pymol.cmd"]
    from . import selector
    from . import colorprinting
    from .cmd import _cmd, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error, \
          is_list, space_sc, safe_list_eval, is_string, loadable
    from .constants import _loadable
    from pymol.creating import unquote

    def incentive_format_not_available_func(format=''):
        raise pymol.IncentiveOnlyException(
                "'%s' format not supported by this PyMOL build" % format)

    from chempy import io

    def filename_to_objectname(fname, _self=cmd):
        oname, _, _, _ = filename_to_format(fname)
        return _self.get_legal_name(oname)

    def filename_to_format(filename):
        filename = os.path.basename(filename)
        pre, delim, ext = filename.rpartition('.')

        if ext in ('gz', 'bz2',):
            zipped = ext
            pre, delim, ext = pre.rpartition('.')
        else:
            zipped = ''

        if not pre:
            pre = ext

        ext = ext.lower()

        if ext in ('brick', 'callback', 'cgo', 'model', 'plugin'):
            # names of special loadables, not accepted as file extensions
            format = ''
        elif ext in ('ent', 'p5m'):
            format = 'pdb'
        elif ext in ('pze',):
            zipped = 'gz'
            format = 'pse'
        elif ext in ('pzw',):
            zipped = 'gz'
            format = 'psw'
        elif ext in ('mmd', 'out', 'dat',):
            format = 'mmod'
        elif ext in ('cc2',):
            format = 'cc1'
        elif ext in ('sd',):
            format = 'sdf'
        elif ext in ('rst7',):
            format = 'rst'
        elif ext in ('o', 'dsn6', 'omap',):
            format = 'brix'
        elif ext in ('maegz',):
            zipped = 'gz'
            format = 'mae'
        elif ext in ('ph4',):
            format = 'moe'
        elif ext in ('spi',):
            format = 'spider'
        elif ext in ('pym', 'pyc',):
            format = 'py'
        elif ext in ('p1m', 'pim',):
            format = 'pml'
        elif ext in ('xml',):
            format = 'pdbml'
        elif re.match(r'pdb\d+$', ext):
            format = 'pdb'
        elif re.match(r'xyz_\d+$', ext):
            format = 'xyz'
        else:
            format = ext

        return pre, ext, format, zipped

    def check_gromacs_trj_magic(filename):
        colorprinting.warning('check_gromacs_trj_magic is deprecated')
        return check_trj_magic(filename)[0] == 'trj'

    def check_trj_magic(filename):
        try:
            magic = open(filename, 'rb').read(4)
            if b'\xc9' in magic and b'\x07' in magic:
                return 'trj', loadable.trj2
            if magic in (b'CDF\001', b'CDF\002'):
                return 'netcdf', loadable.plugin
        except IOError as e:
            print('trj magic test failed: ' + str(e))
        return '', loadable.trj

    def _guess_trajectory_object(candidate, _self):
        # if candidate does not exist as a molecular object, return last
        # structure (presumably most recently added)
        onames = _self.get_object_list()
        if onames and candidate not in onames:
            return onames[-1]
        return candidate

    def auto_zoom(zoom, selection, state=0, _self=cmd):
        if zoom > 0 or zoom < 0 and _self.get_setting_int("auto_zoom"):
            _self.zoom(selection, state=state)

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
            if a is None:
                try:
                    _self.lock(_self)
                    r = _cmd.set_session(_self._COb,session,int(partial),int(quiet))
                finally:
                    _self.unlock(r,_self)
                try:
                    if 'session' in session:
                        if steal:
                            _pymol.session = session['session']
                            del session['session']
                        else:
                            _pymol.session = copy.deepcopy(session['session'])
                    else:
                        _pymol.session = pymol.Session_Storage()
                    if cache:
                        if 'cache' in session:
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
                if not a(session, _self=_self): # don't stop on errors...try to complete anyway
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
        return _self.load_object(*lst, **kw)

    def load_map(*arg,**kw):
        _self = kw.get('_self',cmd)
        '''
    Temporary routine for the Phenix project.
    '''

        lst = [loadable.chempymap]
        lst.extend(list(arg))
        return _self.load_object(*lst, **kw)

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

        tables = { 'cmyk' : "$PYMOL_DATA/pymol/cmyk.png",
                   'pymol' : 'pymol',
                   'rgb' : 'rgb',
                   'greyscale': 'greyscale' }

        space_auto = space_sc.interpret(space)
        if (space_auto is not None) and not is_list(space_auto):
            space = space_auto

        if space=="":
            filename = ""
        else:
            filename = tables.get(space.lower(),"")
            if filename == "":
                print("Error: unknown color space '%s'."%space)
                filename = None
        if filename is not None:
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
        return _self.load_object(*lst)

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
        return _self.load_object(*lst, **kw)

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
        return _self.load_object(*lst, **kw)

    def load_traj(filename,object='',state=0,format='',interval=1,
                      average=1,start=1,stop=-1,max=-1,selection='all',image=1,
                      shift="[0.0,0.0,0.0]",plugin="",_self=cmd):
        '''
DESCRIPTION

    "load_traj" reads trajectory files.

    Most of the trajectory formats listed here are supported:
    http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/

USAGE

    load_traj filename [,object [,state [,format [,interval [,average ]
                             [,start [,stop [,max [,selection [,image [,shift
                             ]]]]]]]]]

ARGUMENTS

    filename = str: path to trajectory file

    object = str: name of the molecular object where the trajectory should be
    appended as states {default: guess from filename or last object in list}

    state = int: first object state to populate, or 0 to append after
    last state {default: 0}

    format = str: file format {default: guess from extension}

    interval = int: interval to take frames from file {default: 1}

    average = int: ? (trj only, possibly broken)

    start = int: first frame to load from file {default: 1}

    stop = int: last frame to load from file, or -1 to load all {default: -1}

    max = int: maximum number of states to load, or 0 to load all {default: 0}

    selection = str: atom selection to only load a subset of coordinates
    {default: all}

    image = 0/1: residue-based period image transformation (trj only)

    shift = float-3: offset for image transformation {default: (0,0,0}

    plugin = str: name of VMD plugin to use {default: guess from magic string
    of from format}

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

            filename = unquote(filename)

            noext, ext, format_guessed, zipped = filename_to_format(filename)
            fname = _self.exp_path(filename)

            if zipped:
                raise pymol.CmdException('zipped (%s) trajectories not supported' % (zipped))

            if not format:
                format = format_guessed

            if not plugin:
                if format == 'trj':
                    plugin, ftype = check_trj_magic(fname)
                else:
                    try:
                        ftype = int(format)
                    except:
                        plugin = _cmd.find_molfile_plugin(_self._COb, format, 0x2)

    # get object name
            oname = object.strip()
            if not oname:
                oname = _guess_trajectory_object(noext, _self)
                if not len(oname): # safety
                    oname = 'obj01'

            if ftype>=0 or plugin:
                r = _cmd.load_traj(_self._COb,str(oname),fname,int(state)-1,int(ftype),
                                         int(interval),int(average),int(start),
                                         int(stop),int(max),str(selection),
                                         int(image),
                                         float(shift[0]),float(shift[1]),
                                         float(shift[2]),str(plugin))
            else:
                raise pymol.CmdException("unknown format '%s'" % format)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

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
                    if key not in legal_dict:
                        seq_order.append(legal_key)
                    legal_dict[key] = legal_key
                    key = legal_key
                    seq = line[16:].strip()
                    if seq != '':
                        seq_dict[key] = seq_dict.get(key,'') + seq
        for key in seq_order:
            raw_seq = seq_dict[key].replace('-','')
            _self.fab(raw_seq, key, quiet=quiet)

    def _processFASTA(fname, oname, quiet=1, _self=cmd):
        legal_dict = {}
        seq_dict = {}
        seq_order = []
        for line in open(fname).readlines():
            line = line.strip()
            if len(line):
                if line[0:1] == '>':
                    key = line[1:].strip()
                    legal_key = _self.get_legal_name(key)
                    if key not in legal_dict:
                        seq_order.append(legal_key)
                    legal_dict[key] = legal_key
                    key = legal_key
                elif key:
                    seq = line

                    if '-' in seq:
                        # sequence alignment
                        from pymol.seqalign import load_aln_multi
                        return load_aln_multi(fname, oname, quiet=quiet,
                                _self=_self)

                    seq_dict[key] = seq_dict.get(key,'') + seq
        for key in seq_order:
            raw_seq = seq_dict[key].replace('-','')
            _self.fab(raw_seq, key, quiet=quiet)

    def _processPWG(fname,_self=cmd):
        r = DEFAULT_ERROR

        if sys.version_info[0] < 3:
            import urllib
        else:
            import urllib.request as urllib

        import shlex

        try:
            from .pymolhttpd import PymolHttpd
            browser_flag = 0
            launch_flag = 0
            report_url = None
            logging = 1
            root = None
            port = 0
            wrap_native = 0
            headers = []
            if '://' in fname:
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
                        elif keyword == 'header':
                            a = shlex.split(input[1])
                            if len(a) != 3 or a[0] != 'add' or a[1].endswith(':'):
                                raise ValueError('header command must be: header add Some-Key "some value"')
                            headers.append(a[1:])
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
                                    print("Error: requested path '%s' does not exist."%root)
                            else:
                                print("Error: missing path to root content")
                        elif keyword == 'browser':
                            # could perhaps interpret a browser name here
                            browser_flag = 1
                        elif keyword == 'launch': # launch the module named in the file (must exist!)
                            if len(input)>1:
                                mod_name = input[1]
                                try:
                                    __import__(mod_name)
                                    mod = sys.modules[mod_name]
                                    if hasattr(mod,'__launch__'):
                                        mod.__launch__(_self)
                                        r = DEFAULT_SUCCESS
                                except:
                                    traceback.print_exc()
                                    print("Error: unable to launch web application'%s'."%mode_name)
                        elif keyword == 'report':
                            if len(input)>1:
                                report_url = input[1]
                        elif keyword == 'delete':
                            os.unlink(fname)
                        elif keyword == 'options':
                            # parsed during invocation
                            pass
                        elif keyword == 'wrap_native_return_types':
                            wrap_native = 1
                        else:
                            print("Error: unrecognized input:  %s"%str(input))
            if launch_flag:
                server = PymolHttpd(port,root,logging,wrap_native,headers=headers)
                if port == 0:
                    port = server.port # get the dynamically assigned port number
                server.start()
                if browser_flag: # fire up a local browser
                    import webbrowser
                    webbrowser.open("http://localhost:%d"%port,new=1)
                    r = DEFAULT_SUCCESS
                else:
                    r = DEFAULT_SUCCESS
                if report_url is not None: # report port back to server url (is this secure?)
                    try:
                        report_url = report_url + str(port)
                        print(" Reporting back pymol port via: '%s'"%report_url)
                        urllib.urlretrieve(report_url)
                    except:
                        print(" Report attempt may have failed.")
        except ImportError:
            traceback.print_exc()

        if is_error(r):
            print("Error: unable to handle PWG file")
        return r

    def _magic_check_cor_charmm(filename):
        # http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/corplugin.html
        # assume at least 2 title/comment lines, starting with *
        with open(filename) as handle:
            if (handle.readline().startswith('*') and
                handle.readline().startswith('*')):
                return True
        return False


    def _eval_func(func):
        '''
        Evaluate a "module:callable" signature, e.g. "os.path:dirname".
        '''
        if not isinstance(func, str):
            return func

        try:
            m = __import__(func.split(':')[0])
        except ImportError as e:
            print(' Warning: ' + str(e))
            return incentive_format_not_available_func

        return eval(func.replace(':', '.'), {m.__name__: m}, {})


    def load(filename, object='', state=0, format='', finish=1,
             discrete=-1, quiet=1, multiplex=None, zoom=-1, partial=0,
             mimic=1, object_props=None, atom_props=None, _self=cmd):
        '''
DESCRIPTION

    "load" can by used to read molecules, crystallographic maps and
    other volumetric data, PyMOL sessions, and some other types of
    content.

USAGE

    load filename [, object [, state [, format [, finish [, discrete [, quiet
            [, multiplex [, zoom [, partial [, mimic [, object_props
            [, atom_props ]]]]]]]]]]]]

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

    save, load_traj, fetch
        '''
        r = DEFAULT_ERROR
        if object_props or atom_props:
            print(' Warning: properties are not supported in Open-Source PyMOL')
        try:
            _self.lock(_self)
            plugin = ''
            state = int(state)
            finish = int(finish)
            zoom = int(zoom)
            discrete = int(discrete)
            if multiplex is None:
                multiplex=-2

            # file format
            try:
                # user specified the type as an int
                ftype = int(format)
                format = loadable._reverse_lookup(format)
            except ValueError:
                format = str(format)
                ftype = getattr(_loadable, format, -1)

            if ftype in _self._load2str.values():
                assert format.endswith('str')
                colorprinting.warning(
                    ' cmd.load(format="{}") is deprecated, use cmd.load_raw(format="{}")'
                    .format(format, format[:-3]))
                return _self.load_raw(filename, format[:-3], object, state,
                        finish, discrete, quiet, multiplex, zoom)

            filename = unquote(filename)

            # analyze filename
            noext, ext, format_guessed, zipped = filename_to_format(filename)

            if ftype == -1:
                if not format:
                    format = format_guessed
                elif format.startswith('plugin'):
                    format, _, plugin = format.partition(':')
                else:
                    ext = format
                if format == 'pkl':
                    format = 'model' # legacy
                ftype = getattr(_loadable, format, -1)

            filename = _self.exp_path(filename)

            # object name
            object = str(object).strip()
            if not object:
                object = noext if noext else _self.get_unused_name('obj')
                if format in ['dcd', 'dtr']:
                    # for trajectories, use most recently added structure
                    object = _guess_trajectory_object(object, _self)

            # molfile plugins
            if (ftype < 0 and format not in loadfunctions or
                    format == 'plugin' and not plugin):
                plugin = _cmd.find_molfile_plugin(_self._COb, ext)
                if not plugin:
                    raise pymol.CmdException('unsupported file type: ' + ext)
                ftype = loadable.plugin

            # special handling for trj files (autodetect AMBER versus GROMACS)
            if ftype == loadable.trj:
                plugin, ftype = check_trj_magic(filename)
                if plugin:
                    plugin += ':2'  # cPlugIOManager_traj

            # special handling for cdr files (autodetect AMBER versus CHARMM)
            if ftype == loadable.crd and _magic_check_cor_charmm(filename):
                ftype = loadable.plugin
                plugin = 'cor'

            # generic forwarding to format specific load functions
            func = loadfunctions.get(format, pymol.internal._load)
            func = _eval_func(func)
            kw = {
                'filename': filename,
                'fname': filename, # alt
                'object': object,
                'prefix': object, # alt
                'state': state,
                'format': format,
                'finish': finish,
                'discrete': discrete,
                'quiet': quiet,
                'multiplex': multiplex,
                'zoom': zoom,
                'partial': partial,
                'mimic': mimic,
                'object_props': object_props,
                'atom_props': atom_props,
                '_self': _self,

                # for _load
                'ftype': ftype,
                'plugin': plugin,
                'finfo': filename, # alt
                'oname': object, # alt
            }

            import inspect
            spec = inspect.getargspec(func)

            if spec.varargs:
                print('FIXME: loadfunctions[%s]: *args' % (format))

            if not spec.keywords:
                kw = dict((n, kw[n]) for n in spec.args if n in kw)

            if 'contents' in spec.args:
                kw['contents'] = _self.file_read(filename)

            return func(**kw)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def load_pse(filename, partial=0, quiet=1, format='pse', _self=cmd):
        try:
            contents = _self.file_read(filename)
            session = io.pkl.fromString(contents)
        except AttributeError as e:
            raise pymol.CmdException('PSE contains objects which cannot be unpickled (%s)' % str(e))

        r = _self.set_session(session, quiet=quiet, partial=partial, steal=1)

        if not partial:
            _self.set("session_file",
                    # always use unix-like path separators
                    filename.replace("\\", "/"), quiet=1)

        if ((format == 'psw' or
            _self.get_setting_boolean("presentation")) and
            _self.get_setting_boolean("presentation_auto_start")):

            # set movie to first frame
            if _self.get_movie_length():
                _self.rewind()

            # go to first scene
            _self.scene("auto", "start", animate=0)

        return r

    def load_embedded(key=None, name=None, state=0, finish=1, discrete=1,
                      quiet=1, zoom=-1, multiplex=-2, object_props=None,
                      atom_props=None, _self=cmd):
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
        if object_props or atom_props:
            print(' Warning: properties are not supported in Open-Source PyMOL')
        list = _self._parser.get_embedded(key)
        if list is None:
            print("Error: embedded data '%s' not found."%key)
            return DEFAULT_ERROR

        if not name:
            name = key if key else _self.cmd._parser.get_default_key()

        return _self.load_raw(''.join(list[1]), list[0], name, state,
                finish, discrete, quiet, multiplex, zoom)

    def load_raw(content, format, object='', state=0, finish=1,
                 discrete=-1, quiet=1, multiplex=None, zoom=-1,_self=cmd):
        '''
DESCRIPTION

    API-only function for loading data from memory.

EXAMPLE

    contents = open('example.mmtf', 'rb').read()
    cmd.load_raw(contents, 'mmtf')
        '''
        if multiplex is None:
            multiplex=-2
        ftype = getattr(loadable, format, None)
        if ftype is None:
            raise pymol.CmdException("unknown raw format '{}'".format(format))
        with _self.lockcm:
            return _cmd.load(_self._COb, str(object), None, content,
                    int(state) - 1, cmd._load2str.get(ftype, ftype),
                    int(finish), int(discrete),
                                  int(quiet),int(multiplex),int(zoom))

    def read_sdfstr(sdfstr,name,state=0,finish=1,discrete=1,quiet=1,
                    zoom=-1,multiplex=-2,object_props=None,_self=cmd):
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
        if object_props:
            print(' Warning: properties are not supported in Open-Source PyMOL')
        with _self.lockcm:
            return _cmd.load(_self._COb, str(name), None, sdfstr, int(state) - 1,
                              loadable.sdf2str,int(finish),int(discrete),
                              int(quiet),int(multiplex),int(zoom))

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
        with _self.lockcm:
            return _cmd.load(_self._COb, str(name), None, molstr,
                          int(state) - 1,
                          loadable.molstr,int(finish),int(discrete),
                          int(quiet),0,int(zoom))

    def read_mmodstr(content, name, state=0, quiet=1, zoom=-1, _self=cmd, **kw):
        '''
DESCRIPTION

    "read_mmodstr" reads a macromodel format structure from a Python
    string.

    '''
        with _self.lockcm:
            return _cmd.load(_self._COb, str(name), None, content,
                    int(state)-1, loadable.mmodstr, 1, 1, int(quiet), 0,
                    int(zoom))

    def read_pdbstr(contents, oname, state=0, finish=1, discrete=0, quiet=1,
            zoom=-1, multiplex=-2, object_props=None, _self=cmd):
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
        with _self.lockcm:
            r = _cmd.load(_self._COb, str(oname), None, contents,
                    int(state) - 1, loadable.pdbstr,
                              int(finish),int(discrete),int(quiet),
                              int(multiplex),int(zoom))
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
        with _self.lockcm:
            return _cmd.load(_self._COb, str(name), None, mol2,
                    int(state) - 1, loadable.mol2str,
                              int(finish),int(discrete),int(quiet),
                              int(multiplex),int(zoom))


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
        with _self.lockcm:
            return _cmd.load(_self._COb, str(name), None, xplor,
                    int(state) - 1, loadable.xplorstr,
                              int(finish),int(discrete),int(quiet),
                              0,int(zoom))

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

    fetchHosts = {
        "pdb"  : "http://ftp.wwpdb.org/pub/pdb",
        "pdbe" : "ftp://ftp.ebi.ac.uk/pub/databases/pdb",
        "pdbj" : "ftp://ftp.pdbj.org/pub/pdb",
    }

    hostPaths = {
        "mmtf" : "http://mmtf.rcsb.org/v1.0/full/{code}.mmtf.gz",
        "bio"  : [
            "http://files.rcsb.org/download/{code}.{type}.gz",
            "/data/biounit/coordinates/divided/{mid}/{code}.{type}.gz",
        ],
        "pdb"  : [
            "http://files.rcsb.org/download/{code}.{type}.gz",
            "/data/structures/divided/pdb/{mid}/pdb{code}.ent.gz",
        ],
        "cif"  : [
            "http://files.rcsb.org/download/{code}.{type}.gz",
            "/data/structures/divided/mmCIF/{mid}/{code}.cif.gz",
            "http://ftp-versioned.wwpdb.org/pdb_versioned/views/latest/coordinates/mmcif/{mid}/pdb_{code:0>8}/pdb_{code:0>8}_xyz.cif.gz",
        ],
        "2fofc" : "https://www.ebi.ac.uk/pdbe/coordinates/files/{code}.ccp4",
        "fofc" : "https://www.ebi.ac.uk/pdbe/coordinates/files/{code}_diff.ccp4",
        "pubchem": [
            "http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?{type}={code}&disopt=3DSaveSDF",
            "http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?{type}={code}&disopt=SaveSDF",
        ],
        "emd": "ftp://ftp.wwpdb.org/pub/emdb/structures/EMD-{code}/map/emd_{code}.map.gz",
        "cc": [
            "http://files.rcsb.org/ligands/download/{code}.cif",
            "ftp://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/files/mmcif/{code}.cif",
        ],
    }

    def _fetch(code, name, state, finish, discrete, multiplex, zoom, type, path,
            file, quiet, _self=cmd):
        '''
        code = str: single pdb identifier
        name = str: object name
        state = int: object state
        finish =
        discrete = bool: make discrete multi-state object
        multiplex = bool: split states into objects (like split_states)
        zoom = int: zoom to new loaded object
        type = str: fofc, 2fofc, pdb, pdb1, ... 
        path = str: fetch_path
        file = str or file: file name or open file handle
        '''
        r = DEFAULT_ERROR

        fetch_host_list = [x if '://' in x else fetchHosts[x]
                for x in _self.get("fetch_host").split()]

        # file types can be: fofc, 2fofc, pdb, pdb1, pdb2, pdb3, etc...
        # bioType is the string representation of the type
        # nameFmt is the file name pattern after download
        bioType = type
        nameFmt = '{code}.{type}'
        if type == 'pdb':
            pass
        elif type in ('fofc', '2fofc'):
            nameFmt = '{code}_{type}.ccp4'
        elif type == 'emd':
            nameFmt = '{type}_{code}.ccp4'
        elif type in ('cid', 'sid'):
            bioType = 'pubchem'
            nameFmt = '{type}_{code}.sdf'
        elif type == 'cif':
            pass
        elif type == 'mmtf':
            pass
        elif type == 'cc':
            nameFmt = '{code}.cif'
        elif re.match(r'pdb\d+$', type):
            bioType = 'bio'
        else:
            raise ValueError('type')

        url = hostPaths[bioType]
        url_list = []
        for url in url if cmd.is_sequence(url) else [url]:
            url_list += [url] if '://' in url else [fetch_host + url
                for fetch_host in fetch_host_list]

        if bioType not in ['cc']:
            code = code.lower()

        fobj = None
        contents = None

        if not file or file in (1, '1', 'auto'):
            file = os.path.join(path, nameFmt.format(code=code, type=type))

        if not is_string(file):
            fobj = file
            file = None
        elif os.path.exists(file):
            # skip downloading
            url_list = []

        for url in url_list:
            url = url.format(mid=code[-3:-1], code=code, type=type)

            try:
                contents = _self.file_read(url)

                # assume HTML content means error on server side without error HTTP code
                if b'<html' in contents[:500].lower():
                    raise pymol.CmdException

            except pymol.CmdException:
                if not quiet:
                    colorprinting.warning(" Warning: failed to fetch from %s" % (url,))
                continue

            if file:
                try:
                    fobj = open(file, 'wb')
                except IOError:
                    colorprinting.warning(' Warning: Cannot write to "%s"' % file)

            if fobj:
                fobj.write(contents)
                fobj.flush()
                if file:
                    fobj.close()

            if not file:
                return DEFAULT_SUCCESS

            break

        if os.path.exists(file):
            r = _self.load(file, name, state, '',
                    finish, discrete, quiet, multiplex, zoom)
        elif contents and bioType in ('pdb', 'bio'):
            r = _self.read_pdbstr(contents, name, state,
                    finish, discrete, quiet, zoom, multiplex)
        elif contents and bioType in ('cif', 'cc'):
            r = _self.load_raw(contents, 'cif', name, state,
                    finish, discrete, quiet, multiplex, zoom)
        elif contents and bioType in ('mmtf',):
            r = _self.load_raw(contents, 'mmtf', name, state,
                    finish, discrete, quiet, multiplex, zoom)

        if not _self.is_error(r):
            return name

        colorprinting.error(" Error-fetch: unable to load '%s'." % code)
        return DEFAULT_ERROR

    def _multifetch(code,name,state,finish,discrete,multiplex,zoom,type,path,file,quiet,_self):
        import string
        r = DEFAULT_SUCCESS
        code_list = code.split()
        name = name.strip()
        if (name!='') and (len(code_list)>1) and (discrete<0):
            discrete = 1 # by default, select discrete  when loading
            # multiple PDB entries into a single object

        all_type = type
        for obj_code in code_list:
            obj_name = name
            type = all_type

            if not type:
                if 1 < len(obj_code) < 4:
                    type = 'cc'
                else:
                    type = _self.get('fetch_type_default')

            # allow fetching codes like EMD-3489 or emd_3489
            if obj_code[3:4] in ('_', '-') and \
                    obj_code[:3].upper() in ('CID', 'SID', 'EMD'):
                if not obj_name:
                    obj_name = obj_code
                type = obj_code[:3].lower()
                obj_code = obj_code[4:]

            if not obj_name:
                obj_name = obj_code
                if type.endswith('fofc'):
                    obj_name += '_' + type
                elif type == 'emd':
                    obj_name = 'emd_' + obj_code

            chain = None
            if len(obj_code) in (5,6,7) and type in ('pdb', 'cif', 'mmtf'):
                obj_code = (
                    obj_code
                        .replace('.', '')
                        .replace('_', '')
                        .replace('-', '')
                        .replace(':', '')
                    )
                obj_code, chain = obj_code[:4], obj_code[4:]

            obj_name = _self.get_legal_name(obj_name)

            r = _fetch(obj_code, obj_name, state, finish,
                    discrete, multiplex, zoom, type, path, file, quiet, _self)

            if chain and isinstance(r, str):
                if _self.count_atoms(r'?%s & c. \%s' % (r, chain)) == 0:
                    _self.delete(r)
                    raise pymol.CmdException('no such chain: ' + chain)
                _self.remove(r'?%s & ! c. \%s' % (r, chain))

        return r

    def fetch(code, name='', state=0, finish=1, discrete=-1,
              multiplex=-2, zoom=-1, type='', async_=0, path='',
              file=None, quiet=1, _self=cmd, **kwargs):

        '''
DESCRIPTION

    "fetch" downloads a file from the internet (if possible)

USAGE

    fetch code [, name [, state [, finish [, discrete [, multiplex
        [, zoom [, type [, async [, path ]]]]]]]]]

ARGUMENTS

    code = a single PDB identifier or a list of identifiers. Supports
    5-letter codes for fetching single chains (like 1a00A).

    name = the object name into which the file should be loaded.

    state = the state number into which the file should loaded.

    type = str: cif, pdb, pdb1, 2fofc, fofc, emd, cid, sid {default: cif
    (default was "pdb" up to 1.7.6)}

    async_ = 0/1: download in the background and do not block the PyMOL
    command line {default: 0 -- changed in PyMOL 2.3}

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
        state, finish, discrete = int(state), int(finish), int(discrete)
        multiplex, zoom = int(multiplex), int(zoom)
        async_, quiet = int(kwargs.pop('async', async_)), int(quiet)

        if kwargs:
            raise pymol.CmdException('unknown argument: ' + ', '.join(kwargs))

        r = DEFAULT_SUCCESS
        if not path:
            # blank paths need to be reset to '.'
            path = _self.get('fetch_path') or '.'
        if async_ < 0: # by default, run asynch when interactive, sync when not
            async_ = not quiet
        args = (code, name, state, finish, discrete, multiplex, zoom, type, path, file, quiet, _self)
        kwargs = { '_self' : _self }
        if async_:
            _self.async_(_multifetch, *args, **kwargs)
        else:
            try:
                _self.block_flush(_self)
                r = _multifetch(*args)
            finally:
                _self.unblock_flush(_self)
        return r

    def load_coordset(coords, object, state=0, quiet=1, _self=cmd):
        '''
DESCRIPTION

    API only. Load object coordinates. Loads them in the original atom
    order (order from PDB file for example), not in the atom property
    sorted order (like cmd.iterate, cmd.load_coods, etc.).

ARGUMENTS

    coords = list: Nx3 float array

    object = str: object name

    state = int: object state, or 0 for append {default: 0}

SEE ALSO

    cmd.load_coords
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.load_coordset(_self._COb, object, coords, int(state)-1)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def load_coords(coords, selection, state=1, quiet=1, _self=cmd):
        '''
DESCRIPTION

    API only. Load selection coordinates.

    CHANGED IN VERSION 1.7.3: This used to be the load_coordset function.
    load_coordset may load coordinates in different order (original order from
    PDB file) than load_coords (atom sorted order).

ARGUMENTS

    coords = list: Nx3 float array

    selection = str: atom selection

    state = int: object state {default: 1}
        '''
        with _self.lockcm:
            r = _cmd.load_coords(_self._COb, selection, coords, int(state)-1)
        return r

    def load_idx(filename, object, state=0, quiet=1, zoom=-1, _self=cmd):
        '''
DESCRIPTION

    Load a Desmond trajectory (including topology) from an IDX file.
        '''
        state, quiet, zoom = int(state), int(quiet), int(zoom)
        dirname = os.path.dirname(filename)
        data = {}

        with open(filename) as handle:
            for line in handle:
                k, v = line.split('=', 1)
                data[k.strip()] = v.strip()

        oname = object
        _self.delete(oname)
        r = _self.load(
                os.path.join(dirname, data['structure']), oname, state, quiet=quiet, zoom=0)

        if not is_error(r):
            _self.load_traj(
                os.path.join(dirname, data['trajectory'], 'clickme.dtr'), oname, state)

            auto_zoom(zoom, oname, _self=_self)

        return r

    def load_mtz(filename, prefix='', amplitudes='', phases='', weights='None',
            reso_low=0, reso_high=0, quiet=1, _self=cmd):
        '''
DESCRIPTION

    Load a MTZ file as two map objects (fofc, 2fofc) or if amplitudes and
    phases column names are given, as one map object.

USAGE

    load_mtz filename [, prefix [, amplitudes, phases [, weights
        [, reso_low [, reso_high ]]]]]

ARGUMENTS

    filename = str: filename

    prefix = str: object name or prefix {default: filename without extension}

    amplitudes = str: amplitudes column name, guess if blank {default: }

    phases = str: phases column name, required if amplitudes are given {default: }

    weights = str: weights column name, optional {default: None}

    reso_low = float: minimum resolution {default: 0, read from file}

    reso_high = float: maximum resolution {default: 0, read from file}

        '''
        raise pymol.IncentiveOnlyException()

    def loadall(pattern, group='', quiet=1, _self=cmd, **kwargs):
        '''
DESCRIPTION

    Load all files matching given globbing pattern

USAGE

    loadall pattern [, group ]

EXAMPLE

    loadall *.pdb
        '''
        import glob

        filenames = glob.glob(_self.exp_path(pattern))

        for filename in filenames:
            if not quiet:
                print(' Loading', filename)
            _self.load(filename, **kwargs)

        if group:
            if kwargs.get('object', '') != '':
                print(' Warning: group and object arguments given')
                members = [kwargs['object']]
            else:
                members = map(filename_to_objectname, filenames)
            _self.group(group, ' '.join(members))


    def load_mmtf(filename, object='', discrete=0, multiplex=0, zoom=-1, quiet=1, _self=cmd):
        '''
DESCRIPTION

    Load an MMTF file or URL.
        '''
        from chempy.mmtf import MmtfReader

        if not object:
            object = filename_to_objectname(filename)

        data = MmtfReader.from_url(filename)
        models = data.to_chempy(_self.get_setting_int('cif_use_auth'))

        if len(models) == 1:
            _self.load_model(models[0], object, discrete=discrete, zoom=zoom, quiet=quiet)

            assembly = _self.get('assembly').encode()
            if not assembly:
                return

            try:
                transformList = next((a[b'transformList']
                    for a in data.get(b'bioAssemblyList', ())
                    if a[b'name'] == assembly))
            except StopIteration:
                raise pymol.CmdException('No such assembly: "%s"' % (assembly))

            chainIdList = data.get(b'chainIdList')

            tmp = _self.get_unused_name('_asu')
            _self.set_name(object, tmp)

            for state, trans in enumerate(transformList, 1):
                chains = [chainIdList[i] for i in trans[b'chainIndexList']]
                _self.create(object,
                        'model %s and segi %s' % (tmp, '+'.join(chains)),
                        1, state, zoom=0, discrete=discrete)
                _self.transform_object(object, trans[b'matrix'], state=state)

            _self.delete(tmp)
            _self.set('all_states', 1, object)

            if int(zoom):
                _self.zoom(object)

            return

        tmp = _self.get_unused_name('_model')

        for i, model in enumerate(models):
            _self.load_model(model, tmp, 1, zoom=0)
            if multiplex > 0:
                _self.set_name(tmp, '%s_%04d' % object)
            else:
                _self.create(object, tmp, 1, -1, discrete=discrete, zoom=zoom)
                _self.delete(tmp)

    def load_ply(filename, object, state=0, zoom=-1, _self=cmd):
        from . import cgo
        obj = cgo.from_plystr(_self.file_read(filename))
        r = _self.load_cgo(obj, object, state, zoom=zoom)
        _self.set('cgo_lighting', 1, object)
        return r

    def load_r3d(filename, object, state=0, zoom=-1, _self=cmd):
        from . import cgo
        obj = cgo.from_r3d(filename)
        return _self.load_cgo(obj, object, state, zoom=zoom)

    def load_cc1(filename, object, state=0, _self=cmd):
        obj = io.cc1.fromFile(filename)
        return _self.load_model(obj, object, state)

    loadfunctions = {
        'mae': incentive_format_not_available_func,
        'pdbml': 'pymol.lazyio:load_pdbml',
        'cml': 'pymol.lazyio:load_cml',
        'mtz': load_mtz,
        'py': lambda filename, _self: _self.do("_ run %s" % filename),
        'pml': lambda filename, _self: _self.do("_ @%s" % filename),
        'pwg': _processPWG,
        'aln': 'pymol.seqalign:load_aln_multi',
        'fasta': _processFASTA,
        'png': 'pymol.viewing:load_png',
        'idx': load_idx,
        'pse': load_pse,
        'psw': load_pse,
        'ply': load_ply,
        'r3d': load_r3d,
        'cc1': load_cc1,
        'pdb': read_pdbstr,
        'stl': 'pymol.lazyio:read_stlstr',

        # Incentive
        'vis': incentive_format_not_available_func,
        'moe': incentive_format_not_available_func,
        'phypo': incentive_format_not_available_func,
    }
