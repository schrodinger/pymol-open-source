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

if __name__=='pymol.exporting':
    import os
    import sys
    if sys.version_info[0] == 2:
        import thread
        import cPickle
    else:
        import _thread as thread
        import pickle as cPickle
    from . import selector
    import re
    import copy
    
    import pymol
    cmd = sys.modules["pymol.cmd"]
    from .cmd import _cmd,lock,unlock,Shortcut,QuietException
    from chempy import io
    from chempy.sdf import SDF,SDFRec
    from .cmd import _feedback,fb_module,fb_mask, \
                     DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error, \
                     is_list, is_dict, is_tuple, loadable
    import traceback

    def copy_image(quiet=1,_self=cmd): # incentive feature / proprietary
        r = DEFAULT_ERROR
        if thread.get_ident() == pymol.glutThread:
            r = _self._copy_image(_self,int(quiet))
        else:
            r = _self.do('cmd._copy_image(quiet=%d,_self=cmd)'%int(quiet))
        if _self._raising(r,_self): raise QuietException         
        return r

    cache_action_dict = {
        'enable'      : 0,
        'disable'     : 1,
        'read_only'   : 2,
        'clear'       : 3,
        'optimize'    : 4,
    }

    cache_action_sc = Shortcut(cache_action_dict.keys())

    def cache(action='optimize', scenes='',state=-1, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "cache" manages storage of precomputed results, such as
    molecular surfaces.

USAGE

    cache action [, scenes [, state ]]

ARGUMENTS

    action = string: enable, disable, read_only, clear, or optimize

    scenes = string: a space-separated list of scene names (default: '')

    state = integer: state index (default: -1)

EXAMPLES

    cache enable
    cache optimize
    cache optimize, F1 F2 F5

NOTES

    "cache optimize" will iterate through the list of scenes provided
    (or all defined scenes), compute any missing surfaces, and store
    them in the cache for later reuse.
    
PYMOL API

    cmd.cache(string action, string scenes, int state, int quiet)

    '''
        
        r = DEFAULT_ERROR
        action = cache_action_dict[cache_action_sc.auto_err(str(action),'action')]
        quiet = int(quiet)
        if action == 0: # enable
            r = _self.set('cache_mode',2,quiet=quiet)
        elif action == 1: # disable
            r =_self.set('cache_mode',0,quiet=quiet)
        elif action == 2: # read_only
            r =_self.set('cache_mode',1,quiet=quiet)
        elif action == 3: # clear
            r =_self._cache_clear(_self=_self)
        elif action == 4: # optimize
            r = DEFAULT_SUCCESS
            _self._cache_mark(_self=_self)
            cur_scene = _self.get('scene_current_name')
            cache_max = int(_self.get('cache_max'))
            if cache_max>0:
                # allow double memory for an optimized cache
                _self.set('cache_max',cache_max*2) 
            scenes = str(scenes)
            scene_list = scenes.split()
            cache_mode = int(_self.get('cache_mode'))
            _self.set('cache_mode',2)
            if not len(scene_list):
                scene_list = _self.get_scene_list()
            for scene in scene_list:
                scene = scene.strip()
                if not quiet:
                    print(" cache: optimizing scene '%s'."%scene)
                cmd.scene(scene,animate=0)
                cmd.rebuild()                
                cmd.refresh()
            if len(cur_scene):
                cmd.scene(cur_scene,animate=0)
            else:
                scene_list = _self.get_scene_list()
                if len(scene_list):
                    cmd.scene(scene_list[0],animate=0)
                else:
                    if not quiet:
                        print(" cache: no scenes defined -- optimizing current display.")
                    cmd.rebuild() 
                    cmd.refresh()
            usage = _self._cache_purge(-1,_self=_self)
            if cache_mode:
                _self.set('cache_mode',cache_mode)
            else:
                _self.set('cache_mode',2) # hmm... could use 1 here instead.
            _self.set('cache_max',cache_max) # restore previous limits
            if not quiet:            
                print(" cache: optimization complete (~%0.1f MB)."%(usage*4/1000000.0))
        try:
            _self.lock(_self)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException         
        return r

    _resn_to_aa =  {
            'ALA' : 'A',
            'CYS' : 'C',
            'ASP' : 'D',
            'GLU' : 'E',
            'PHE' : 'F',
            'GLY' : 'G',
            'HIS' : 'H',
            'ILE' : 'I',
            'LYS' : 'K',
            'LEU' : 'L',
            'MET' : 'M',
            'ASN' : 'N',
            'PRO' : 'P',
            'GLN' : 'Q',
            'ARG' : 'R',
            'SER' : 'S',
            'THR' : 'T',
            'VAL' : 'V',
            'TRP' : 'W',
            'TYR' : 'Y',
            }
    
    def get_fastastr(selection="all", state=-1, quiet=1, _self=cmd):
        dict = { 'seq' : {} }
        # we use (alt '' or alt 'A') because 'guide' picks up 
        # non-canonical structures: eg, 1ejg has residue 22 as a SER and 
        # PRO, which guide will report twice
        _self.iterate("("+selection+") and polymer and name CA and alt +A",
                    "seq[model]=seq.get(model,[]);seq[model].append(resn)",space=dict)
        seq = dict['seq']
        result = []
        for obj in _self.get_names("objects",selection='('+selection+')'):
            if obj in seq:
                cur_seq = [_resn_to_aa.get(x,'?') for x in seq[obj]]
                result.append(">%s"%obj)
                cur_seq = ''.join(cur_seq)
                while len(cur_seq):
                    if len(cur_seq)>=70:
                        result.append(cur_seq[0:70])
                        cur_seq=cur_seq[70:]
                    else:
                        result.append(cur_seq)
                        break
        result = '\n'.join(result)
        if len(result):
            result = result + '\n'
        return result

    def get_pdbstr(selection="all", state=-1, ref='', ref_state=-1, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "get_pdbstr" in an API-only function which returns a pdb
    corresponding to the atoms in the selection provided and that are
    present in the indicated state

PYMOL API

    cmd.get_pdbstr(string selection, int state)

NOTES

    "state" is a 1-based state index for the object.

    if state is -1, then current state is used.

    if state is 0, then all states are saved.
    
    '''
        return get_str('pdb', selection, state, ref, ref_state, -1, quiet, _self)

    def _get_dump_str(obj):
        if is_list(obj):
            list = map(_get_dump_str,obj)
            result = "[ " + ",\n".join(list) + " ] "
        elif is_dict(obj):
            list = []
            for key in obj.keys():
                list.append( _get_dump_str(key)+" : "+_get_dump_str(obj[key]) )
            result = "{ " + ",\n".join(list) + " } "
        elif is_tuple(obj):
            list = map(_get_dump_str,obj)
            result = "( " + ",\n".join(list) + " ) "
        else:
            result = str(obj)
        return result
    
    def _session_convert_legacy(session, version, _self=cmd):
        '''
        Convert session contents to be compatible with previous PyMOL versions
        '''
        if version >= _self.get_version()[1]:
            return

        print(" Applying pse_export_version=%.3f compatibility" % (version))

        def bitmaskToList(mask):
            if not isinstance(mask, int):
                return mask
            r = []
            while (mask >> len(r)):
                r.append(1 if ((mask >> len(r)) & 1) else 0)
            return r

        def convert_settings(settings_list):
            if not settings_list:
                return

            # get index -> setting dict
            settings = dict((s[0], s) for s in settings_list if s is not None)

            # changed type
            cast_fn = {
                3: float,
                4: lambda v: list(_self.get_color_tuple(v)),
            }
            for (i, old, new) in [
                    (6, 4, 5),          # bg_rgb
                    (254, 3, 1),        # scenes_changed
                    ]:
                if i in settings and settings[i][1] == new:
                    settings[i][1] = old
                    settings[i][2] = cast_fn[old](settings[i][2])

            # changed defaults
            changed = []
            if version < 1.76:
                changed += [
                    (91, 7, -1),        # cartoon_sampling
                    (93, 6., -1.),      # cartoon_loop_quality
                    (102, 10., -1.),    # cartoon_oval_quality
                    (104, 9., -1.),     # cartoon_tube_quality
                    (378, 11., -1.),    # cartoon_putty_quality
                ]
            if version < 1.5:
                changed += [
                    (421, -1, 9),       # sphere_mode
                ]
            for (i, old, new) in changed:
                if i in settings and settings[i][2] == new:
                    settings[i][2] = old

        if 'settings' in session:
            convert_settings(session['settings'])

        # objects
        for name in  session['names']:
            if name is None:
                continue

            # spec repOn
            if version < 1.76 and name[3] is None:
                name[3] = [0] * 20

            # only continue for objects
            if name[1] != 0:
                continue

            if name[4] == 8: # gadget
                cobject = name[5][0][0]
            else:
                cobject = name[5][0]

            # object visRep
            if version < 1.76:
                cobject[3] = bitmaskToList(cobject[3])

            # object settings
            convert_settings(cobject[8])

            # molecule
            if name[4] == 1:
                # atoms visRep
                if version < 1.76:
                    for atom in name[5][7]:
                        atom[20] = bitmaskToList(atom[20])
                # state settings
                for cset in name[5][4]:
                    if cset:
                        convert_settings(cset[7])

            # map
            elif name[4] == 2:
                if version < 1.5:
                    # crystal -> [crystal, spacegroup]
                    for state in name[5][2]:
                        if state and state[1]:
                            state[1] = state[1][0]

            # cgo
            elif name[4] == 6:
                if version < 1.506:
                    for state in name[5][2]:
                        if len(state) == 1:
                            # std and ray CGO
                            state.append(state[0])

    def get_session(names='', partial=0, quiet=1, compress=-1, cache=-1, _self=cmd):
        session = {}
        r = DEFAULT_SUCCESS
        cache = int(cache)
        compress = int(compress)
        partial = int(partial)

        pse_export_version = round(_self.get_setting_float('pse_export_version'), 4)
        legacyscenes = (0 < pse_export_version < 1.76) and _self.get_scene_list()

        if legacyscenes:
            _self.pymol._scene_dict = {}

            scene_current_name = _self.get('scene_current_name')
            tempname = '_scene_db96724c3cef00875c3bebb4f348f711'
            _self.scene(tempname, 'store')

            for name in legacyscenes:
                _self.scene(name, 'recall', animate=0)
                wizard = _self.get_wizard()
                message = wizard.message if getattr(wizard, 'from_scene', 0) else None
                pymol.viewing._legacy_scene(name, 'store', message, _self=_self)

            _self.scene(tempname, 'recall', animate=0)
            _self.scene(tempname, 'clear')
            _self.set('scene_current_name', scene_current_name)

        if cache:
            cache_opt = int(_self.get('session_cache_optimize'))
            if cache != 0:
                cache_mode = int(_self.get('cache_mode'))
                if ((cache_mode > 0) and (cache_opt != 0)) or (cache_opt==1):
                    _self.cache('optimize')
        for a in _self._pymol._session_save_tasks:
            if a==None:
                try:
                    _self.lock(_self)
                    r = _cmd.get_session(_self._COb,session,str(names),
                                         int(partial),int(quiet))
                finally:
                    _self.unlock(r,_self)
                try:
                    session['session'] = copy.deepcopy(_self._pymol.session)
                    if cache and hasattr(_self._pymol,'_cache'):
                        session['cache'] = _self._pymol._cache
                except:
                    traceback.print_exc()
            else:
                try:
                    if is_error(a(*(session,), **{'_self':_self})):
                        r = DEFAULT_ERROR
                except:
                    traceback.print_exc()
                    print("Error: An error occurred when trying to generate session.")
                    print("Error: The resulting session file may be incomplete.")

        if legacyscenes:
            del session['moviescenes']
            session['scene_dict'] = _self.pymol._scene_dict
            session['scene_order'] = legacyscenes
            pymol.viewing._legacy_scene('*', 'clear', _self=_self)
            del _self.pymol._scene_dict

        if is_ok(r):
            if pse_export_version > 0.0:
                try:
                    _session_convert_legacy(session, pse_export_version, _self)
                except Exception as e:
                    print(' Warning: failed to backport session:', e)

            if(compress<0):
                compress = _self.get_setting_boolean('session_compression')
            if(compress):
                import zlib
                session = zlib.compress(cPickle.dumps(session, 1))
            return session
        elif _self._raising(r,_self):
            raise QuietException                  
        return r
        
    def _unit2px(value, dpi, unit=''):
        '''API only. Returns pixel units given a string representation in other units'''
        if cmd.is_string(value):
            m = re.search(r'[a-z].*', value, re.I)
            if m:
                value, unit = value[:m.start()], m.group(0).lower()

        if not unit or unit == 'px':
            return int(value)

        upi = {'in': 1.0, 'mm': 25.4, 'cm': 2.54}
        if unit not in upi:
            raise pymol.CmdException('unknown unit, supported units are: ' +
                    ', '.join(upi))

        if dpi < 1:
            raise pymol.CmdException('dpi > 0 required with unit "%s" '
                    '(hint: set the "image_dots_per_inch" setting)' % unit)

        return float(value) * dpi / upi[unit] + 0.5

    def png(filename, width=0, height=0, dpi=-1.0, ray=0,
            quiet=1, prior=0, format=0, _self=cmd):
        '''
DESCRIPTION

    "png" saves a PNG format image file of the current display.

USAGE

    png filename [, width [, height [, dpi [, ray]]]]

ARGUMENTS

    filename = string: file path to be written
    
    width = integer or string: width in pixels (without units), inches (in)
    or centimeters (cm). If unit suffix is given, dpi argument is required
    as well. If only one of width or height is given, the aspect ratio of
    the viewport is preserved. {default: 0 (current)}

    height = integer or string: height (see width) {default: 0 (current)}

    dpi = float: dots-per-inch {default -1.0 (unspecified)}

    ray = 0 or 1: should ray be run first {default: 0 (no)}

EXAMPLES

    png image.png
    png image.png, dpi=300
    png image.png, 10cm, dpi=300, ray=1

NOTES

    PNG is the only image format supported by PyMOL.

SEE ALSO

    mpng, save
    
PYMOL API

    cmd.png(string filename, int width, int height, float dpi,
            int ray, int quiet)
        '''
        r = DEFAULT_ERROR
        prior = int(prior)

        if format == 'png':
            format = 0

        if prior: 
            # fetch the prior image, without doing any work (fast-path / non-GLUT thread-safe)
            r = _self._png(str(filename),0,0,float(dpi),0,int(quiet),1,
                           int(format),_self)            
            if r != 1: # no prior image available -- revert to default behavior
                if prior < 0: # default is to fall back to actual rendering
                    prior = 0
        if not prior:
            dpi = float(dpi)
            if dpi < 0:
                dpi = cmd.get_setting_float('image_dots_per_inch')
            width = _unit2px(width, dpi)
            height = _unit2px(height, dpi)

            if thread.get_ident() == pymol.glutThread:
                r = _self._png(str(filename),int(width),int(height),float(dpi),
                               int(ray),int(quiet),0,int(format),_self)
            else:
                r = _self._do("cmd._png('''%s''',%d,%d,%1.6f,%d,%d,%d,%d)"%
                              (filename,width,height,dpi,
                               ray,int(quiet),0,int(format)),_self=_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def export_coords(obj,state,_self=cmd): # experimental
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.export_coords(_self._COb,str(obj),int(state)-1)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def multisave(filename, pattern="all", state=-1,
                  append=0, format='', quiet=1, _self=cmd): 
        '''
DESCRIPTION

    "multisave" will save a multi-entry PDB file.

    Every object in the given selection (pattern) will have a HEADER and a
    CRYST (if symmetry is defined) record, and is terminated with END.
    Loading such a multi-entry PDB file into PyMOL will load each entry
    as a separate object.

    This behavior is different to the "save" command, where a multi-object
    selection is written "flat" to a PDB file, without HEADER or CRYST
    records.

ARGUMENTS

    filename = string: file path to be written

    pattern = str: atom selection (before 1.8.4: object name pattern)

    state = int: object state (-1=current, 0=all) {default: -1}

    append = 0/1: append to existing file {default: 0}

    format = str: file format {default: guess from extension, or 'pdb'}
    '''
        from pymol.importing import filename_to_format
        _, _, format_guessed, zipped = filename_to_format(filename)

        if zipped:
            raise pymol.CmdException(zipped + ' not supported with multisave')

        if not format:
            format = format_guessed or 'pdb'

        if format == 'pmo':
            raise pymol.CmdException('pmo format not supported anymore')

        if format not in ('pdb', 'cif'):
            raise pymol.CmdException(format + ' format not supported with multisave')

        s = get_str(format, pattern, state, '', -1, 1, quiet, _self)

        if s is None:
            if _self._raising(): raise QuietException
            return DEFAULT_ERROR

        filename = _self.exp_path(filename)

        with open(filename, 'a' if int(append) else 'w') as handle:
            handle.write(s)

        return DEFAULT_SUCCESS

    def assign_atom_types( selection, format = "mol2", state=1, quiet=1, _self=cmd):
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            # format : mol2/sybyl = 1, macromodel/mmd = 2, global setting atom_type_format = 0
            r = _cmd.assign_atom_types(_self._COb, selection, int(1), int(state-1), quiet)
        finally:
            _self.unlock(r,_self)
        return r

    def get_str(format, selection='(all)', state=-1, ref='',
             ref_state=-1, multi=-1, quiet=1, _self=cmd):
        '''
DESCRIPTION

    API-only function which exports the selection to a molecular file
    format and returns it as a string.

ARGUMENTS

    format = str: pdb|cif|sdf|mol|mol2|mae|pqr|xyz

    selection = str: atom selection {default: all}

    state = int: object state (-1=current, 0=all) {default: -1}

    ref = str: object name which defines reference frame {default: }

    ref_state = int: state of ref object {default: -1}

    multi = int: for multi-entry file formats, 0 = single entry,
    1 = by object, 2 = by object-state, -1 = format default {default: -1}
        '''
        with _self.lockcm:
            return _cmd.get_str(_self._COb, str(format), str(selection),
                    int(state) - 1, str(ref), int(ref_state),
                    int(multi), int(quiet))

    def save(filename, selection='(all)', state=-1, format='', ref='',
             ref_state=-1, quiet=1, partial=0,_self=cmd):
        '''
DESCRIPTION

    "save" writes content to a file.
    
USAGE

    save filename [, selection [, state [, format ]]]

ARGUMENTS

    filename = string: file path to be written

    selection = string: atoms to save {default: (all)}

    state = integer: state to save {default: -1 (current state)}
    
PYMOL API

    cmd.save(string file, string selection, int state, string format)

NOTES

    The file format is automatically chosen if the extesion is one of
    the supported output formats: pdb, pqr, mol, sdf, pkl, pkla, mmd, out,
    dat, mmod, cif, pov, png, pse, psw, aln, fasta, obj, mtl, wrl, dae, idtf,
    or mol2.

    If the file format is not recognized, then a PDB file is written
    by default.

    For molecular files and where applicable and supported:
    
    * if state = -1 (default), then only the current state is written.

    * if state = 0, then a multi-state output file is written.
    
SEE ALSO

    load, get_model
        '''
        quiet = int(quiet)

        # preprocess selection
        selection = selector.process(selection)
        #   
        r = DEFAULT_ERROR

        # analyze filename
        from pymol.importing import filename_to_format
        _, _, format_guessed, zipped = filename_to_format(filename)
        filename = _self.exp_path(filename)

        # file format
        if not format:
            if not format_guessed:
                raise pymol.CmdException('Unrecognized file format')
            format = format_guessed

        # PyMOL session
        if format in ('pse', 'psw',):
            _self.set("session_file",
                    # always use unix-like path separators
                    filename.replace("\\", "/"), quiet=1)
            if not quiet:
                print(" Save: Please wait -- writing session file...")

        func_type4 = {
            'mmod': io.mmd.toFile,
            'pkl': io.pkl.toFile, # binary pickle
            'pkla': lambda model, filename: io.pkl.toFile(model, filename, bin=0), # ascii pickle
        }

        contents = None

        if format in savefunctions:
            # generic forwarding to format specific save functions
            func = savefunctions[format]
            kw = {
                'filename': filename,
                'selection': selection,
                'name': selection,      # alt (get_ccp4str)
                'state': state,
                'format': format,
                'ref': ref,
                'ref_state': ref_state,
                'quiet': quiet,
                'partial': partial,
                '_self': _self,
            }

            import inspect
            spec = inspect.getargspec(func)

            if spec.varargs:
                print('FIXME: savefunctions[%s]: *args' % (format))

            if not spec.keywords:
                kw = dict((n, kw[n]) for n in spec.args if n in kw)

            contents = func(**kw)

            if 'filename' in spec.args:
                # assume function wrote directly to file and returned a status
                return contents

        elif format in func_type4:
            func_type4[format](_self.get_model(selection, state, ref, ref_state), filename)
            r = DEFAULT_SUCCESS
        else:
            raise pymol.CmdException('File format not supported for export')

        # function returned sequence of strings or bytes
        if isinstance(contents, (tuple, list)) and contents:
            contents = contents[0][:0].join(contents)

        if cmd.is_string(contents):
            if not isinstance(contents, bytes):
                contents = contents.encode()

            if zipped == 'gz':
                import gzip
                fopen = gzip.open
            else:
                fopen = open
                if zipped == 'bz2':
                    import bz2
                    contents = bz2.compress(contents)

            with fopen(filename, 'wb') as handle:
                handle.write(contents)
            r = DEFAULT_SUCCESS
            
        if _self._raising(r,_self): raise QuietException

        if not quiet:
            if r == DEFAULT_SUCCESS:
                print(' Save: wrote "' + filename + '".')
            else:
                print(' Save-Error: no file written')

        return r

    def get_cifstr(selection="all", state=-1, quiet=1, _self=cmd):
        '''
DESCRIPTION

    API-only function which returns a mmCIF string.

SEE ALSO

    get_pdbstr
        '''
        return get_str('cif', selection, state, '', -1, -1, quiet, _self)

    def get_xyzstr(selection, state=-1, quiet=1, _self=cmd):
        return get_str('xyz', selection, state, '', -1, -1, quiet, _self)

    def get_sdfstr(selection, state=-1, ref='', ref_state=-1, quiet=1, _self=cmd):
        return get_str('sdf', selection, state, ref, ref_state, -1, quiet, _self)

    def get_mol2str(selection, state=-1, ref='', ref_state=-1, quiet=1, _self=cmd):
        return get_str('mol2', selection, state, ref, ref_state, -1, quiet, _self)

    def get_alnstr(selection, state=-1, quiet=1, _self=cmd):
        with _self.lockcm:
            return _cmd.get_seq_align_str(_self._COb, str(selection),
                    int(state)-1, 0, int(quiet))

    def get_pqrstr(selection, state=-1, ref='', ref_state=-1, quiet=1, _self=cmd):
        return get_str('pqr', selection, state, ref, ref_state, -1, quiet, _self)

    def get_maestr(selection, state=-1, ref='', ref_state=-1, quiet=1, _self=cmd):
        return get_str('mae', selection, state, ref, ref_state, -1, quiet, _self)

    def get_ccp4str(name, state=1, quiet=1, _self=cmd):
        with _self.lockcm:
            return _cmd.get_ccp4str(_self._COb, str(name),
                    int(state) - 1, int(quiet))

    def get_psestr(selection, partial, quiet, _self):
        if '(' in selection: # ignore selections
            selection = ''
        return cPickle.dumps(_self.get_session(selection, partial, quiet), 1)

    def _get_mtl_obj(format, _self):
        # TODO mtl not implemented, always returns empty string
        i = {'mtl': 0, 'obj': 1}.get(format)
        return _self.get_mtl_obj()[i]

    from . import querying

    savefunctions = {
        'cif': get_str, # mmCIF
        'xyz': get_str,
        'pdb': get_str,
        'pqr': get_str,
        'sdf': get_str,
        'mol2': get_str,
        'mae': get_str,
        'mol': get_str,

        'pse': get_psestr,
        'psw': get_psestr,

        'fasta': get_fastastr,
        'aln': get_alnstr,
        'ccp4': get_ccp4str,

        'png': png,

        # no arguments (some have a "version" argument)
        'dae': querying.get_collada,
        'wrl': querying.get_vrml,
        'pov': querying.get_povray,
        'idtf': querying.get_idtf,
        'mtl': _get_mtl_obj, # TODO not implemented
        'obj': _get_mtl_obj,
    }
