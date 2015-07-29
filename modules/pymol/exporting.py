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

if __name__=='pymol.exporting':
    import os
    import thread
    import selector
    import string
    import re
    import copy
    import cPickle
    
    import pymol
    import cmd
    from cmd import _cmd,lock,unlock,Shortcut,QuietException
    from chempy import io
    from chempy.sdf import SDF,SDFRec
    from cmd import _feedback,fb_module,fb_mask, \
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
            scene_list = string.split(scenes)
            cache_mode = int(_self.get('cache_mode'))
            _self.set('cache_mode',2)
            if not len(scene_list):
                scene_list = _self.get_scene_list()
            for scene in scene_list:
                scene = string.strip(scene)
                if not quiet:
                    print " cache: optimizing scene '%s'."%scene
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
                        print " cache: no scenes defined -- optimizing current display."
                    cmd.rebuild() 
                    cmd.refresh()
            usage = _self._cache_purge(-1,_self=_self)
            if cache_mode:
                _self.set('cache_mode',cache_mode)
            else:
                _self.set('cache_mode',2) # hmm... could use 1 here instead.
            _self.set('cache_max',cache_max) # restore previous limits
            if not quiet:            
                print " cache: optimization complete (~%0.1f MB)."%(usage*4/1000000.0)
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
            if seq.has_key(obj):
                cur_seq = map(lambda x:_resn_to_aa.get(x,'?'),seq[obj])
                result.append(">%s"%obj)
                cur_seq = string.join(cur_seq,'')
                while len(cur_seq):
                    if len(cur_seq)>=70:
                        result.append(cur_seq[0:70])
                        cur_seq=cur_seq[70:]
                    else:
                        result.append(cur_seq)
                        break
        result = string.join(result,'\n')
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
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.get_pdb(_self._COb,str(selection),int(state)-1,0,
                             str(ref),int(ref_state),int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException         
        return r

    def _get_dump_str(obj):
        if is_list(obj):
            list = map(_get_dump_str,obj)
            result = "[ "+string.join(list,",\n")+" ] "
        elif is_dict(obj):
            list = []
            for key in obj.keys():
                list.append( _get_dump_str(key)+" : "+_get_dump_str(obj[key]) )
            result = "{ "+string.join(list,",\n")+" } "
        elif is_tuple(obj):
            list = map(_get_dump_str,obj)
            result = "( "+string.join(list,",\n")+" ) "
        else:
            result = str(obj)
        return result
    
    def _session_convert_legacy(session, version, _self=cmd):
        '''
        Convert session contents to be compatible with previous PyMOL versions
        '''
        if version >= _self.get_version()[1]:
            return

        print " Applying pse_export_version=%.3f compatibility" % (version)

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
                    if is_error(apply(a,(session,),{'_self':_self})):
                        r = DEFAULT_ERROR
                except:
                    traceback.print_exc()
                    print "Error: An error occurred when trying to generate session."
                    print "Error: The resulting session file may be incomplete."
        if is_ok(r):
            pse_export_version = round(_self.get_setting_float('pse_export_version'), 4)
            if pse_export_version > 0.0:
                try:
                    _session_convert_legacy(session, pse_export_version, _self)
                except Exception as e:
                    print ' Warning: failed to backport session:', e

            if(compress<0):
                compress = _self.get_setting_boolean('session_compression')
            if(compress):
                import zlib
                session = zlib.compress(cPickle.dumps(session, 1))
            return session
        elif _self._raising(r,_self):
            raise QuietException                  
        return r
        
    def _unit2px(value, dpi):
        '''API only. Returns pixel units given a string representation in other units'''
        if isinstance(value, str):
            m = re.search(r'[a-z].*', value, re.I)
            if m:
                if dpi <= 0:
                    raise pymol.CmdException('need dpi if units are given')
                value, unit = value[:m.start()], m.group(0).lower()
                upi = {'in': 1.0, 'mm': 25.4, 'cm': 2.54}
                if unit not in upi:
                    raise pymol.CmdException('unknown unit, supported units are: ' +
                            ', '.join(upi))
                value = float(value) * dpi / upi[unit] + 0.5
        return int(value)

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

    "multisave" is an unsupported command.
    
    '''
        r = DEFAULT_ERROR
        filename = _self.exp_path(filename)        
        lc_filename=string.lower(filename)
        if format=='':
            # refactor following if/elif cascade 
            # with a dictionary lookup
            if re.search("\.pdb$|\.ent$",lc_filename):
                format = 'pdb'
            elif re.search("\.pmo$",lc_filename):
                format = 'pmo'
        if format == 'pdb':
            ftype = loadable.pdb
        elif format == 'pmo':
            ftype = loadable.pmo
        else:
            ftype = loadable.pdb # default
        try:
            _self.lock(_self)
            r = _cmd.multisave(_self._COb, str(filename), str(pattern),
                               int(state)-1, int(append), int(ftype),int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def assign_atom_types( selection, format = "mol2", state=1, quiet=1, _self=cmd):
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            # format : mol2/sybyl = 1, macromodel/mmd = 2, global setting atom_type_format = 0
            r = _cmd.assign_atom_types(_self._COb, selection, int(1), int(state-1), quiet)
        finally:
            _self.unlock(r,_self)
        return r

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
        import gzip

        quiet = int(quiet)
        do_gzip = False

        # preprocess selection
        selection = selector.process(selection)
        #   
        r = DEFAULT_ERROR

        if format=='':
            ext_list = filename.lower().rsplit('.', 2)
            if ext_list[-1] == 'gz':
                do_gzip = True
                ext = ext_list[-2]
            else:
                ext = ext_list[-1]

            if ext in ['cif', 'pqr', 'mol', 'sdf', 'pkl', 'xyz', 'pov',
                    'png', 'aln', 'fasta', 'obj', 'mtl', 'wrl', 'dae', 'idtf',
                    'mol2']:
                format = ext
            elif ext in ["pdb", "ent"]:
                format = 'pdb'
            elif ext in ["mmod", "mmd", "out", "dat"]:
                format = 'mmod'
            elif ext in ["pse", "psw"]:
                format = 'pse'
            elif ext in ["pze", "pzw"]:
                do_gzip = True
                format = 'pse'
            else:
                if not quiet:
                    print " Save-Warning: Unrecognized file type -- defaulting to PDB format."
                format='pdb'

        filename = _self.exp_path(filename)

        func_type1 = {
            'cif': get_cifstr, # mmCIF
            'xyz': get_xyzstr,
            'fasta': get_fastastr,
            'aln': get_alnstr,
        }

        func_type2 = {
            'pdb': get_pdbstr,
            'pqr': get_pqrstr,
            'sdf': get_sdfstr,
            'mol2': get_mol2str,
        }

        func_type3 = {
            'dae': _self.get_collada,
            'wrl': _self.get_vrml,
            'mtl': lambda: _self.get_mtl_obj()[0],
            'obj': lambda: _self.get_mtl_obj()[1],
            'pov': lambda: ''.join(_self.get_povray()),
            'idtf': lambda: ''.join(_self.get_idtf()),
        }

        func_type4 = {
            'mmod': io.mmd.toFile,
            'pkl': io.pkl.toFile, # binary pickle
            'pkla': lambda model, filename: io.pkl.toFile(model, filename, bin=0), # ascii pickle
            'mol': io.mol.toFile,
        }

        contents = None

        if format in func_type1:
            contents = func_type1[format](selection, state, quiet, _self=_self)
        elif format in func_type2:
            contents = func_type2[format](selection, state, ref, ref_state, quiet, _self=_self)
        elif format in func_type3:
            contents = func_type3[format]()
        elif format in func_type4:
            func_type4[format](_self.get_model(selection, state, ref, ref_state), filename)
            r = DEFAULT_SUCCESS
        elif format=='pse': # PyMOL session
            filename = filename.replace("\\","/") # always use unix-like path separators	
            _self.set("session_file",filename,quiet=1)
            if '(' in selection: # ignore selections
                selection = ''
            if not quiet:
                print " Save: Please wait -- writing session file..."
            contents = cPickle.dumps(_self.get_session(selection, partial, quiet), 1)
        elif format=='png':
            return _self.png(filename, quiet=quiet)

        if isinstance(contents, basestring):
            with (gzip.open if do_gzip else open)(filename, 'wb') as handle:
                handle.write(contents)
            r = DEFAULT_SUCCESS
            
        if _self._raising(r,_self): raise QuietException

        if not quiet:
            if r == DEFAULT_SUCCESS:
                print ' Save: wrote "' + filename + '".'
            else:
                print ' Save-Error: no file written'

        return r

    # mmCIF export

    re_cifsimpledatavalue_match = re.compile(r'[^_#\$\'"\[\];]\S*$').match
    re_cifendsinglequote_search = re.compile(r"'\s").search
    re_cifenddoublequote_search = re.compile(r'"\s').search

    def cifisreserved(s):
        '''return true if s is a reserved cif keyword'''
        return s in ('loop_', 'stop_', 'global_') or s[:5] in ('data_', 'save_')

    def cifrepr(s):
        '''returns s, if s is a simple data value, or some quoted version of s'''
        if not s:
            return '.'
        if not cifisreserved(s) and re_cifsimpledatavalue_match(s) is not None:
            return s
        if '\n' not in s:
            if re_cifendsinglequote_search(s) is None:
                return "'" + s + "'"
            if re_cifenddoublequote_search(s) is None:
                return '"' + s + '"'
        if '\n;' in s:
            print ' Warning: CIF data value contains <newline><semicolon>'
            s = s.replace('\n;', '\n ;')
        return '\n;' + s + '\n;'

    def get_cifstr(selection="all", state=-1, quiet=1, _self=cmd):
        '''
DESCRIPTION

    API-only function which returns a mmCIF string.

SEE ALSO

    get_pdbstr
        '''
        tmp = _self.get_unused_name('_sele')
        _self.select(tmp, selection, 0)

        buf = ['# generated by PyMOL %s' % (_self.get_version()[0],)]

        def callback(type_, ID, elem, name, alt, resn, segi, chain, resi,
                x, y, z, q, b, formal_charge, state, entity_id):
            resv, ins = (resi[:-1], resi[-1]) if resi[-1].isalpha() else (resi, '')
            buf.append('%-6s %-3d %s %-3s '
                    '%s %-3s %s %s '
                    '%-2s %s %6.3f %6.3f %6.3f '
                    '%4.2f %6.2f %d %s %d\n' % (
                type_, ID, cifrepr(elem), cifrepr(name),
                cifrepr(alt), cifrepr(resn), cifrepr(segi),
                cifrepr(entity_id),
                resv, cifrepr(ins), x, y, z,
                q, b, formal_charge, cifrepr(chain), state,
            ))

        for model in _self.get_object_list('?' + tmp):
            buf.append('''
data_%s
_entry.id %s
''' % (model, cifrepr(model)))

            symmetry = _self.get_symmetry(model)
            if symmetry:
                buf.append('''#
_cell.entry_id %s
_cell.length_a %f
_cell.length_b %f
_cell.length_c %f
_cell.angle_alpha %f
_cell.angle_beta %f
_cell.angle_gamma %f
_symmetry.entry_id %s
_symmetry.space_group_name_H-M %s
''' % (
                    cifrepr(model),
                    symmetry[0],
                    symmetry[1],
                    symmetry[2],
                    symmetry[3],
                    symmetry[4],
                    symmetry[5],
                    cifrepr(model),
                    cifrepr(symmetry[6]),
                ))

            buf.append('''#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
''')

            loop_i = len(buf) - 1

            _self.iterate_state(state, '?%s & ?%s' % (model, tmp),
                'callback(type, ID, elem, name, alt, resn, segi, chain, resi, '
                'x, y, z, q, b, formal_charge, state, "")',
                space={'callback': callback})

            # no loop header for zero rows
            if loop_i == len(buf) - 1:
                buf[loop_i] = '#'

        _self.delete(tmp)

        return ''.join(buf)

    def get_xyzstr(selection, state=-1, quiet=1, _self=cmd):
        state = int(state)
        buf = []

        for i, (osele, ostate) in enumerate(
                pymol.selecting.objsele_state_iter(selection, state)):
            n_atoms_i = len(buf)
            buf.append('') # natoms (deferred)
            buf.append('') # comment
            n = _self.iterate_state(ostate, osele,
                    r'_buf.append("%s %f %f %f\n" % (elem, x, y, z))',
                    space={'_buf': buf})
            buf[n_atoms_i] = str(n)

        if not quiet:
            print " Save-XYZ: %d object-state(s) in selection." % (i + 1)

        return ''.join(buf)

    def get_sdfstr(selection, state=-1, ref='', ref_state=-1, quiet=1, _self=cmd):
        state = int(state)
        buf = []

        for i, (osele, ostate) in enumerate(
                pymol.selecting.objsele_state_iter(selection, state)):
            rec = SDFRec(io.mol.toList(_self.get_model(osele, ostate, ref, ref_state)))
            buf.extend(rec.toList())
            buf.append('$$$$\n')

        if not quiet:
            print " Save-SDF: %d object-state(s) in selection." % (i + 1)

        return ''.join(buf)

    def get_mol2str(selection, state=-1, ref='', ref_state=-1, quiet=1, _self=cmd):
        state = int(state)
        buf = []

        for i, (osele, ostate) in enumerate(
                pymol.selecting.objsele_state_iter(selection, state)):
            assign_atom_types(osele, "mol2", ostate, 1, _self)
            buf.extend(io.mol2.toList(_self.get_model(osele,
                ostate, ref, ref_state), selection=osele, state=ostate))

        if not quiet:
            print " Save-MOL2: %d object-state(s) in selection." % (i + 1)

        return ''.join(buf)

    def get_alnstr(selection, state=-1, quiet=1, _self=cmd):
        with _self.lockcm:
            return _cmd.get_seq_align_str(_self._COb, str(selection),
                    int(state)-1, 0, int(quiet))

    def get_pqrstr(selection, state=-1, ref='', ref_state=-1, quiet=1, _self=cmd):
        with _self.lockcm:
            return _cmd.get_pdb(_self._COb, str(selection), int(state)-1, 1,
                    str(ref), int(ref_state)-1, int(quiet))
