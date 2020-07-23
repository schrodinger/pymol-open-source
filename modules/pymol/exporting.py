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

from . import colorprinting

if True:
    import os
    import sys
    import _thread as thread
    from . import selector
    import re
    import copy

    import pymol
    cmd = sys.modules["pymol.cmd"]
    from .cmd import _cmd,Shortcut,QuietException
    from chempy import io
    from chempy.pkl import cPickle
    from .cmd import _feedback,fb_module,fb_mask, \
                     DEFAULT_ERROR, DEFAULT_SUCCESS, is_error, \
                     is_list, is_dict, is_tuple, loadable

    def copy_image(quiet=1, *, _self=cmd): # incentive feature / proprietary
        if _self.is_gui_thread():
            _self._copy_image(_self, int(quiet))
        else:
            _self.do('cmd._copy_image(quiet=%d)' % int(quiet))

    cache_action_dict = {
        'enable'      : 0,
        'disable'     : 1,
        'read_only'   : 2,
        'clear'       : 3,
        'optimize'    : 4,
    }

    cache_action_sc = Shortcut(cache_action_dict.keys())

    def cache(action='optimize', scenes='',state=-1, quiet=1, *, _self=cmd):
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

        action = cache_action_dict[cache_action_sc.auto_err(str(action),'action')]
        quiet = int(quiet)
        if action == 0: # enable
            _self.set('cache_mode', 2, quiet=quiet)
        elif action == 1:  # disable
            _self.set('cache_mode', 0, quiet=quiet)
        elif action == 2:  # read_only
            _self.set('cache_mode', 1, quiet=quiet)
        elif action == 3: # clear
            _self._cache_clear()
        elif action == 4: # optimize
            _self._cache_mark()
            cur_scene = _self.get('scene_current_name')
            cache_max = _self.get_setting_int('cache_max')
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
                _self.scene(scene,animate=0)
                _self.rebuild()
                _self.refresh()
            if len(cur_scene):
                _self.scene(cur_scene,animate=0)
            else:
                scene_list = _self.get_scene_list()
                if len(scene_list):
                    _self.scene(scene_list[0],animate=0)
                else:
                    if not quiet:
                        print(" cache: no scenes defined -- optimizing current display.")
                    _self.rebuild()
                    _self.refresh()
            usage = _self._cache_purge(-1)
            if cache_mode:
                _self.set('cache_mode',cache_mode)
            else:
                _self.set('cache_mode',2) # hmm... could use 1 here instead.
            _self.set('cache_max',cache_max) # restore previous limits
            if not quiet:
                print(" cache: optimization complete (~%0.1f MB)."%(usage*4/1000000.0))
        else:
            raise ValueError('action')

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
            # RNA
            'A'   : 'A',
            'U'   : 'U',
            'G'   : 'G',
            'C'   : 'C',
            # DNA
            'DA'  : 'A',
            'DT'  : 'T',
            'DG'  : 'G',
            'DC'  : 'C',
            }

    def get_fastastr(selection="all", state=-1, quiet=1, key='', *, _self=cmd):
        '''
DESCRIPTION

    API only. Get protein and nucleic acid sequences in fasta format.

    Used for saving:
    PyMOL> save foo.fasta

    New in PyMOL 2.2:
    - chain specific keys (key argument)
    - nucleic acid support

ARGUMENTS

    selection = str: atom selection (reduced to "guide & alt +A") {default: all}

    state = int: (only used if state > 0)

    quiet = 0/1: UNUSED

    key = str: python expression {default: model + "_" + chain}
    Use key=model to get the old (non-chain specific) behavior
        '''
        import textwrap
        import collections

        seq = collections.OrderedDict()
        lines = []

        selection = "(" + selection + ") & guide & alt +A"
        if int(state) > 0:
            selection += ' & state {}'.format(state)

        if not key:
            key = 'model + "_" + chain'
            # for discrete objects: state specific keys
            key += ' + (":{}".format(state) if state else "")'

        _self.iterate(selection,
            "seq.setdefault((" + key + "),[]).append(resn)",
            space={'seq': seq, 'str': str})

        for key, resn_list in seq.items():
            cur_seq = ''.join(_resn_to_aa.get(x, '?') for x in resn_list)
            lines.append(">{}".format(key))
            lines.extend(textwrap.wrap(cur_seq, 70))

        if lines:
            lines.append('')  # final newline
        return '\n'.join(lines)

    def get_pdbstr(selection="all", state=-1, ref='', ref_state=-1, quiet=1, *, _self=cmd):
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
        return _self.get_str('pdb', selection, state, ref, ref_state, -1, quiet)

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

    def get_session(names='', partial=0, quiet=1, compress=-1, cache=-1,
                    binary=-1, version=-1,
                    *, _self=cmd):
        '''
        :param names: Names of objects to export, or the empty string to export all objects.
        :param partial: If true, do not store selections, settings, view, movie.
        :param compress: DEPRECATED Pickle the dictionary to a string and
        compress that with zlib {default: session_compression}
        :param cache: ?
        :param binary: Use efficient binary format {default: pse_binary_dump}
        :param version: {default: pse_export_version}
        '''
        session = {}
        cache = int(cache)
        compress = int(compress)
        partial = int(partial)
        quiet = int(quiet)

        if version < 0:
            version = _self.get_setting_float('pse_export_version')

        pse_export_version = round(float(version), 4)
        legacyscenes = (0 < pse_export_version < 1.76) and _self.get_scene_list()

        if True:
            legacypickle = (0 < pse_export_version < 1.9)
            if legacypickle:
                print(' Using Python 2 compatible legacy pickler')
            cPickle.configure_legacy_dump(legacypickle)

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
            cache_opt = _self.get_setting_int('session_cache_optimize')
            if cache_opt == 1 or (cache_opt != 0
                                  and _self.get_setting_int('cache_mode')):
                _self.cache('optimize')

        with _self.lockcm:
            _cmd.get_session(_self._COb, session, str(names), int(partial),
                             int(quiet), binary, pse_export_version)

        if True:
                try:
                    session['session'] = copy.deepcopy(_self._pymol.session)
                    if cache and hasattr(_self._pymol,'_cache'):
                        session['cache'] = _self._pymol._cache
                except:
                    colorprinting.print_exc()

        for a in _self._pymol._session_save_tasks:
            assert a is not None
            error = False
            try:
                # TODO _session_restore_tasks use `not <return-value>`
                # instead of `is_error(<return-value>)`
                error = is_error(a(session, _self=_self))
            except:
                if not quiet:
                    colorprinting.print_exc([__file__])
                error = True
            if error:
                msg = f'Warning: session-save-task "{a.__name__}" failed.'
                if not quiet:
                    msg += '\n         The resulting session file may be incomplete.'
                colorprinting.warning(msg)

        if legacyscenes:
            del session['moviescenes']
            session['scene_dict'] = _self.pymol._scene_dict
            session['scene_order'] = legacyscenes
            pymol.viewing._legacy_scene('*', 'clear', _self=_self)
            del _self.pymol._scene_dict

        if True:
            if pse_export_version > 0.0:
                try:
                    _session_convert_legacy(session, pse_export_version, _self)
                except Exception as e:
                    print(' Warning: failed to backport session:', e)

            if(compress<0):
                compress = _self.get_setting_boolean('session_compression')
            if(compress):
                colorprinting.warning(
                    ' Warning: `session_compression` is deprecated. Save '
                    'to ".pze" or ".pse.gz" files instead, and/or use '
                    '`pse_binary_dump`')
                import zlib
                session = zlib.compress(cPickle.dumps(session, 1))
            return session

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
            quiet=1, prior=0, format=0, *, _self=cmd):
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
                           int(format))
            if r != 1: # no prior image available -- revert to default behavior
                if prior < 0: # default is to fall back to actual rendering
                    prior = 0
        if not prior:
            dpi = float(dpi)
            if dpi < 0:
                dpi = _self.get_setting_float('image_dots_per_inch')
            width = _unit2px(width, dpi)
            height = _unit2px(height, dpi)

            if _self.is_gui_thread():
                r = _self._png(str(filename),int(width),int(height),float(dpi),
                               int(ray),int(quiet),0,int(format))
            else:
                r = _self._do("cmd._png('''%s''',%d,%d,%1.6f,%d,%d,%d,%d)"%
                              (filename,width,height,dpi,
                               ray,int(quiet),0,int(format)),_self=_self)
        if _self._raising(r): raise QuietException
        return r

    def multisave(filename, pattern="all", state=-1,
                  append=0, format='', quiet=1, *, _self=cmd):
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

        s = _self.get_str(format, pattern, state, '', -1, 1, quiet)

        if s is None:
            raise QuietException

        filename = _self.exp_path(filename)

        with open(filename, 'a' if int(append) else 'w') as handle:
            handle.write(s)

        return DEFAULT_SUCCESS

    def assign_atom_types( selection, format = "mol2", state=1, quiet=1, *, _self=cmd):
        fmt = {'mol2': 1, 'mmd': 2}[format]
        with _self.lockcm:
            # format : mol2/sybyl = 1, macromodel/mmd = 2, global setting atom_type_format = 0
            _cmd.assign_atom_types(_self._COb, selection, fmt, int(state - 1),
                                   quiet)

    def get_str(format, selection='(all)', state=-1, ref='',
             ref_state=-1, multi=-1, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    Like "get_bytes" but return a unicode string.
        '''
        assert format not in ('mmtf',), 'binary format, use get_bytes'
        b = _self.get_bytes(format, selection, state, ref, ref_state, multi, quiet)
        if b is None:
            return None
        return b.decode('utf-8')

    def get_bytes(format, selection='(all)', state=-1, ref='',
             ref_state=-1, multi=-1, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    API-only function which exports the selection to a molecular file
    format and returns it as a binary ("bytes") string.

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

    def multifilesave(filename, selection='*', state=-1, format='', ref='',
             ref_state=-1, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    For a selection that spans multiple molecular objects and/or states,
    save each object and/or state to a separate file. Takes a filename
    argument with placeholders:

    {name}  : object name
    {state} : state number
    {title} : state title
    {num}   : file number
    {}      : object name (first) or state (second)

EXAMPLES

    multifilesave /tmp/{name}.pdb
    multifilesave /tmp/{name}-{state}.cif, state=0
    multifilesave /tmp/{}-{}.cif, state=0
    multifilesave /tmp/{}-{title}.sdf, state=0
        '''
        for (fname, osele, ostate) in multifilenamegen(
                filename, selection, int(state), _self):
            r = _self.save(fname, osele, ostate, format, ref, ref_state, quiet)
        return r


    def multifilenamegen(filename, selection, state, _self=cmd):
        '''Given a filename pattern, atom selection and state argument,
        Generate object-state specific filenames and selections.
        '''
        import string

        filename = _self.exp_path(filename)

        fmt_keys = [v[1]
                for v in string.Formatter().parse(filename)
                if v[1] is not None]

        nindexed = fmt_keys.count('')
        multiobject = nindexed > 0 or 'name' in fmt_keys or 'num' in fmt_keys
        multistate = nindexed > 1 or 'state' in fmt_keys or 'title' in fmt_keys

        if not (multiobject or multistate):
            raise ValueError('need one or more of {name}, {num}, {state}, {title}')

        odata = []
        for oname in _self.get_object_list(selection):
            osele = '(%s) & ?%s' % (selection, oname)
            first = last = state

            if multistate:
                if state < 0:
                    first = last = _self.get_object_state(oname)

                if first == 0:
                    first = 1
                    last = _self.count_states('%' + oname)

            for ostate in range(first, last + 1):
                odata.append((oname, osele, ostate))

        # pad {state} and {num} with zeros
        swidth = len(str(max(v[2] for v in odata)))
        nwidth = len(str(len(odata)))
        filename = filename.replace('{state}', '{state:0%d}' % swidth)
        filename = filename.replace('{num}', '{num:0%d}' % nwidth)

        for num, (oname, osele, ostate) in enumerate(odata, 1):
            title = _self.get_title(oname, ostate) if multistate else ''
            fname = filename.format(oname, ostate,
                    name=oname, state=ostate, num=num,
                    title=title)
            yield fname, osele, ostate


    def save(filename, selection='(all)', state=-1, format='', ref='',
             ref_state=-1, quiet=1, partial=0, *, _self=cmd):
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
        from pymol.importing import filename_to_format, _eval_func
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
            func = _eval_func(func)
            kw_all = {
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
            sig = inspect.signature(func, follow_wrapped=False)
            kw = {}

            for n, param in sig.parameters.items():
                if param.kind == inspect.Parameter.VAR_KEYWORD:
                    kw = kw_all
                    break

                if param.kind == inspect.Parameter.VAR_POSITIONAL:
                    print('FIXME: savefunctions[%s]: *args' % (format))
                elif param.kind == inspect.Parameter.POSITIONAL_ONLY:
                    raise Exception('positional-only arguments not supported')
                elif n in kw_all:
                    kw[n] = kw_all[n]

            contents = func(**kw)

            if 'filename' in sig.parameters:
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

        if _self._raising(r): raise QuietException

        if not quiet:
            if r == DEFAULT_SUCCESS:
                print(' Save: wrote "' + filename + '".')
            else:
                print(' Save-Error: no file written')

        return r

    def get_cifstr(selection="all", state=-1, quiet=1, *, _self=cmd):
        '''
DESCRIPTION

    API-only function which returns a mmCIF string.

SEE ALSO

    get_pdbstr
        '''
        return _self.get_str('cif', selection, state, '', -1, -1, quiet)

    def get_xyzstr(selection, state=-1, quiet=1, *, _self=cmd):
        return _self.get_str('xyz', selection, state, '', -1, -1, quiet)

    def get_sdfstr(selection, state=-1, ref='', ref_state=-1, quiet=1, *, _self=cmd):
        return _self.get_str('sdf', selection, state, ref, ref_state, -1, quiet)

    def get_mol2str(selection, state=-1, ref='', ref_state=-1, quiet=1, *, _self=cmd):
        return _self.get_str('mol2', selection, state, ref, ref_state, -1, quiet)

    def get_alnstr(selection, state=-1, quiet=1, *, _self=cmd):
        with _self.lockcm:
            return _cmd.get_seq_align_str(_self._COb, str(selection),
                    int(state)-1, 0, int(quiet))

    def get_pqrstr(selection, state=-1, ref='', ref_state=-1, quiet=1, *, _self=cmd):
        return _self.get_str('pqr', selection, state, ref, ref_state, -1, quiet)

    def get_maestr(selection, state=-1, ref='', ref_state=-1, quiet=1, *, _self=cmd):
        return _self.get_str('mae', selection, state, ref, ref_state, -1, quiet)

    def get_ccp4str(name, state=1, quiet=1, format='ccp4', *, _self=cmd):
        ftype = getattr(loadable, format, -1)
        with _self.lockcm:
            return _cmd.get_ccp4str(_self._COb, str(name),
                    int(state) - 1, int(quiet), ftype)

    def get_psestr(selection, partial, quiet, _self):
        if '(' in selection: # ignore selections
            selection = ''
        session = _self.get_session(selection, partial, quiet)
        return cPickle.dumps(session, 1)

    def _get_mtl_obj(format, _self):
        # TODO mtl not implemented, always returns empty string
        if format == 'mtl':
            raise pymol.CmdException('.MTL export not implemented')
        i = {'mtl': 0, 'obj': 1}.get(format)
        return _self.get_mtl_obj()[i]

    savefunctions = {
        'cif': get_str, # mmCIF
        'xyz': get_str,
        'pdb': get_str,
        'pqr': get_str,
        'sdf': get_str,
        'mol2': get_str,
        'mae': get_str,
        'mol': get_str,
        'mmtf': get_bytes,

        'pse': get_psestr,
        'psw': get_psestr,

        'fasta': get_fastastr,
        'aln': get_alnstr,
        'ccp4': get_ccp4str,
        'mrc': get_ccp4str,
        'map': get_ccp4str,

        'png': png,

        # no arguments (some have a "version" argument)
        'dae': 'pymol.querying:get_collada',
        'gltf': 'pymol.querying:get_gltf',
        'wrl': 'pymol.querying:get_vrml',
        'pov': 'pymol.querying:get_povray',
        'idtf': 'pymol.querying:get_idtf',
        'mtl': _get_mtl_obj, # TODO not implemented
        'obj': _get_mtl_obj,
        'stl': 'pymol.lazyio:get_stlstr',
    }
