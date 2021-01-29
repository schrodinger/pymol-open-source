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

# must match layer1/Setting.h
cSetting_tuple = -1
cSetting_blank = 0
cSetting_boolean = 1
cSetting_int = 2
cSetting_float = 3
cSetting_float3 = 4
cSetting_color = 5
cSetting_string = 6

if True:

    import traceback
    from . import selector
    from .shortcut import Shortcut
    cmd = __import__("sys").modules["pymol.cmd"]
    from .cmd import _cmd,lock,lock_attempt,unlock,QuietException, \
          is_string, \
          _feedback,fb_module,fb_mask, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error

    # name -> index mapping
    index_dict = _cmd.get_setting_indices()

    # index -> name mapping
    name_dict = dict((v,k) for (k,v) in index_dict.items())

    name_list = list(index_dict.keys())
    setting_sc = Shortcut(name_list)

    # legacy
    index_dict['ray_shadows'] =     index_dict['ray_shadow']

    boolean_dict = {
        "true" : 1,
        "false": 0,
        "on"   : 1,
        "off"  : 0,
        "1"    : 1,
        "0"    : 0,
        "1.0"  : 1,
        "0.0"  : 0,
        }

    boolean_sc = Shortcut(boolean_dict.keys())

    def _get_index(name):
        '''Get setting index for given name. `name` may be abbreviated.
        Raises QuietException for unknown names or ambiguous abbreviations.'''
        if isinstance(name, int) or name.isdigit():
            return int(name)
        if name not in index_dict:
            name = setting_sc.auto_err(name, 'Setting')
        return index_dict[name]

    def _get_name(index):
        # legacy, in case someone used that in a script
        return name_dict.get(index, "")

    def get_index_list():
        # legacy, in case someone used that in a script (e.g. grepset)
        return list(name_dict.keys())

    def get_name_list():
        return name_list

    def _validate_value(type, value):
        if type == cSetting_boolean:  # (also support non-zero float for truth)
            try: # number, non-zero, then interpret as TRUE
                return 1 if float(value) else 0
            except:
                pass
            return boolean_dict[boolean_sc.auto_err(str(value), "boolean")]
        if type in (cSetting_int, cSetting_float):
            if is_string(value) and boolean_sc.has_key(value):
                value = boolean_dict[boolean_sc.auto_err(str(value), "boolean")]
            if type == cSetting_int:
                return int(value)
            return float(value)
        if type == cSetting_float3:  # some legacy handling req.
            if not is_string(value):
                v = value
            elif ',' in value:
                v = cmd.safe_eval(value)
            else:
                v = value.split()
            return (float(v[0]), float(v[1]), float(v[2]))
        if type == cSetting_color:
            return str(value)
        if type == cSetting_string:
            v = str(value)
            # strip outermost quotes (cheesy approach)
            if len(v) > 1 and v[0] == v[-1] and v[0] in ('"', "'"):
                v = v[1:-1]
            return v
        raise Exception

    ###### API functions

    def set_bond(name, value, selection1, selection2=None,
                 state=0, updates=1, log=0, quiet=1, _self=cmd):
        ''' 
DESCRIPTION

    "set_bond" changes per-bond settings for all bonds which exist
    between two selections of atoms.

USAGE

    set_bond name, value, selection1 [, selection2 ]

ARGUMENTS

    name = string: name of the setting

    value = string: new value to use

    selection1 = string: first set of atoms

    selection2 = string: seconds set of atoms {default: (selection1)}

EXAMPLE

    set_bond stick_transparency, 0.7, */n+c+ca+o


NOTES

    The following per-bond settings are currently implemented.  Others
    may seem to be recognized but will currently have no effect when
    set at the per-bond level.
    
    * valence
    * line_width
    * line_color
    * stick_radius
    * stick_color
    * stick_transparency

    Note that if you attempt to use the "set" command with a per-bond
    setting over a selection of atoms, the setting change will appear
    to take, but no change will be observed.
    
PYMOL API

    cmd.set_bond ( string name, string value,
                   string selection1,
                   string selection2,
                   int state, int updates, log=0, quiet=1)

       '''
        selection1 = selector.process(selection1)
        selection2 = selector.process(selection2) if selection2 else selection1
        index = _get_index(str(name))
        type = _cmd.get_setting_type(index)
        v = (type, _validate_value(type, value))
        if log:
            name = name_dict.get(index, name)
            _self.log(
                '',
                f"cmd.set_bond({name!r},{value!r},{selection1!r},{selection2!r},{state})\n"
            )
        with _self.lockcm:
            r = _cmd.set_bond(_self._COb, index, v, selection1, selection2,
                              int(state) - 1, int(quiet), int(updates))
        return r


    def set(name, value=1, selection='', state=0, updates=1, log=0,
            quiet=1,_self=cmd):

        '''
DESCRIPTION

    "set" changes global, object, object-state, or per-atom settings.

USAGE

    set name [,value [,selection [,state ]]]

ARGUMENTS

    name = string: setting name

    value = string: a setting value {default: 1}

    selection = string: name-pattern or selection-expression
    {default:'' (global)}

    state = a state number {default: 0 (per-object setting)}

EXAMPLES

    set orthoscopic

    set line_width, 3

    set surface_color, white, 1hpv
    
    set sphere_scale, 0.5, elem C

NOTES

    The default behavior (with a blank selection) is global.  If the
    selection is "all", then the setting entry in each individual
    object will be changed.  Likewise, for a given object, if state is
    zero, then the object setting will be modified.  Otherwise, the
    setting for the indicated state within the object will be
    modified.

    If a selection is provided as opposed to an object name, then the
    atomic setting entries are modified.

    The following per-atom settings are currently implemented.  Others
    may seem to be recognized but will have no effect when set on a
    per-atom basis.
    
    * sphere_color
    * surface_color
    * mesh_color
    * label_color
    * dot_color
    * cartoon_color
    * ribbon_color
    * transparency (for surfaces)
    * sphere_transparency
    
    Note that if you attempt to use the "set" command with a per-bond
    setting over a selection of atoms, the setting change will appear
    to take, but no change will be observed.  Please use the
    "set_bond" command for per-bond settings.
    

PYMOL API

    cmd.set(string name, string value, string selection, int state,
            int updates, int quiet)

SEE ALSO

    get, set_bond
    
'''
        selection = selector.process(selection)
        index = _get_index(name)
        type = _cmd.get_setting_type(index)
        v = (type, _validate_value(type, value))
        if log:
            name = name_dict.get(index, name)
            _self.log('',
                      f"cmd.set({name!r},{value!r},{selection!r},{state})\n")
        with _self.lockcm:
            r = _cmd.set(_self._COb, index, v, selection,
                         int(state) - 1, int(quiet), int(updates))
        return r

    def unset(name, selection='', state=0, updates=1, log=0, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "unset" clears a setting and restores its default value.

    WARNING: The behavior for global settings changed in PyMOL 2.5.
    Previously, "unset settingname" would set the global value of
    "settingname" to zero/off instead of the default value.
    To set a setting to zero, do "set settingname, 0".

USAGE

    unset name [,selection [,state ]]

EXAMPLE

    unset orthoscopic

    unset surface_color, 1hpv

    unset sphere_scale, elem C
    
NOTES

    If selection is not provided, unset changes the named global
    setting to its default value.

    If a selection is provided, then "unset" undefines per-object,
    per-state, or per-atom settings.

PYMOL API

    cmd.unset(string name, string selection, int state, int updates,
                int log)

SEE ALSO

    unset_deep, set, set_bond
    
        '''
        selection = selector.process(selection)
        index = _get_index(str(name))
        if log:
            name = name_dict.get(index, name)
            _self.log('', f"cmd.unset({name!r},{selection!r},{state})\n")
        with _self.lockcm:
            r = _cmd.unset(_self._COb, index, selection,
                           int(state) - 1, int(quiet), int(updates))
        return r

    def unset_bond(name,selection1,selection2=None,state=0,updates=1,log=0,quiet=1,_self=cmd):
        '''
DESCRIPTION

    "unset_bond" removes a per-bond setting for a given set of bonds.
    
USAGE

    unset name [,selection [, selection [,state ]]]

        '''
        selection1 = selector.process(selection1)
        selection2 = selector.process(selection2) if selection2 else selection1
        index = _get_index(str(name))
        if log:
            name = name_dict.get(index, name)
            _self.log(
                '',
                f"cmd.unset_bond({name!r},{selection1!r},{selection2!r},{state})\n"
            )
        with _self.lockcm:
            r = _cmd.unset_bond(_self._COb,int(index),selection1,selection2,
                                   int(state)-1,int(quiet),
                                   int(updates))
        return r

    def get_setting(name,object='',state=0,_self=cmd): # INTERNAL
        return get_setting_tuple_new(name, object, state, _self)[1]

    def get(name, selection='', state=0, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "get" prints out the current value of a setting.

USAGE

    get name [, selection [, state ]]
    
EXAMPLE

    get line_width

ARGUMENTS

    name = string: setting name

    selection = string: object name (selections not yet supported)

    state = integer: state number
    
NOTES

    "get" currently only works with global, per-object, and per-state
    settings.  Atom level settings get be queried with "iterate" (e.g.
    iterate all, print s.line_width)
    
PYMOL API

    cmd.get(string name, string object, int state, int quiet)

SEE ALSO

    set, set_bond, get_bond

    '''

        state = int(state)
        i = _get_index(name)
        r = get_setting_text(i, str(selection), state, _self)
        if is_ok(r) and (r is not None):
            if not int(quiet):
                name = name_dict.get(i, name)
                r_str = str(r)
                if len(r_str) > 200:
                    r_str = r_str[:185] + '... (truncated)'
                if(selection==''):
                    print(" get: %s = %s"%(name,r_str))
                elif state<=0:
                    print(" get: %s = %s in object %s"%(name,r_str,selection))
                else:
                    print(" get: %s = %s in object %s state %d"%(name,r_str,selection,state))
        return r

    def get_setting_tuple_new(name,object='',state=0,_self=cmd): # INTERNAL
        i = _get_index(name)
        with _self.lockcm:
            return _cmd.get_setting_of_type(_self._COb, i, str(object), int(state) - 1, cSetting_tuple)

    def get_setting_tuple(name,object='',state=0,_self=cmd): # INTERNAL
        r = get_setting_tuple_new(name, object, state, _self)
        if r[0] != cSetting_float3:
            # legacy API
            r = (r[0], (r[1],))
        return r

    def get_setting_boolean(name,object='',state=0,_self=cmd): # INTERNAL
        i = _get_index(name)
        with _self.lockcm:
            return _cmd.get_setting_of_type(_self._COb, i, str(object), int(state) - 1, cSetting_boolean)

    def get_setting_int(name,object='',state=0,_self=cmd): # INTERNAL
        i = _get_index(name)
        with _self.lockcm:
            return _cmd.get_setting_of_type(_self._COb, i, str(object), int(state) - 1, cSetting_int)

    def get_setting_float(name,object='',state=0,_self=cmd): # INTERNAL
        i = _get_index(name)
        with _self.lockcm:
            return _cmd.get_setting_of_type(_self._COb, i, str(object), int(state) - 1, cSetting_float)

    def get_setting_text(name,object='',state=0,_self=cmd):  # INTERNAL
        i = _get_index(name)
        with _self.lockcm:
            return _cmd.get_setting_of_type(_self._COb, i, str(object), int(state) - 1, cSetting_string)

    def get_setting_updates(object='', state=0, _self=cmd): # INTERNAL
        r = []
        if lock_attempt(_self):
            try:
                r = _cmd.get_setting_updates(_self._COb, object, state-1)
            finally:
                _self.unlock(r,_self)
        return r

    def get_bond(name, selection1, selection2=None,
                 state=0, updates=1, quiet=1, _self=cmd):
        ''' 
DESCRIPTION

    "get_bond" gets per-bond settings for all bonds which exist
    between two selections of atoms.

USAGE

    get_bond name, selection1 [, selection2 ]

ARGUMENTS

    name = string: name of the setting

    selection1 = string: first set of atoms

    selection2 = string: seconds set of atoms {default: (selection1)}

EXAMPLE

    get_bond stick_transparency, */n+c+ca+o


NOTES

    The following per-bond settings are currently implemented.  Others
    may seem to be recognized but will currently have no effect when
    set at the per-bond level.
    
    * valence
    * line_width
    * line_color
    * stick_radius
    * stick_color
    * stick_transparency

PYMOL API

    cmd.get_bond ( string name,
                   string selection1,
                   string selection2,
                   int state, int updates, quiet=1)

       '''
        state, quiet = int(state), int(quiet)
        selection1 = selector.process(selection1)
        selection2 = selector.process(selection2) if selection2 else selection1

        index = _get_index(str(name))
        with _self.lockcm:
            r = _cmd.get_bond(_self._COb, index, selection1, selection2,
                              state - 1, quiet, int(updates))

        if not quiet:
            name = name_dict.get(index, name)
            suffix = ' state %d' % state if state > 0 else ''
            for model, vlist in r:
                print(' %s = %s for object %s' % (name, _self.get(name, model), model))
                for idx1, idx2, value in vlist:
                    if value is None:
                        continue
                    print(' %s = %s between (%s`%d)-(%s`%d%s)' % (name,
                            value, model, idx1, model, idx2, suffix))
        return r

    def unset_deep(settings='', object='*', updates=1, quiet=1, _self=cmd):
        '''
DESCRIPTION

    Unset all object, object-state, atom, and bond level settings.

    Note: Does currently NOT unset atom-state level settings. Check for
    atom-state level settings with:
    PyMOL> iterate_state 1, *, print(list(s))
    Unset e.g. atom-state level "label_screen_point" (index 728) with:
    PyMOL> alter_state 1, *, del s[728]

ARGUMENTS

    settings = str: space separated list of setting names or empty string
    for all settings {default: }

    object = str: name of one object or * for all objects {default: *}
        '''
        quiet = int(quiet)
        kwargs = {'quiet': quiet, 'updates': 0, '_self': _self}

        if not settings:
            settings = iter(name_dict) # index iterator
        elif _self.is_string(settings):
            settings = settings.split()

        if object in ['all', '*']:
            object = '*'
            selection = '(*)'
        else:
            selection = None
            try:
                if _self.get_type(object) in (
                        'object:group', 'object:molecule'):
                    selection = '(' + object + ')'
            except:
                pass

        # 0 (object-level) and 1-N (object-state-level)
        states = range(_self.count_states(object) + 1)

        for setting in settings:
            try:
                for state in states:
                    unset(setting, object, state=state, **kwargs)
                if selection:
                    unset(setting, selection, **kwargs)
            except:
                if not quiet:
                    print(' Setting: %s unset failed' % setting)

        if int(updates):
            _self.rebuild(object)
