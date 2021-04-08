import os
import sys
cmd = sys.modules["pymol.cmd"]
from pymol import _cmd
import threading
import traceback

import _thread as thread
import urllib.request as urllib2

import re
import time
import pymol

import chempy.io

from .cmd import DEFAULT_ERROR, DEFAULT_SUCCESS, loadable, _load2str, Shortcut, \
   is_string, is_ok

# cache management:

def _cache_validate(_self=cmd):
    r = DEFAULT_SUCCESS
    with _self.lock_api_data:
        _pymol = _self._pymol
        if not hasattr(_pymol,"_cache"):
            _pymol._cache = []
        if not hasattr(_pymol,"_cache_memory"):
            _pymol._cache_memory = 0

def _cache_clear(_self=cmd):
    r = DEFAULT_SUCCESS
    with _self.lock_api_data:
        _pymol = _self._pymol
        _pymol._cache = []
        _pymol._cache_memory = 0
    return r

def _cache_mark(_self=cmd):
    r = DEFAULT_SUCCESS
    with _self.lock_api_data:
        _pymol = _self._pymol
        _cache_validate(_self)
        for entry in _self._pymol._cache:
            entry[5] = 0.0
    return r

def _cache_purge(max_size, _self=cmd):
    r = DEFAULT_SUCCESS
    with _self.lock_api_data:
        _pymol = _self._pymol
        _cache_validate(_self)
        if len(_pymol._cache):
            cur_size = sum(x[0] for x in _pymol._cache)
            if max_size>=0: # purge to reduce size
                now = time.time()
                # sort by last access time
                new_cache = [[(now-x[5])/x[4],x] for x in _pymol._cache]
                new_cache.sort()
                new_cache = [x[1] for x in new_cache]
                # remove oldest entries one by one until size requirement is met
                while (cur_size>max_size) and (len(new_cache)>1):
                    entry = new_cache.pop()
                    cur_size = cur_size - entry[0]
                _pymol._cache = new_cache
                _pymol._cache_memory = cur_size
            else: # purge to eliminate unused entries
                new_cache = []
                for entry in _pymol._cache:
                    if entry[5] == 0.0:
                        cur_size = cur_size - entry[0]
                    else:
                        new_cache.append(entry)
                _pymol._cache = new_cache
                _pymol._cache_memory = cur_size
        result = _pymol._cache_memory
    return result

def _cache_get(target, hash_size = None, _self=cmd):
    result = None
    with _self.lock_api_data:
        try:
            if hash_size is None:
                hash_size = len(target[1])
            key = target[1][0:hash_size]
            # should optimize this with a dictionary lookup, key -> index in _cache
            for entry in _self._pymol._cache:
                if entry[1][0:hash_size] == key:
                    if entry[2] == target[2]:
                        while len(entry)<6:
                            entry.append(0)
                        entry[4] = entry[4] + 1 # access count
                        entry[5] = time.time() # timestamp
                        result = entry[3]
                        break
        except:
            traceback.print_exc()
    return result

def _cache_set(new_entry, max_size, _self=cmd):
    r = DEFAULT_SUCCESS
    with _self.lock_api_data:
        _pymol = _self._pymol
        _cache_validate(_self)
        try:
            hash_size = len(new_entry[1])
            key = new_entry[1][0:hash_size]
            count = 0
            found = 0
            new_entry[4] = new_entry[4] + 1 # incr access count
            new_entry[5] = time.time() # timestamp
            for entry in _pymol._cache:
                if entry[1][0:hash_size] == key:
                    if entry[2] == new_entry[2]: # dupe (shouldn't happen)
                        entry[3] = new_entry[3]
                        found = 1
                        break
                count = count + 1
            if not found:
                _pymol._cache.append(new_entry)
                _pymol._cache_memory = _pymol._cache_memory + new_entry[0]
                if max_size > 0:
                    if _pymol._cache_memory > max_size:
                        _cache_purge(max_size, _self)
        except:
            traceback.print_exc()
    return r

# ray tracing threads

def _ray_anti_spawn(thread_info,_self=cmd):
    # WARNING: internal routine, subject to change
    # internal routine to support multithreaded raytracing
    thread_list = []
    for a in thread_info[1:]:
        t = threading.Thread(target=_cmd.ray_anti_thread,
                                    args=(_self._COb,a))
        t.setDaemon(1)
        thread_list.append(t)
    for t in thread_list:
        t.start()
    _cmd.ray_anti_thread(_self._COb,thread_info[0])
    for t in thread_list:
        t.join()

def _ray_hash_spawn(thread_info,_self=cmd):
    # WARNING: internal routine, subject to change
    # internal routine to support multithreaded raytracing
    thread_list = []
    for a in thread_info[1:]:
        if a is not None:
            t = threading.Thread(target=_cmd.ray_hash_thread,
                                 args=(_self._COb,a))
            t.setDaemon(1)
            thread_list.append(t)
    for t in thread_list:
        t.start()
    if thread_info[0] is not None:
        _cmd.ray_hash_thread(_self._COb,thread_info[0])
    for t in thread_list:
        t.join()

def _ray_spawn(thread_info,_self=cmd):
    # WARNING: internal routine, subject to change
    # internal routine to support multithreaded raytracing
    thread_list = []
    for a in thread_info[1:]:
        t = threading.Thread(target=_cmd.ray_trace_thread,
                                    args=(_self._COb,a))
        t.setDaemon(1)
        thread_list.append(t)
    for t in thread_list:
        t.start()
    _cmd.ray_trace_thread(_self._COb,thread_info[0])
    for t in thread_list:
        t.join()

def _coordset_update_thread(list_lock,thread_info,_self=cmd):
    # WARNING: internal routine, subject to change
    while 1:
        list_lock.acquire()
        if not len(thread_info):
            list_lock.release()
            break
        else:
            info = thread_info.pop(0)
            list_lock.release()
        _cmd.coordset_update_thread(_self._COb,info)

def _coordset_update_spawn(thread_info,n_thread,_self=cmd):
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

def _object_update_thread(list_lock,thread_info,_self=cmd):
    # WARNING: internal routine, subject to change
    while 1:
        list_lock.acquire()
        if not len(thread_info):
            list_lock.release()
            break
        else:
            info = thread_info.pop(0)
            list_lock.release()
        _cmd.object_update_thread(_self._COb,info)

def _object_update_spawn(thread_info,n_thread,_self=cmd):
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

# do command (while API already locked)

def _do(cmmd,log=0,echo=1,_self=cmd):
    return _cmd.do(_self._COb,cmmd,log,echo)

# movie rendering

def _mpng(prefix, first=-1, last=-1, preserve=0, modal=0,
          format=-1, mode=-1, quiet=1,
          width=0, height=0,
          _self=cmd): # INTERNAL
    format = int(format)
    # WARNING: internal routine, subject to change
    try:
        _self.lock(_self)
        fname = prefix
        if re.search("[0-9]*\.png$",fname): # remove numbering, etc.
            fname = re.sub("[0-9]*\.png$","",fname)
        if re.search("[0-9]*\.ppm$",fname):
            if format<0:
                format = 1 # PPM
            fname = re.sub("[0-9]*\.ppm$","",fname)
        if format<0:
            format = 0 # default = PNG
        fname = cmd.exp_path(fname)
        r = _cmd.mpng_(_self._COb,str(fname),int(first),
                       int(last),int(preserve),int(modal),
                       format,int(mode),int(quiet),
                       int(width), int(height))
    finally:
        _self.unlock(-1,_self)
    return r

# copy image

def _copy_image(_self=cmd,quiet=1):
    # cmd._copy_image may be monkey-patched by GUI implementations
    raise NotImplementedError


# loading

def file_read(finfo, _self=cmd):
    '''
    Read a file, possibly gzipped or bzipped, and return the
    uncompressed file contents as a string.

    finfo may be a filename, URL or open file handle.
    '''
    try:
        if not is_string(finfo):
            handle = finfo
        elif '://' in finfo:
            req = urllib2.Request(finfo,
                    headers={'User-Agent': 'PyMOL/' + _self.get_version()[0]})
            handle = urllib2.urlopen(req)
        else:
            handle = open(finfo, 'rb')
        contents = handle.read()
        handle.close()
    except IOError:
        raise pymol.CmdException('failed to open file "%s"' % finfo)

    if contents[:2] == b'\x1f\x8b': # gzip magic number
        import io, gzip
        fakestream = io.BytesIO(contents)
        return gzip.GzipFile(fileobj=fakestream).read()

    if contents[:2] == b'BZ' and contents[4:10] == b'1AY&SY': # bzip magic
        import bz2
        return bz2.decompress(contents)

    return contents

def download_chem_comp(resn, quiet=1, _self=cmd):
    '''
    WARNING: internal routine, subject to change

    Download the chemical components CIF for the given residue name
    and return its local filename, or an empty string on failure.
    '''
    filename = os.path.join(_self.get('fetch_path'), resn + ".cif")
    if os.path.exists(filename):
        return filename

    url = "ftp://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/files/mmcif/" + resn + ".cif"
    url = "http://files.rcsb.org/ligands/download/" + resn + ".cif"
    if not quiet:
        print(' Downloading ' + url)

    try:
        contents = _self.file_read(url)
        if not contents: raise
    except:
        print(' Error: Download failed')
        return ''

    try:
        with open(filename, 'wb') as handle:
            handle.write(contents)
    except IOError as e:
        print(e)
        print('Your "fetch_path" setting might point to a read-only directory')
        return ''

    if not quiet:
        print('  ->' + filename)

    return filename

def _load(oname,finfo,state,ftype,finish,discrete,
          quiet=1,multiplex=0,zoom=-1,mimic=1,
          plugin='',
          object_props=None,
          atom_props=None, _self=cmd):
    # WARNING: internal routine, subject to change
    # caller must already hold API lock
    # NOTE: state index assumes 1-based state
    r = DEFAULT_ERROR
    contents = None
    size = 0
    if ftype not in (loadable.model,loadable.brick):
        if True:
            if ftype in _load2str:
                contents = _self.file_read(finfo)
                ftype = _load2str[ftype]
        return _cmd.load(_self._COb, str(oname), str(finfo), contents,
                          int(state) - 1, int(ftype),
                          int(finish),int(discrete),int(quiet),
                          int(multiplex),int(zoom), plugin,
                          object_props, atom_props, int(mimic))
    else:
        try:
            x = chempy.io.pkl.fromFile(finfo)
            if isinstance(x, (list, tuple)):
                for a in x:
                    r = _cmd.load_object(_self._COb,str(oname),a,int(state)-1,
                                                int(ftype),0,int(discrete),int(quiet),
                                                int(zoom))
                    if(state>0):
                        state = state + 1
                _cmd.finish_object(_self._COb,str(oname))
            else:
                r = _cmd.load_object(_self._COb,str(oname),x,
                                            int(state)-1,int(ftype),
                                            int(finish),int(discrete),
                                            int(quiet),int(zoom))
        except:
#            traceback.print_exc()
            print("Load-Error: Unable to load file '%s'." % finfo)
    return r

# function keys and other specials

modifier_keys = [
    '',
    'SHFT',
    'CTRL',
    'CTSH',
    'ALT',
]

special_key_codes = {
    # GLUT special key codes (see glutSpecialFunc)

    1        :  'F1',
    2        :  'F2',
    3        :  'F3',
    4        :  'F4',
    5        :  'F5',
    6        :  'F6',
    7        :  'F7',
    8        :  'F8',
    9        :  'F9',
    10       :  'F10',
    11       :  'F11',
    12       :  'F12',

    100      :  'left',
    101      :  'up',
    102      :  'right',
    103      :  'down',
    104      :  'pgup',
    105      :  'pgdn',
    106      :  'home',
    107      :  'end',
    108      :  'insert',
}

special_key_names = set(special_key_codes.values())

def _invoke_key(key, quiet=0, _self=cmd):
    '''Invoke a function that was mapped with cmd.set_key()'''
    try:
        mapping = _self.key_mappings[key]
    except KeyError:
        mapping = None

    if not mapping:
        if not quiet:
            print(" No key mapping for '%s'" % (key))
        return False

    if is_string(mapping):
        _self.do(mapping)
    else:
        fn, args, kwargs = mapping
        fn(*args, **kwargs)

    return True

def _special(k,x,y,m=0,_self=cmd): # INTERNAL (invoked when special key is pressed)
    pymol=_self._pymol
    # WARNING: internal routine, subject to change
    k=int(k)
    m=int(m)

    # convert numeric codes to string key

    try:
        key = special_key_codes[k]

        if m:
            key = modifier_keys[m] + '-' + key
    except KeyError:
        return False

    # check for explicit mapping

    if _invoke_key(key, 1, _self):
        return True

    # check for scenes and views

    for (fn, sc) in [
            (_self.scene, Shortcut(_self.get_scene_list())),
            (_self.view,  pymol._view_dict_sc),
            ]:
        if key in sc.keywords:
            fn(key)
            return True

        autocomp = sc.interpret(key + '-')
        if is_string(autocomp):
            fn(autocomp)
            return True

    print(" No key mapping and no scene or view for '%s'" % (key))
    return False

# control keys

def _ctrl(k,_self=cmd):
    # WARNING: internal routine, subject to change
    _invoke_key('CTRL-' + k, 0, _self)

# alt keys

def _alt(k,_self=cmd):
    # WARNING: internal routine, subject to change
    _invoke_key('ALT-' + k.upper(), 0, _self)

# command (apple) keys

def _cmmd(k,_self=cmd):
    # WARNING: internal routine, subject to change
    # command-key on macs
    if k in _self.cmmd:
        ak = _self.cmmd[k]
        if ak[0] is not None:
            ak[0](*ak[1], **ak[2])
    return None

def _ctsh(k,_self=cmd):
    # WARNING: internal routine, subject to change
    _invoke_key('CTSH-' + k, 0, _self)


# quitting (thread-specific)

def _quit(code=0, _self=cmd):
    pymol=_self._pymol
    # WARNING: internal routine, subject to change
    _self.interrupt()
    try:
        _self.lock(_self)
        try: # flush and close log if possible to avoid threading exception
            if pymol._log_file is not None:
                try:
                    pymol._log_file.flush()
                except:
                    pass
                pymol._log_file.close()
                del pymol._log_file
        except:
            pass
        if _self.reaper is not None:
            try:
                _self.reaper.join()
            except:
                pass
        r = _cmd.quit(_self._COb, int(code))
    finally:
        _self.unlock(-1,_self)
    return r

# screen redraws (thread-specific)

def _refresh(swap_buffers=1,_self=cmd):  # Only call with GLUT thread!
    # WARNING: internal routine, subject to change
    r = None
    if _self.is_gui_thread():
        def func():
            with _self.lockcm:
                if swap_buffers:
                    r = _cmd.refresh_now(_self._COb)
                else:
                    r = _cmd.refresh(_self._COb)
                return r
        r = _self._call_with_opengl_context(func)
    else:
        with _self.lockcm:
            r = _cmd.refresh_later(_self._COb)
    return r

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

def _validate_color_sc(_self=cmd):
    # WARNING: internal routine, subject to change
    if _self.color_sc is None: # update color shortcuts if needed
        lst = _self.get_color_indices()
        names = [x[0] for x in lst]
        names.extend(['default', 'auto', 'current', 'atomic'])
        names.extend(_self.get_names_of_type('object:ramp'))
        _self.color_sc = Shortcut(names)

def _invalidate_color_sc(_self=cmd):
    # WARNING: internal routine, subject to change
    _self.color_sc = None

def _get_color_sc(_self=cmd):
    # WARNING: internal routine, subject to change
    _validate_color_sc(_self=_self)
    return _self.color_sc

def _get_feedback(_self=cmd): # INTERNAL
    # WARNING: internal routine, subject to change
    l = []
    if _self.lock_attempt(_self):
        try:
            r = _cmd.get_feedback(_self._COb)
            while r:
                l.append(r)
                r = _cmd.get_feedback(_self._COb)
        finally:
            _self.unlock(-1,_self)
    else:
        l = None
    return l

def _fake_drag(_self=cmd): # internal
    _self.lock(_self)
    try:
        _cmd.fake_drag(_self._COb)
    finally:
        _self.unlock(-1,_self)
    return 1

def _sdof(tx,ty,tz,rx,ry,rz,_self=cmd):
    _cmd._sdof(_self._COb,tx,ty,tz,rx,ry,rz)

# testing tools

# for comparing floating point numbers calculated using
# different FPUs and which may show some wobble...

def _dump_floats(lst,format="%7.3f",cnt=9):
    # WARNING: internal routine, subject to change
    c = cnt
    for a in lst:
        print(format%a, end=' ')
        c = c -1
        if c<=0:
            print()
            c=cnt
    if c!=cnt:
        print()

def _dump_ufloats(lst,format="%7.3f",cnt=9):
    # WARNING: internal routine, subject to change
    c = cnt
    for a in lst:
        print(format%abs(a), end=' ')
        c = c -1
        if c<=0:
            print()
            c=cnt
    if c!=cnt:
        print()
