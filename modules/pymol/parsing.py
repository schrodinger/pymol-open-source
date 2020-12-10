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

# parser2.py
# An improved command parser for PyMOL, but still a terrible kluuge
# PyMOL needs to migrate to a real parser soon!

# === Goals:
#  1. improved 1to1 mapping between pymol "cmd" API and command language
#  2. support for named arguments
#  3. support for calling arbitrary python functions via this mapping

# === Syntatic Examples

# * simple commands
# command

# * commands with arguments
#
# command value1
# command value1,value2
# command value1,value2,value3

# * commands with named arguments
#
# command argument1=value1
# command argument1=value1,argument2=value2,argument3=value3
# command argument3=value1,argument2=value2,argument1=value1

# * mixed...
#
# command value1,argument3=value3,arg

# * commands with legacy '=' support for first argument
#
# command string1=value1
# * which should map to
# command string1,value1

# * legacy '=' combined as above
#
# command string1=value1,value2
# command string1=value1,value2,value3,...
# command string1=value1,argument3=value3

# === Burdens placed on API functions...
#
# None. However, function must have real arguments for error checking.
#

if True:

    import re
    import sys
    import threading
    import types
    import inspect
    from . import colorprinting

    class QuietException(Exception):
        pass

    # constants for keyword modes

    SIMPLE      = 0  # original pymol parsing (deprecated)
    MOVIE       = 1  # ignore ";", treat entire line as a single command
    RUN         = 2  # run command
    SPAWN       = 3  # for spawn and fork commands
    ABORT       = 4  # terminates command script
    PYTHON      = 5  # pass entire line to python
    EMBED       = 6  # embedded data
    PYTHON_BLOCK = 7 # embedded python block
    SKIP        = 8  # skipping commands
    NO_CHECK    = 10 # no error checking
    STRICT      = 11 # strict name->argument checking
    SECURE      = 12 # command not available in "secure" mode
    LEGACY      = 13 # support legacy construct str1=val1,... -> str1,val1,...
    LITERAL     = 20 # argument is to be treated as a literal string
    LITERAL1    = 21 # one regular argument, followed by literal string
    LITERAL2    = 22 # two regular argument, followed by literal string

    # key regular expressions

    arg_name_re = re.compile(r"[A-Za-z0-9_]+\s*\=")
    nester_char_re = re.compile(r"\(|\)|\[|\]")
    nester_re = re.compile(r"[^,;]*[\(\[]")
    arg_pre_nester_re = re.compile(r"([^,;\(\[]+)[\(\[]")
    arg_post_nester_re = re.compile(r"[^,;\(\[]*")
    arg_easy_nester_re = re.compile(r"\([^,]*\)|\[[^,]*\]")
    arg_hard_nester_re = re.compile(r"\(.*\)|\[.*\]")
    # NOTE '''sdf'sdfs''' doesn't work in below.
    arg_value_re = re.compile(r"'''[^']*'''|'[^']*'|"+r'"[^"]*"|[^,;]+')
    def trim_nester(st):
        # utility routine, returns single instance of a nested string
        # should be modified to handle quotes too
        pc = 1
        l = len(st)
        c = 1
        while c<l:
            if st[c] in ('(','['):
                pc = pc + 1
            if st[c] in (')',']'):
                pc = pc - 1
            c = c + 1
            if not pc:
                break
        if pc:
            return None
        return st[0:c]

    def apply_arg(inp_arg,par=(),def_dict={}):
        n_inp = len(inp_arg)
        n_req = n_inp - len(def_dict)
        result = []
        inp_dict = {}
        for a in inp_arg:
            if a[0] is not None:
                inp_dict[a[0]] = a[1];
        c = 0
        for p in par:
            if c<n_inp:
                a = inp_arg[c]
                if a[0] is None:
                    result.append(a[1])
                    c = c + 1
                    continue
            if p in inp_dict:
                result.append(inp_dict[p])
                del inp_dict[p]
            elif p in def_dict:
                result.append(def_dict[p])
            elif c<n_req:
                raise QuietException("Error: invalid argument(s).")
            c = c + 1
        if len(inp_dict):
            raise QuietException("Error: invalid argument(s).")
        return result

    def parse_arg(st,mode=STRICT,_self=None):
        '''
    parse_arg(st)

    expects entire command to be passed in

    returns list of tuples of strings: [(None,value),(name,value)...]
    '''
        result = []
        # current character
        cc = 0
        a = st.split(None, 1)
        if len(a) == 2:
            st = a[1]
            while 1:
                if mode>=LITERAL: # LITERAL argument handling
                    if (mode-LITERAL)==len(result):
                        result.append((None, st[cc:].strip()))
                        return result
                # clean whitespace
                st = st.lstrip()
                if st == '':
                    break
                # read argument name, if any
                mo = arg_name_re.match(st)
                if mo:
                    nam = mo.group(0)[:-1].strip()
                    st = st[mo.end(0):].lstrip()
                else:
                    nam = None
                # is one or more nesters present?
                skip_flag = 0
                if nester_re.match(st[cc:]):
                    skip_flag = 1
                    nest_flag = 1
                    nest_str = ''
                    while nest_flag: # parse all the nesters
                        nest_flag = 0
                        # text before nester?
                        mo = arg_pre_nester_re.match(st[cc:])
                        if mo:
                            nest_str = nest_str + mo.group(1)
                            cc=cc+mo.end(1)
                        # special handling for nesters (selections, lists, tuples, etc.)
                        mo = arg_easy_nester_re.match(st[cc:]) # no internal commas
                        if mo:
                            cnt = len(nester_char_re.findall(mo.group(0)))
                            if cnt % 2 == 1: # make sure nesters are matched in count
                                mo = None
                        if mo:
                            nest_str = nest_str + mo.group(0)
                            cc=cc+mo.end(0)
                            # text after nester?
                            mo = arg_post_nester_re.match(st[cc:])
                            if mo:
                                post_nester = mo.group(0)
                                cc=cc+mo.end(0)
                            nest_str = nest_str + post_nester
                            nest_flag = 1 # one more cycle
                        else:
                            mo = arg_hard_nester_re.match(st[cc:])
                            if mo:
                                se = trim_nester(mo.group(0))
                                if se is None:
                                    colorprinting.error("Error: "+st)
                                    colorprinting.error("Error: "+" "*cc+"^ syntax error (type 1).")
                                    raise QuietException
                                else:
                                    cc = cc + len(se)
                                    nest_str = nest_str + se
                                    # text after nester?
                                    mo = arg_post_nester_re.match(st[cc:])
                                    if mo:
                                        nest_str = nest_str + mo.group(0)
                                        cc=cc+mo.end(0)
                                    nest_flag = 1 # one more cycle
                    if not len(nest_str): # we must have failed to parse...
                        skip_flag = 0
                    else:
                        result.append((nam, nest_str.strip()))
                if not skip_flag:
                    # no nester, so just read normal argument value
                    argval = None
                    mo = arg_value_re.match(st[cc:])
                    if not mo:
                        if(st[cc:cc+1]!=','):
                            colorprinting.error("Error: "+st)
                            colorprinting.error("Error: "+" "*cc+"^ syntax error (type 2).")
                            raise QuietException
                        else:
                            # allow blank arguments
                            result.append((nam,None))
                    else:
                        argval = mo.group(0)
                        cc=cc+mo.end(0)
                        while 1: # pickup unqouted characters after quotes
                            mo = arg_value_re.match(st[cc:])
                            if not mo:
                                break
                            argval = argval + mo.group(0)
                            cc=cc+mo.end(0)
                        if argval is not None:
                            result.append((nam, argval.strip()))
                # clean whitespace
                st = st[cc:].lstrip()
                cc = 0
                # skip over comma
                if st != '':
                    if st.startswith(','):
                        st = st[1:].lstrip()
                    else:
                        colorprinting.error("Error: "+st)
                        colorprinting.error("Error: "+" "*cc+"^ syntax error (type 3).")
                        raise QuietException
        if __name__!='__main__':
            if _self._feedback(_self.fb_module.parser, _self.fb_mask.debugging):
                _self.fb_debug.write(" parsing-DEBUG: tup: "+str(result)+"\n")
        return result

    def dump_str_list(list):
        lst = list_to_str_list(list)
        for a in lst:
            print(a)

    def list_to_str_list(list,width=77,margin=2): # format strings into a list
        result = []
        ll=len(list)
        if ll>0:
            mxln = 1
            for a in list:
                if len(a)>mxln:
                    mxln = len(a)
            n_col = width//mxln
            width = width - margin
            while (n_col * mxln + n_col*2)>width:
                n_col = n_col - 1
            if n_col < 1:
                n_col = 1
            ll = len(list)
            n_row = len(list)//n_col
            while (n_row*n_col)<ll:
                n_row = n_row + 1
            rows = []
            for a in range(n_row):
                rows.append([])
            row = 0
            pad_list = []
            for a in list:
                pad_list.append(("%-"+str(mxln)+"s")%a)
            for a in pad_list:
                rows[row].append(a)
                row = row + 1
                if row >= n_row:
                    row = 0
            for a in rows:
                st = margin*' '
                row = 0
                st = st + '  '.join(a)
                result.append(st)
        return result

    def dump_arg(name,arg_lst,nreq):
        ac = 0
        pc = 0
        st = "Usage: "+name
        if '_self' in arg_lst:
            arg_lst = list(arg_lst)
            arg_lst.remove('_self')
        for a in arg_lst:
            if ac>=nreq:
                st = st + " ["
                pc = pc + 1
            if ac:
                st = st + ", " + a
            else:
                st = st + " " + a
            ac = ac + 1
        print(st + " " + "]"*pc)

    def prepare_call(fn, lst, mode=STRICT, name=None, _self=None):
        """
        Validates the command arguments. Handles legacy syntax
        (set name=value) and the "?" usage argument.

        @param fn Function
        @param lst Argument list with (key, value) items, key may be None for positional arguments
        @param mode Argument parsing mode
        @param name Name of the command (only for error reporting)
        @return Tuple of (args: list, kwargs: dict)
        @raises QuietException on error
        """
        if name is None:
            name = fn.__name__

        arg = []
        kw = {}
        co = fn.__code__

        # disable error checking for *arg or **kw functions
        if (co.co_flags & 0xC):
            mode = NO_CHECK

        offset = 1 if inspect.ismethod(fn) else 0
        arg_nam = co.co_varnames[offset:co.co_argcount + co.co_kwonlyargcount]
        npositional = co.co_argcount - offset
        nreq = max(0, npositional - len(fn.__defaults__ or ()))

        # co_posonlyargcount is new in Python 3.8
        nposonly = max(0, getattr(co, 'co_posonlyargcount', 0) - offset)

        assert nposonly <= npositional

        if lst == [(None, '?')]:
            dump_arg(name, arg_nam, nreq)
            raise QuietException

        if mode == NO_CHECK:
            # no error checking
            for (key, value) in lst:
                if key is None:
                    arg.append(value)
                else:
                    kw[key] = value

            # always send _self in the dictionary
            if "_self" not in kw:
                kw["_self"] = _self
        else:
            # error checking enabled

            # handle legacy string=value transformation
            if mode == LEGACY:
                tmp_lst = []
                for a in lst:
                    if a[0] is not None and a[0] not in arg_nam:
                        tmp_lst.append((None, a[0]))
                        tmp_lst.append((None, a[1]))
                    else:
                        tmp_lst.append(a)
                lst = tmp_lst

            # match names to unnamed arguments to create argument dictionary
            for ac, (key, value) in enumerate(lst):
                if key is None:
                    if ac < nposonly:
                        arg.append(value)
                        continue

                    if ac >= npositional:
                        dump_arg(name, arg_nam, nreq)
                        raise QuietException(
                            f"Error: too many positional arguments for {name}")
                    key = arg_nam[ac]
                kw[key] = value

            # now check to make sure we don't have any missing arguments
            for key in set(arg_nam[nposonly:nreq]).difference(kw):
                raise QuietException("Parsing-Error: missing required "
                                     f"argument in function {name} : {key}")

            # make sure command knows which PyMOL instance to message
            if "_self" not in kw and "_self" in arg_nam:
                kw["_self"] = _self

        # set feedback argument (quiet), if extant, results enabled, and not overridden
        if "quiet" not in kw and "quiet" in arg_nam and _self._feedback(
                _self.fb_module.cmd, _self.fb_mask.results):
            kw["quiet"] = 0

        return arg, kw


    # launching routines

    def run(filename, namespace='global', _spawn=0, _self=None):
        '''
DESCRIPTION

    "run" executes an external Python script in a local name space,
    the main Python namespace, the global PyMOL namespace, or in its
    own namespace (as a module).

USAGE

    run file [, namespace ]

ARGUMENTS

    file = string: a Python program, typically ending in .py or .pym.

    namespace = local, global, module, main, or private {default: global}

NOTES

    Due to an idiosyncracy in Pickle, you can not pickle objects
    directly created at the main level in a script run as "module",
    (because the pickled object becomes dependent on that module).
    Workaround: delegate construction to an imported module.

SEE ALSO

    spawn
        '''
        from __main__ import __dict__ as ns_main
        from pymol    import __dict__ as ns_pymol

        if not _self:
            from pymol import cmd as _self

        if filename.endswith('.pml'):
            return _self.load(filename)

        path = _self.exp_path(filename)
        spawn = int(_spawn)
        run_ = spawn_file if spawn else run_file

        if namespace == 'global':
            run_(path, ns_pymol, ns_pymol)
        elif namespace == 'local':
            run_(path, ns_pymol, {})
        elif namespace == 'main':
            run_(path, ns_main, ns_main)
        elif namespace == 'private':
            run_(path, ns_main, {})
        elif namespace == 'module':
            run_file_as_module(path, spawn=spawn)
        else:
            raise ValueError('invalid namespace "%s"' % namespace)

    def spawn(filename, namespace='module', _self=None):
        '''
DESCRIPTION

    "spawn" launches a Python script in a new thread which will run
    concurrently with the PyMOL interpreter. It can be run in its own
    namespace (like a Python module, default), a local name space, or
    in the global namespace.

USAGE

    spawn file [, namespace ]

NOTES

    The default namespace for spawn is "module".

    The best way to spawn processes at startup is to use the -l option
    (see "help launching").

SEE ALSO

    run
        '''
        return run(filename, namespace, 1, _self)

    def _print_exc():
        colorprinting.print_exc([__file__])

    def execfile(filename, global_ns, local_ns):
        import pymol.internal as pi
        co = compile(pi.file_read(filename), filename, 'exec')
        exec(co, global_ns, local_ns)

    def run_file(file,global_ns,local_ns):
        global_ns['__script__'] = file
        try:
            execfile(file,global_ns,local_ns)
        except pymol.CmdException:
            # so the idea here is to print the traceback here and then
            # cascade all the way back up to the interactive level
            # without any further output
            _print_exc()
            raise QuietException

    def run_file_as_module(file,spawn=0):
        name = re.sub('[^A-Za-z0-9]','_',file)
        mod = types.ModuleType(name)
        mod.__file__ = file
        mod.__script__ = file
        sys.modules[name]=mod
        if spawn:
            t = threading.Thread(target=execfile,
                args=(file,mod.__dict__,mod.__dict__))
            t.setDaemon(1)
            t.start()
        else:
            try:
                execfile(file,mod.__dict__,mod.__dict__)
            except pymol.CmdException:
                _print_exc()
                raise QuietException
            del sys.modules[name]
            del mod

    def spawn_file(args,global_ns,local_ns):
        local_ns['__script__'] = args
        t = threading.Thread(target=execfile,args=(args,global_ns,local_ns))
        t.setDaemon(1)
        t.start()

    def split(str, tok, mx=0):
        '''
    split(string,token[,count]) -> list of strings

    UTILITY FUNCTION, NOT PART OF THE API
    Breaks strings up by tokens but preserves quoted strings and
    parenthetical groups (such as atom selections).

    USAGE OF THIS FUNCTION IS DISCOURAGED - THE GOAL IS TO
    MAKE IT UNNECESSARY BY IMPROVING THE BUILT-IN PARSER
    '''
        pair = { '(':')','[':']','{':'}',"'":"'",'"':'"' }
        plst = list(pair.keys())
        stack = []
        lst = []
        c = 0
        nf = 0
        l = len(str)
        wd = ""
        while str[c]==tok:
            c = c + 1
        while c<l:
            ch = str[c]
            if (ch in tok) and (len(stack)==0):
                lst.append(wd.strip())
                nf = nf + 1
                if mx:
                    if nf==mx:
                        wd = str[c+1:].strip()
                        break;
                wd = ''
                w = 0
            else:
                if len(stack):
                    if ch==stack[0]:
                        stack = stack[1:]
                    elif (ch in plst):
                        stack[:0]=[pair[ch]]
                elif (ch in plst):
                    stack[:0]=[pair[ch]]
                wd = wd + ch
            c = c + 1
        if len(wd):
            lst.append(wd.strip())
        return lst

    import pymol
