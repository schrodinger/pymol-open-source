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

# parser.py
# Python parser module for PyMol
#

from __future__ import absolute_import

# Don't import __future__.print_function


class SecurityException(Exception):
    pass


SCRIPT_TOPLEVEL = 'toplevel'


if True:

    import pymol
    import traceback
    import collections
    import re
    import glob
    import sys
    import os

    from . import parsing
    from . import colorprinting
    from .cmd import _feedback,fb_module,fb_mask,exp_path

    QuietException = parsing.QuietException
    CmdException = pymol.CmdException

    py_delims = { '=' : 1, '+='  : 1, '-='  : 1, '*=' : 1,
                      '/=' :1, '//=' : 1, '%='  : 1, '&=' : 1,
                      '|=' :1, '^='  : 1, '>>=' : 1,'<<=' : 1,
                      '**=':1 }

    remove_lists_re = re.compile("\[[^\]]*\]")

    def complete_sc(st,sc,type_name,postfix, mode=0):
        result = None

        try:
            sc=sc() # invoke lambda functions (if any)
        except:
            traceback.print_exc()
        amb = sc.interpret(st, mode)
        if amb is None:
            colorprinting.warning(" parser: no matching %s."%type_name)
        elif isinstance(amb, str):
            result = amb+postfix
        else:
            amb.sort()
            colorprinting.suggest(" parser: matching %s:"%type_name)
            flist = [x for x in amb if x[0]!='_']
            lst = parsing.list_to_str_list(flist)
            for a in lst:
                colorprinting.suggest(a)
            # now append up to point of ambiguity
            if not len(flist):
                css = []
            else:
                css = list(flist[0]) # common sub-string (css)
                for a in flist:
                    ac = list(a)
                    tmp = css
                    css = []
                    for c in range(len(tmp)):
                        if tmp[c]!=ac[c]:
                            break
                        css.append(tmp[c])
            css = [_f for _f in css if _f]
            css = ''.join(css)
            if len(css)>len(st):
                result = css
        return result

    class NestLayer:

        def __init__(self):

            self.cont = ""
            self.com0 = ""
            self.sc_path = SCRIPT_TOPLEVEL
            self.lineno = 0
            self.literal_python_fallback = False
            self.embed_sentinel = None
            self.embed_dict = {}
            self.next = []

    class Parser:

        def __init__(self,cmd):
            cmd = cmd._weakrefproxy

            self.cmd = cmd
            self.nest = 0
            self.layer = collections.defaultdict(NestLayer)
            self.pymol_names = self.cmd._pymol.__dict__

            # parsing state implemented with dictionaries to enable safe recursion
            # to arbitrary depths

            #com0 = {}    # verbose line, as read in
            #com1 = {}    # line w/o trailing whitespace
            #com2 = {}    # non-compound command
            #cont = {}    # continued characters from previous lines (i.e., before \ )
            #script = {}  # file handles
            #sc_path = {}   # file paths
            #kw = {}      # row in the keyword table for the current command
            #input = {}   # list of length two - command and unparsed arguments string
            #next = {}    # characters for next command (i.e., after ; )
            #args = {}    # parsed non-keyword argument string
            #kw_args = {} # parser keyword argument string
            #embed_dict = {}
            #embed_list = {}
            #embed_sentinel = {}
            #embed_type = {}
            #embed_line = {}

            # The resulting value from a pymol command (if any) is stored in the
            # parser.result global variable.  However, script developers will
            # geerally want to switch to the Python API for any of this kind of
            # stuff.

            self.result = None

            # initialize parser

            self.cmd._pymol.__script__ = SCRIPT_TOPLEVEL

        def exec_python(self, s, secure=False, fallback=False):
            if secure:
                raise SecurityException('Python expressions disallowed in this file')

            layer = self.layer[self.nest]
            layer.literal_python_fallback = fallback

            # for meaningful line number in error messages
            blanklines = layer.lineno - 1 - s.count('\n')

            s = '\n' * blanklines + s + '\n'
            s = compile(s, layer.sc_path, 'exec')
            exec(s, self.pymol_names, self.pymol_names)

        # main parser routine

        def parse(self,s,secure=0):
            try:
                self.nest += 1
                return self._parse(s, secure)
            finally:
                self.nest -= 1

        def _parse(self, s, secure):
            layer = self.layer[self.nest]
            self.result = None
# report any uncaught errors...
# WLD: this is problematic if parse is called inside an exception...removed.
#            if sys.exc_info()!=(None,None,None):
#                traceback.print_exc()
#                sys.exc_clear()
            def parse_embed():
                if s.strip() == layer.embed_sentinel:
                    etn = layer.embed_type
                    if etn == 0: # embedded data
                        colorprinting.parrot(" Embed: read %d lines."%(len(layer.embed_list)))
                        layer.embed_sentinel=None
                    elif etn == 1: # python block
                        colorprinting.parrot("PyMOL>"+s.rstrip())
                        py_block = ''.join(layer.embed_list)
                        del layer.embed_list
                        layer.embed_sentinel=None
                        self.exec_python(py_block)
                    elif etn == 2: # skip block
                        colorprinting.parrot(" Skip: skipped %d lines."%(layer.embed_line))
                        layer.embed_sentinel=None
                else:
                    etn = layer.embed_type
                    if etn == 0: # normal embedded data
                        layer.embed_list.append(s.rstrip()+"\n")
                    elif etn == 1: # python block
                        el = layer.embed_line + 1
                        colorprinting.parrot("%5d:%s"%(el,s.rstrip()))
                        layer.embed_line = el
                        layer.embed_list.append(s.rstrip()+"\n")
                    elif etn == 2:
                        layer.embed_line = layer.embed_line + 1

            p_result = 1
            layer.com0 = s
            try:
                if layer.embed_sentinel is not None:
                    parse_embed()
                    return 1

                layer.com1 = layer.com0.rstrip() # strips trailing whitespace
                if len(layer.com1) > 0:
                    if str(layer.com1[-1]) == "\\":
                        # prepend leftovers
                        if layer.cont != '':
                            layer.cont = layer.cont + "\n" + layer.com1[:-1]
                        else:
                            layer.cont = layer.com1[:-1]
                    else:
                        # prepend leftovers
                        if layer.cont != '':
                            layer.com1 = layer.cont + "\n" + layer.com1
                            layer.cont = ''
        # this routine splits up the line first based on semicolon

                        layer.next = parsing.split(layer.com1,';',1) + layer.next[1:]

        # layer.com2 now a full non-compound command
                        layer.com2 = layer.next[0]
                        layer.input = layer.com2.split(' ',1)
                        lin = len(layer.input)
                        if lin:
                            layer.input[0] = layer.input[0].strip()
                            com = layer.input[0]
                            if (com[0:1]=='/'):
                                # explicit literal python
                                layer.com2 = layer.com2[1:].strip()
                                if len(layer.com2)>0:
                                    self.exec_python(layer.com2, secure)
                            elif lin>1 and layer.input[-1:][0].split(' ',1)[0] in py_delims:
                                self.exec_python(layer.com2, secure)
                            else:
                                # try to find a keyword which matches
                                if com in self.cmd.kwhash:
                                    amb = self.cmd.kwhash.interpret(com)
                                    if amb is None:
                                        com = self.cmd.kwhash[com]
                                    elif not isinstance(amb, str):
                                        colorprinting.warning('Error: ambiguous command: ')
                                        amb.sort()
                                        amb = parsing.list_to_str_list(amb)
                                        for a in amb:
                                            colorprinting.warning(a)
                                        raise QuietException
                                    com = amb
                                if com in self.cmd.keyword:
                                    # here is the command and argument handling section
                                    layer.kw = self.cmd.keyword[com]
                                    if layer.kw[4]>=parsing.NO_CHECK:
                                        # stricter, Python-based argument parsing
                                        # remove line breaks (only important for Python expressions)
                                        layer.com2=layer.com2.replace('\n','')

                                        if layer.kw[4]>=parsing.LITERAL: # treat literally
                                            layer.next = []
                                            if not secure:
                                                layer.com2=layer.com1
                                            else:
                                                raise SecurityException('Python expressions disallowed in this file')
                                        if secure and (layer.kw[4]==parsing.SECURE):
                                            layer.next = []
                                            raise SecurityException('Command disallowed in this file')
                                        else:
                                           (layer.args, layer.kw_args) = \
                                            parsing.prepare_call(
                                             layer.kw[0],
                                             parsing.parse_arg(layer.com2,mode=layer.kw[4],_self=self.cmd),
                                             layer.kw[4], com, _self=self.cmd) # will raise exception on failure
                                        self.result=layer.kw[0](*layer.args, **layer.kw_args)
                                    elif layer.kw[4]==parsing.PYTHON:
                                            # handle python keyword
                                            layer.com2 = layer.com2.strip()
                                            if len(layer.com2)>0:
                                                self.exec_python(layer.com2, secure)
                                    else:
                                        # remove line breaks (only important for Python expressions)
                                        layer.com2=layer.com2.replace('\n','')
                                        # old parsing style, being phased out
                                        if layer.kw[4]==parsing.ABORT:
                                            return None # SCRIPT ABORT EXIT POINT
                                        if layer.kw[4]==parsing.MOVIE: # copy literal single line, no breaks
                                            layer.next = []
                                            if not secure:
                                                layer.input = layer.com1.split(' ',1)
                                            else:
                                                raise SecurityException('Movie commands disallowed in this file')
                                        if len(layer.input)>1:
                                            layer.args = parsing.split(layer.input[1],layer.kw[3])
                                            while 1:
                                                nArg = len(layer.args) - 1
                                                c = 0
                                                while c < nArg:
                                                    if ( layer.args[c].count('(')!=
                                                         layer.args[c].count(')')):
                                                        tmp=layer.args[c+1]
                                                        layer.args.remove(tmp)
                                                        layer.args[c]=layer.args[c].strip()+\
                                                                          ','+tmp.strip()
                                                        nArg = nArg-1
                                                        break;
                                                    c = c + 1
                                                if c == nArg:
                                                    break;
                                            if len(layer.args)==1 and len(layer.args[0])==0:
                                                layer.args = []
                                        else:
                                            layer.args = []
                                        if layer.kw[1]<= len(layer.args) <= layer.kw[2]:
                                            layer.args = [a.strip() for a in layer.args]
                                            if layer.kw[4]<parsing.RUN:
                                                #
                                                # this is where old-style commands are invoked
                                                #
                                                self.result=layer.kw[0](*layer.args)
                                                #
                                            elif (layer.kw[4]==parsing.EMBED):
                                                layer.next = []
                                                if secure or self.nest==0: # only legal on top level and p1m files
                                                    l = len(layer.args)
                                                    if l>0:
                                                        key = layer.args[0]
                                                    else:
                                                        key = self.get_default_key()
                                                    if l>1:
                                                        format = layer.args[1]
                                                    else:
                                                        format = 'pdb'
                                                    if l>2:
                                                        layer.embed_sentinel = layer.args[2]
                                                    else:
                                                        layer.embed_sentinel = "embed end"
                                                    list = []
                                                    layer.embed_dict[key] = ( format, list )
                                                    layer.embed_list = list
                                                    layer.embed_type = 0 # not a python block
                                                else:
                                                    print('Error: embed only legal in special files (e.g. p1m)')
                                                    raise None
                                            elif (layer.kw[4]==parsing.SKIP):
                                                layer.next = []
                                                arg = parsing.apply_arg(
                                                    parsing.parse_arg(layer.com2,_self=self.cmd),
                                                    ('sentinel',),
                                                    {'sentinel':'skip end'})
                                                print(arg) # ???
                                                if len(layer.args):
                                                    if layer.args[0]=='end': # probable 'skip end' to ignore
                                                        arg = []
                                                if len(arg):
                                                    layer.embed_sentinel = arg[0]
                                                    layer.embed_type = 2 # skip block
                                                    layer.embed_line = 0
                                            elif (layer.kw[4]==parsing.PYTHON_BLOCK):
                                                layer.next = []
                                                if not secure:
                                                    arg = parsing.apply_arg(
                                                        parsing.parse_arg(layer.com2,_self=self.cmd),
                                                        ('sentinel','skip'),
                                                        {'sentinel':'python end','skip':0})
                                                    layer.embed_sentinel = arg[0]
                                                    list = []
                                                    layer.embed_list = list
                                                    if arg[1]:
                                                        layer.embed_type = 2 # skip block
                                                    else:
                                                        layer.embed_type = 1 # python block
                                                    layer.embed_line = 0
                                                else:
                                                    print('Error: Python blocks disallowed in this file.')
                                                    raise None
                                            else:
                                                print('Error: unknown keyword mode: '+str(layer.kw[4]))
                                                raise QuietException
                                        else:
                                            print('Error: invalid arguments for %s command.' % com)
        #
        # non-keyword command handling
        #
                                elif len(layer.input[0]):
                                    if layer.input[0][0]=='@':
                                        path = exp_path(layer.com2[1:].strip())
                                        if path[-3:].lower()=='p1m':
                                            nest_securely = 1
                                        else:
                                            nest_securely = secure
                                        if re.search("\.py$|\.pym$",path) is not None:
                                            if self.cmd._feedback(fb_module.parser,fb_mask.warnings):
                                                print("Warning: use 'run' instead of '@' with Python files?")
                                        layer.script = open(path,'r')
                                        self.cmd._pymol.__script__ = path
                                        self.nest=self.nest+1
                                        self.layer[self.nest] = NestLayer()
                                        layer = self.layer[self.nest]
                                        layer.cont=''
                                        layer.sc_path=path
                                        layer.embed_sentinel=None
                                        while 1:
                                            layer.com0  = self.layer[self.nest-1].script.readline()
                                            self.layer[self.nest].lineno += 1
                                            if not layer.com0: break
                                            inp_cmd = layer.com0
                                            tmp_cmd = inp_cmd.strip()
                                            if len(tmp_cmd):
                                                if tmp_cmd[0] not in ['#','_','/']: # suppress comments, internals, python
                                                    if layer.embed_sentinel is None:
                                                        colorprinting.parrot("PyMOL>"+tmp_cmd)
                                                elif tmp_cmd[0]=='_' and \
                                                      tmp_cmd[1:2] in [' ','']: # "_ " remove echo suppression signal
                                                    inp_cmd=inp_cmd[2:]
                                            pp_result = self.parse(inp_cmd,nest_securely)
                                            if pp_result is None: # RECURSION
                                                break # abort command gets us out
                                            elif pp_result==0: # QuietException
                                                if self.cmd.get_setting_boolean("stop_on_exceptions"):
                                                    p_result = 0 # signal an error occurred
                                                    colorprinting.error("PyMOL: stopped on exception.")
                                                    break;
                                        self.nest=self.nest-1
                                        layer=self.layer[self.nest]

                                        layer.script.close()
                                        self.cmd._pymol.__script__ = layer.sc_path
                                    else: # nothing found, try literal python
                                        layer.com2 = layer.com2.strip()
                                        if len(layer.com2)>0:
                                            if not secure:
                                                self.exec_python(layer.com2, fallback=True)
                                            elif layer.input[0][0:1]!='#':
                                                colorprinting.error('Error: unrecognized keyword: '+layer.input[0])
                        if (len(layer.next)>1) and p_result:
                            # continue parsing if no error or break has occurred
                            self.nest=self.nest+1
                            self.layer[self.nest] = NestLayer()
                            layer=self.layer[self.nest]
                            layer.com0 = self.layer[self.nest-1].next[1]
                            self.layer[self.nest-1].next=[]
                            layer.cont=''
                            layer.embed_sentinel=None
                            p_result = self.parse(layer.com0,secure) # RECURSION
                            self.nest=self.nest-1
                            layer=self.layer[self.nest]
            except (QuietException, CmdException) as e:
                if e.args:
                    colorprinting.error(e)
                if self.cmd._feedback(fb_module.parser,fb_mask.blather):
                    print("Parser: caught " + type(e).__name__)
                p_result = 0
            except SecurityException as e:
                colorprinting.error('Error: %s' % (e,))
                p_result = None
            except:
                exc_type, exc_value, tb = colorprinting.print_exc(
                        [__file__, SCRIPT_TOPLEVEL])

                p_result = 0 # notify caller that an error was encountered
            if not p_result and self.cmd._pymol.invocation.options.exit_on_error:
                self.cmd.quit(1)
            return p_result  # 0 = Exception, None = abort, 1 = ok

        def get_embedded(self,key=None):
            layer = self.layer[self.nest]
            dict = layer.embed_dict
            if key is None:
                key = self.get_default_key()
            return dict.get(key,None)

        def get_default_key(self):
            layer = self.layer[self.nest]
            return os.path.splitext(os.path.basename(layer.sc_path))[0]

        def stdin_reader(self): # dedicated thread for reading standard input
            import sys
            while 1:
                try:
                    l = sys.stdin.readline()
                except IOError:
                    continue
                if l!="":
                    if self.nest==0:
                        # if we're reading embedded input on stdin
                        # then bypass PyMOL C code altogether
                        if self.layer[0].embed_sentinel is not None:
                            self.parse(l)
                        else:
                            self.cmd.do(l, flush=True)
                    else:
                        self.cmd.do(l, flush=True)
                elif not self.cmd._pymol.invocation.options.keep_thread_alive:
                    self.cmd.quit()
                else:
                    import time
                    time.sleep(.1)

            self.cmd._pymol._stdin_reader_thread = None

        def complete(self,st):
            with self.cmd.lockcm:
                return self._complete(st)

        def _complete(self,st):
            result = None
            pre = ''
            flag = 0
            if not (' ' in st or '@' in st):
                try:
                    result = complete_sc(st, self.cmd.kwhash, 'commands',' ', 1)
                except:
                    traceback.print_exc()
            else:
                full = self.cmd.kwhash.interpret(re.sub(r" .*","",st))
                st_no_lists = remove_lists_re.sub("",st)
                count = st_no_lists.count(',') # which argument are we on
                if self.cmd.is_string(full):
                    try:
                        if count<len(self.cmd.auto_arg):
                            if full in self.cmd.auto_arg[count]: # autocomplete arguments
                                flag = 1
                                pre = re.sub(r"^[^ ]* ",' ',st,count=1) # trim command
                                if re.search(r",",pre) is not None:
                                    pre = re.sub(r"[^\, ]*$","",pre,count=1)
                                    pre = re.sub(r",\s*[^\, ]*$",", ",pre,count=1) # trim 1 arg
                                else:
                                    pre = re.sub("[^ ]*$","",pre,count=1) # trim 1 arg
                                pre = re.sub(r"^ *",'',pre)
                                pre = full+' '+pre
                                pat = re.sub(r".*[\, ]",'',st)
            #               print ":"+pre+":"+pat+":"
#                                print tuple([pat] + self.cmd.auto_arg[count][full])
                                result = complete_sc(*tuple([pat] + self.cmd.auto_arg[count][full]), **{})
                    except:
                        traceback.print_exc()
                if not flag: # otherwise fallback onto filename completion
                    st = self.cmd.as_pathstr(st)
                    loc = 1 + max(map(st.rfind, ',@'))
                    if not loc:
                        loc = 1 + st.find(' ')
                    pre = st[:loc]
                    st3 = st[loc:].lstrip()
                    flist = glob.glob(exp_path(st3)+"*")

                    # environment variable completion
                    if not flist and st3.startswith('$'):
                        flist = ['$' + var for var in os.environ
                                if var.startswith(st3[1:])]

                    lf = len(flist)
                    if lf == 0:
                        print(" parser: no matching files.")
                    elif lf==1:
                        result = flist[0]
                        if os.path.isdir(flist[0]):
                            result += '/' # do not use os.path.sep here
                    else:
                        flist.sort()
                        print(" parser: matching files:")
                        lst = parsing.list_to_str_list(flist)
                        for a in lst:
                            print(a)
                        # now append as much up to point of ambiguity
                        css = os.path.commonprefix(flist)
                        if len(css)>len(st3):
                            result = css
            if result is not None:
                result = pre+result
            return result

    def new_parse_closure(self_cmd): # create parser and return an instance-specific parse function closure
        try:
            p = Parser(self_cmd)
        except:
            traceback.print_exc()
        self_cmd._parser = p
        return lambda s,secure,p=p:p.parse(s,secure)

    def new_complete_closure(self_cmd): # return an instance-specific complete function closure
        return lambda st,p=self_cmd._parser:p.complete(st)

# unused code?
#
#    def _same_(a,b):
#        if a==b:
#            return a
#        else:
#            return None
