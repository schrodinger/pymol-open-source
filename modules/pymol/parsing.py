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

if __name__=='pymol.parsing':

    import re
    import string
    import sys
    import threading
    import new
    import traceback
    import copy
    
    class QuietException:
        def __init__(self,args=None):
            self.args = args

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

    command_re = re.compile(r"[^\s]+")
    whitesp_re = re.compile(r"\s+")
    comma_re = re.compile(r"\s*,\s*")
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
            if a[0] != None:
                inp_dict[a[0]] = a[1];
        c = 0
        for p in par:
            if c<n_inp:
                a = inp_arg[c]
                if a[0] == None:
                    result.append(a[1])
                    c = c + 1
                    continue
            if inp_dict.has_key(p):
                result.append(inp_dict[p])
                del inp_dict[p]
            elif def_dict.has_key(p):
                result.append(def_dict[p])
            elif c<n_req:
                print "Error: invalid argument(s)."
                raise QuietException
            c = c + 1
        if len(inp_dict):
            print "Error: invalid argument(s)."            
            raise QuietException
        return result
    
    def parse_arg(st,mode=STRICT):
        '''
    parse_arg(st)

    expects entire command to be passed in

    returns list of tuples of strings: [(None,value),(name,value)...]
    '''
        result = [] 
        # current character
        cc = 0
        mo = command_re.match(st[cc:])
        if mo:
            cc=cc+mo.end(0)
            while 1:
                if mode>=LITERAL: # LITERAL argument handling
                    if (mode-LITERAL)==len(result):
                        result.append((None,string.strip(st[cc:])))
                        return result
                # clean whitespace
                mo = whitesp_re.match(st[cc:])
                if mo:
                    cc=cc+mo.end(0)            
                if not len(st[cc:]):
                    break
                # read argument name, if any         
                mo = arg_name_re.match(st[cc:])
                if mo:
                    nam = string.strip(mo.group(0)[:-1])
                    cc=cc+mo.end(0)
                else:
                    nam = None
                # clean whitespace
                mo = whitesp_re.match(st[cc:])
                if mo:
                    cc=cc+mo.end(0)
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
                            if (2*(cnt/2))!=cnt: # make sure nesters are matched in count
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
                                if se==None:
                                    print "Error: "+st
                                    print "Error: "+" "*cc+"^ syntax error (type 1)."
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
                        result.append((nam,string.strip(nest_str)))
                if not skip_flag:
                    # no nester, so just read normal argument value
                    argval = None
                    mo = arg_value_re.match(st[cc:])
                    if not mo:
                        if(st[cc:cc+1]!=','):
                            print "Error: "+st
                            print "Error: "+" "*cc+"^ syntax error (type 2)."
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
                        if argval!=None:
                            result.append((nam,string.strip(argval)))
                # clean whitespace
                mo = whitesp_re.match(st[cc:])
                if mo:
                    cc=cc+mo.end(0)
                # skip over comma
                if len(st[cc:]):
                    mo = comma_re.match(st[cc:])
                    if mo:
                        cc=cc+mo.end(0)
                    else:
                        print "Error: "+st
                        print "Error: "+" "*cc+"^ syntax error (type 3)."
                        raise QuietException
        if __name__!='__main__':
            if cmd._feedback(cmd.fb_module.parser,cmd.fb_mask.debugging):
                cmd.fb_debug.write(" parsing-DEBUG: tup: "+str(result)+"\n")
        return result

    def dump_str_list(list):
        lst = list_to_str_list(list)
        for a in lst:
            print a

    def list_to_str_list(list,width=77,margin=2): # format strings into a list
        result = []
        ll=len(list)
        if ll>0:
            mxln = 1
            for a in list:
                if len(a)>mxln:
                    mxln = len(a)
            n_col = width/mxln
            width = width - margin
            while (n_col * mxln + n_col*2)>width:
                n_col = n_col - 1
            if n_col < 1:
                n_col = 1
            ll = len(list)
            n_row = len(list)/n_col
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
                st = st + string.join(a,'  ')
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
        print st,"]"*pc

    def prepare_call(fn,lst,mode=STRICT,name=None): # returns tuple of arg,kw or excepts if error
        if name==None:
            name=fn.__name__
        result = (None,None)
        arg = []
        kw = {}
        co = fn.func_code
        if (co.co_flags & 0xC): # disable error checking for *arg or **kw functions
            mode = NO_CHECK
        arg_nam = co.co_varnames[0:co.co_argcount]
        narg = len(arg_nam)
        if fn.func_defaults:
            ndef = len(fn.func_defaults)
        else:
            ndef = 0
        nreq = narg-ndef
        if len(lst)==1:
            if lst[0]==(None,'?'):
                dump_arg(name,arg_nam,nreq)         
                raise QuietException

        if mode==NO_CHECK:
            # no error checking
            for a in lst:
                if a[0]==None:
                    arg.append(a[1])
                else:
                    kw[a[0]]=a[1]
            # set feedback argument (quiet), if extant, results enabled, and not overridden
            if "quiet" in arg_nam:
                if not kw.has_key("quiet"):
                    if __name__!='__main__':
                        if cmd._feedback(cmd.fb_module.cmd,cmd.fb_mask.results):
                            kw["quiet"] = 0
        else:
            # error checking enabled

            # build name dictionary, with required flag
            arg_dct={}
            c = 0
            for a in arg_nam:
                arg_dct[a]=c<nreq
                c = c + 1
            if mode==LEGACY:
                # handle legacy string=value transformation
                tmp_lst = []
                for a in lst:
                    if(a[0]!=None):
                        if not arg_dct.has_key(a[0]):
                            tmp_lst.extend([(None,a[0]),(None,a[1])])
                        else:
                            tmp_lst.append(a)
                    else:
                        tmp_lst.append(a)
                lst = tmp_lst
            # make sure we don't have too many arguments
            if len(lst)>narg:
                if not narg:
                    print "Error: too many arguments for %s; None expected."%(name)
                elif narg==nreq:
                    print "Error: too many arguments for %s; %d expected, %d found."%(
                        name,nreq,len(lst))
                    dump_arg(name,arg_nam,nreq)
                else:
                    print "Error: too many arguments for %s; %d to %d expected, %d found."%(
                        name,nreq,narg,len(lst))
                    dump_arg(name,arg_nam,nreq)            
                raise QuietException
            # match names to unnamed arguments to create argument dictionary
            ac = 0
            val_dct = {}
            for a in lst:
                if a[0]==None:
                    if ac>=narg:
                        print "Parsing-Error: ambiguous argument: '"+str(a[1])+"'"
                        raise QuietException
                    else:
                        val_dct[arg_nam[ac]]=a[1]
                else:
                    val_dct[a[0]]=a[1]
                ac = ac + 1
            # now check to make sure we don't have any missing arguments
            for a in arg_nam:
                if arg_dct[a]:
                    if not val_dct.has_key(a):
                        print "Parsing-Error: missing required argument:",a
                        raise QuietException
            # return all arguments as keyword arguments
            kw = val_dct
            # set feedback argument (quiet), if extant, results enabled, and not overridden
            if arg_dct.has_key("quiet"):
                if not kw.has_key("quiet"):
                    if cmd._feedback(cmd.fb_module.cmd,cmd.fb_mask.results):
                        kw["quiet"] = 0
        if __name__!='__main__':
            if cmd._feedback(cmd.fb_module.parser,cmd.fb_mask.debugging):
                cmd.fb_debug.write(" parsing-DEBUG: kw: "+str(kw)+"\n")      
        return (arg,kw)


    # launching routines

    def run_file(file,global_ns,local_ns):
        pymol.__script__ = file
        try:
            execfile(file,global_ns,local_ns)
        except pymol.CmdException:
            # so the idea here is to print the traceback here and then
            # cascade all the way back up to the interactive level
            # without any further output
            traceback.print_exc()
            raise QuietException
    
    def run_file_as_module(file,spawn=0):
        name = re.sub('[^A-Za-z0-9]','_',file)
        mod = new.module(name)
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
                traceback.print_exc()
                raise QuietException
            del sys.modules[name]
            del mod

    def spawn_file(args,global_ns,local_ns):
        local_ns['__script__'] = args
        t = threading.Thread(target=execfile,args=(args,global_ns,local_ns))
        t.setDaemon(1)
        t.start()

    def split(*arg,**kw): # custom split-and-trim
        '''
    split(string,token[,count]) -> list of strings

    UTILITY FUNCTION, NOT PART OF THE API
    Breaks strings up by tokens but preserves quoted strings and
    parenthetical groups (such as atom selections).

    USAGE OF THIS FUNCTION IS DISCOURAGED - THE GOAL IS TO
    MAKE IT UNNECESSARY BY IMPROVING THE BUILT-IN PARSER
    '''
        str = arg[0]
        tok = arg[1]
        if len(arg)>2:
            mx=arg[2]
        else:
            mx=0
        pair = { '(':')','[':']','{':'}',"'":"'",'"':'"' }
        plst = pair.keys()
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
                lst.append(string.strip(wd))
                nf = nf + 1
                if mx:
                    if nf==mx:
                        wd = string.strip(str[c+1:])
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
            lst.append(string.strip(wd))
        return lst

    import pymol
    
    if __name__=='__main__':

    # regular expressions

        mo = arg_name_re.match("earth=testing,hello")
        tv = mo.group(0)
        print tv == "earth=",tv

        mo = arg_value_re.match("testing,hello")
        tv = mo.group(0)
        print tv == "testing",tv

        mo = arg_value_re.match("'testing,\"hello'")
        tv = mo.group(0)
        print tv == "'testing,\"hello'",tv

        mo = arg_value_re.match('"testing,\'hello"')
        tv = mo.group(0)
        print tv == '"testing,\'hello"',tv

        mo = arg_value_re.match("\'\'\'testing,h\"ello\'\'\'")
        tv = mo.group(0)
        print tv=="'''testing,h\"ello'''",tv



    # argument parsing

        tv = parse_arg("command val")
        print tv==[(None, 'val')],tv
        tv = parse_arg("command val1,val2")
        print tv==[(None, 'val1'), (None, 'val2')],tv
        tv=parse_arg("command val1,val2,val3")
        print tv== [(None, 'val1'), (None, 'val2'), (None, 'val3')],tv      

        tv = parse_arg("command arg=val")
        print tv == [('arg', 'val')],tv
        tv = parse_arg("command arg1=val1,arg2=val2")
        print tv == [('arg1', 'val1'), ('arg2', 'val2')],tv
        tv = parse_arg("command val1,val2,val3")
        print tv == [(None, 'val1'), (None, 'val2'), (None, 'val3')],tv

        tv = parse_arg("command arg=val")
        print tv == [('arg', 'val')] ,tv
        tv = parse_arg("command arg_1=val1,arg2=val2")
        print tv == [('arg_1', 'val1'), ('arg2', 'val2')],tv
        tv = parse_arg("command val1,val2,val3")
        print tv == [(None, 'val1'), (None, 'val2'), (None, 'val3')],tv

        tv = parse_arg("command val1, str2=val2, str3=val3")
        print tv == [(None, 'val1'), ('str2', 'val2'), ('str3', 'val3')],tv
        tv = parse_arg("command val1, str2 = val2,str3= val3")
        print tv == [(None, 'val1'), ('str2', 'val2'), ('str3', 'val3')],tv
        tv = parse_arg("command val1, str2 =val2 ,str3 = val3")
        print tv == [(None, 'val1'), ('str2', 'val2'), ('str3', 'val3')],tv
        tv = parse_arg("command val1, str2 =val2 , str3= val3   ")
        print tv == [(None, 'val1'), ('str2', 'val2'), ('str3', 'val3')],tv

        tv = parse_arg("command multi word 1, str2 =multi word 2")
        print tv == [(None, 'multi word 1'), ('str2', 'multi word 2')],tv

        tv = parse_arg("command multi word 1  , str2 =  multi word 2 ")
        print tv == [(None, 'multi word 1'), ('str2', 'multi word 2')],tv

        tv = parse_arg("command ( name;ca,c,n ), sel1= (name c,n) ")
        print tv == [(None, '( name;ca,c,n )'), ('sel1', '(name c,n)')],tv

        tv = parse_arg("command ( byres (name;ca,c,n) ), sel1= (name c,n) ")
        print tv==[(None, '( byres (name;ca,c,n) )'), ('sel1', '(name c,n)')],tv

        tv = parse_arg("command ( byres (name;ca,c,n ), sel1= (name c,n)) ")
        print tv==[(None, '( byres (name;ca,c,n ), sel1= (name c,n))')],tv

        tv = parse_arg("command test,")
        print tv==[(None,'test')],tv

        tv = parse_arg("command =sdf,") # desired behavior?
        print tv==[(None,'=sdf')],tv

        tv = parse_arg("command 'hello\",bob','''hel\"lo,dad'''") 
        print tv==[(None, '\'hello",bob\''), (None, '\'\'\'hel"lo,dad\'\'\'')],tv

        tv = parse_arg("command \"'hello''bob'\"")
        print tv==[(None, '"\'hello\'\'bob\'"')],tv


        tv = parse_arg("command this,is,a command;b command",mode=LITERAL)
        print tv==[(None, 'this,is,a command;b command')], tv

        tv = parse_arg("command this,a command;b command",mode=LITERAL1)
        print tv==[(None, 'this'), (None, 'a command;b command')], tv

        tv = parse_arg("command this,is,a command;b command",mode=LITERAL2)
        print tv==[(None, 'this'), (None, 'is'), (None, 'a command;b command')],tv

        tv = parse_arg("command this,is,a=hello;b command",mode=LITERAL2)
        print tv==[(None, 'this'), (None, 'is'), (None, 'a=hello;b command')],tv

    # expected exceptions

        try:
            tv = parse_arg("command ( byres (name;ca,c,n ), sel1= (name c,n) ")      
            print 0, "exception missed"
        except QuietException:
            print 1, "exception raised (as expected)"

        try:
            tv = parse_arg("command ,")
            print 0, "exception missed"
        except QuietException:
            print 1, "exception raised (as expected)"

        try:
            tv = parse_arg("command berf=,")
            print 0, "exception missed"
        except QuietException:
            print 1, "exception raised (as expected)"

        try:
            tv = parse_arg("command 'hello''bob'")
            print 0, "exception missed"
        except QuietException:
            print 1, "exception raised"

    # function call preparation

        def fn1(req1,req2,opt1=1,opt2=2):
            pass

        tv = prepare_call(fn1,parse_arg("dummy hello,world"))
        print tv==([], {'req2': 'world', 'req1': 'hello'}),tv

        tv = prepare_call(fn1,parse_arg("dummy hello,world,opt2=hi"))
        print tv==([], {'opt2': 'hi', 'req2': 'world', 'req1': 'hello'}),tv

        tv = prepare_call(fn1,parse_arg("dummy req1=hello,req2=world,opt2=hi")) 
        print tv==([], {'opt2': 'hi', 'req2': 'world', 'req1': 'hello'}),tv

        tv = prepare_call(fn1,parse_arg("dummy req2=world,req1=hello,opt2=hi")) 
        print tv==([], {'opt2': 'hi', 'req2': 'world', 'req1': 'hello'}),tv

        tv = prepare_call(fn1,parse_arg("dummy hello,world,give,tea"))
        print tv==([], {'opt2': 'tea', 'req2': 'world', 'req1': 'hello', 'opt1': 'give'}),tv

        tv = prepare_call(fn1,parse_arg("dummy hello,world,give,tea"),10)
        print tv==(['hello', 'world', 'give', 'tea'], {}),tv

        tv = prepare_call(fn1,parse_arg("dummy hello=world,give,tea"),12) # test legacy support
        print tv==([], {'opt2': 'tea', 'req2': 'world', 'req1': 'hello', 'opt1': 'give'}),tv

        def fn2(*arg):
            pass

        tv = prepare_call(fn2,parse_arg("dummy req1=hello,req2=world,opt2=hi")) 
        print tv==([], {'opt2': 'hi', 'req2': 'world', 'req1': 'hello'}),tv

        def fn3(*arg,**kw):
            pass

        tv = prepare_call(fn3,parse_arg("dummy req1=hello,req2=world,opt2=hi")) 
        print tv==([], {'opt2': 'hi', 'req2': 'world', 'req1': 'hello'}),tv

        def fn4(req1,req2,**kw):
            pass

        tv = prepare_call(fn4,parse_arg("dummy req1=hello,req2=world,opt2=hi")) 
        print tv==([], {'opt2': 'hi', 'req2': 'world', 'req1': 'hello'}),tv

        tv = prepare_call(fn4,parse_arg("dummy req1=hi and (r; 10) and r. 10,req2=world,opt2=hi")) 
        print tv==([], {'opt2': 'hi', 'req2': 'world', 'req1': 'hi and (r; 10) and r. 10'}),tv

        tv = prepare_call(fn4,parse_arg("dummy req1=hi and (r. 10,12) and r. 10,req2=world,opt2=hi")) 
        print tv==([], {'opt2': 'hi', 'req2': 'world', 'req1': 'hi and (r. 10,12) and r. 10'}),tv

        tv = prepare_call(fn4,parse_arg("dummy req1=hi & (r. 10,12) & (r;10),req2=world,opt2=hi")) 
        print tv==([], {'opt2': 'hi', 'req2': 'world', 'req1': 'hi & (r. 10,12) & (r;10)'}),tv

        tv = prepare_call(fn4,parse_arg("dummy req1=hi & (r. 10,12) & (r;10) & resi 10,req2=world,opt2=hi")) 
        print tv==([], {'opt2': 'hi', 'req2': 'world', 'req1': 'hi & (r. 10,12) & (r;10) & resi 10'}),tv

        tv = prepare_call(fn4,parse_arg("dummy req1=(r. 10,12) & (r;10) & resi 10,req2=world,opt2=hi")) 
        print tv==([], {'opt2': 'hi', 'req2': 'world', 'req1': '(r. 10,12) & (r;10) & resi 10'}),tv

    #   tv = list_to_str_list(['hello','world','this-long-string','hi','dude'])
    #   print tv
    else:
        import cmd

        
