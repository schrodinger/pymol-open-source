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

# parser.py
# Python parser module for PyMol
#

if __name__=='pymol.parser':
    
    import pymol
    import traceback
    import string
    import cmd
    import exceptions
    import re
    import parsing
    import types
    import glob
    import sys
    import os
    import setting
    import __main__

    from cmd import _feedback,fb_module,fb_mask,exp_path

    QuietException = parsing.QuietException
    CmdException = pymol.CmdException
    
    pymol_names = pymol.__dict__

    py_delims = { '=' : 1, '+='  : 1, '-='  : 1, '*=' : 1,
                      '/=' :1, '//=' : 1, '%='  : 1, '&=' : 1,
                      '|=' :1, '^='  : 1, '>>=' : 1,'<<=' : 1,
                      '**=':1 }

    # parsing state implemented with dictionaries to enable safe recursion
    # to arbitrary depths

    com0 = {}    # verbose line, as read in
    com1 = {}    # line w/o trailing whitespace
    com2 = {}    # non-compound command
    cont = {}    # continued characters from previous lines (i.e., before \ )
    script = {}  # file handles
    sc_path = {}   # file paths
    kw = {}      # row in the cmd.keyword table for the current command
    input = {}   # list of length two - command and unparsed arguments string
    next = {}    # characters for next command (i.e., after ; )
    args = {}    # parsed non-keyword argument string
    kw_args = {} # parser keyword argument string
    embed_dict = {}
    embed_list = {}
    embed_sentinel = {}
    embed_type = {}
    embed_line = {}
    
    # The resulting value from a pymol command (if any) is stored in the
    # parser.result global variable.  However, script developers will
    # geerally want to switch to the Python API for any of this kind of
    # stuff.

    result = None

    # initialize parser

    nest=0
    com0[nest]=""
    cont[nest]=""
    sc_path[nest]="default"
    embed_sentinel[nest]=None

    def stdin_reader(): # dedicated thread for reading standard input
        global nest, embed_sentinel
        import sys
        from pymol import cmd
        while 1:
            l = sys.stdin.readline()
            if l!="":
                if nest==0:
                    # if we're reading embedded input on stdin
                    # then bypass PyMOL C code altogether
                    if embed_sentinel[nest]!=None:
                        parse(l)
                    else:
                        cmd.do(l)
                else:
                    cmd.do(l)
            elif not pymol.invocation.options.keep_thread_alive:
                cmd.quit()
        pymol._stdin_reader_thread = None
        
    def get_default_key():
        return os.path.splitext(os.path.basename(sc_path[nest]))[0]

    def get_embedded(key=None):
        if key==None:
            key = get_default_key()
        dict = embed_dict.get(nest,{})
        return dict.get(key,None)
        
    # main parser routine

    def parse(s,secure=0):
        global nest,cmd,result
    # report any uncaught errors...
        if sys.exc_info()!=(None,None,None):
            traceback.print_exc()
            sys.exc_clear()
        if embed_sentinel[nest]!=None:
            if string.strip(s) == embed_sentinel[nest]:
                etn = embed_type[nest]
                if etn == 0: # embedded data
                    print " Embed: read %d lines."%(len(embed_list[nest]))
                    embed_sentinel[nest]=None
                elif etn == 1: # python block
                    print "PyMOL>"+string.rstrip(s)                    
                    py_block = string.join(embed_list[nest],'')
                    del embed_list[nest]
                    embed_sentinel[nest]=None
                    exec(py_block,pymol_names,pymol_names)
                elif etn == 2: # skip block
                    print " Skip: skipped %d lines."%(embed_line[nest])
                    embed_sentinel[nest]=None
            else:
                etn = embed_type[nest]
                if etn == 0: # normal embedded data
                    embed_list[nest].append(string.rstrip(s)+"\n")
                elif etn == 1: # python block
                    el = embed_line[nest] + 1
                    print "%5d:%s"%(el,string.rstrip(s))
                    embed_line[nest] = el
                    embed_list[nest].append(string.rstrip(s)+"\n")
                elif etn == 2:
                    embed_line[nest] = embed_line[nest] + 1
            return 1
        p_result = 1
        com0[nest] = s
        try:
            com1[nest] = string.rstrip(com0[nest]) # strips trailing whitespace
            if len(com1[nest]) > 0:
                if str(com1[nest][-1]) == "\\":
                    # prepend leftovers
                    if cont[nest] != '':
                        cont[nest] = cont[nest] + "\n" + com1[nest][:-1]
                    else:
                        cont[nest] = com1[nest][:-1]
                else:
                    # prepend leftovers
                    if cont[nest] != '':
                        com1[nest] = cont[nest] + "\n" + com1[nest]
                        cont[nest] = ''
    # this routine splits up the line first based on semicolon 
                    next[nest] = parsing.split(com1[nest],';',1)
    # com2[nest] now a full non-compound command            
                    com2[nest] = next[nest][0]
                    input[nest] = string.split(com2[nest],' ',1)
                    lin = len(input[nest])
                    if lin:
                        input[nest][0] = string.strip(input[nest][0])
                        com = input[nest][0]
                        if (com[0:1]=='/'):
                            # explicit literal python 
                            com2[nest] = string.strip(com2[nest][1:])
                            if len(com2[nest])>0:
                                if not secure:
                                    exec(com2[nest]+"\n",pymol_names,pymol_names)
                                else:
                                    print 'Error: Python expressions disallowed in this file.'
                                    return None
                        elif lin>1 and py_delims.has_key(string.split(input[nest][-1:][0],' ',1)[0]):
                            if not secure:
                                exec(com2[nest]+"\n",pymol_names,pymol_names)
                            else:
                                print 'Error: Python expressions disallowed in this file.'
                                return None
                        else:
                            # try to find a keyword which matches
                            if cmd.kwhash.has_key(com):
                                amb = cmd.kwhash.interpret(com)
                                if amb == None:
                                    com = cmd.kwhash[com]
                                elif type(amb)!=types.StringType:
                                    print 'Error: ambiguous command: '
                                    amb.sort()
                                    amb = parsing.list_to_str_list(amb)
                                    for a in amb:
                                        print a
                                    raise QuietException
                                com = amb
                            if cmd.keyword.has_key(com):
    # here is the command and argument handling section
                                kw[nest] = cmd.keyword[com]
                                if kw[nest][4]>=parsing.NO_CHECK:
    # stricter, Python-based argument parsing
                                    # remove line breaks (only important for Python expressions)
                                    com2[nest]=string.replace(com2[nest],'\n','')

                                    if kw[nest][4]>=parsing.LITERAL: # treat literally
                                        next[nest] = ()
                                        if not secure:
                                            com2[nest]=com1[nest]
                                        else:
                                            print 'Error: Python expressions disallowed in this file.  '
                                            return 0
                                    if secure and (kw[nest][4]==parsing.SECURE):
                                        next[nest] = ()
                                        print 'Error: Command disallowed in this file.'
                                        return None
                                    else:
                                       (args[nest],kw_args[nest]) = \
                                                                        parsing.prepare_call(
                                        kw[nest][0],
                                        parsing.parse_arg(com2[nest],mode=kw[nest][4]),
                                        kw[nest][4]) # will raise exception on failure
                                    result=apply(kw[nest][0],args[nest],kw_args[nest])
                                elif kw[nest][4]==parsing.PYTHON:
                                        # handle python keyword
                                        com2[nest] = string.strip(com2[nest])
                                        if len(com2[nest])>0:
                                            if not secure:
                                                exec(com2[nest]+"\n",pymol_names,pymol_names)
                                            else:
                                                next[nest] = ()                                    
                                                print 'Error: Python expressions disallowed in this file.'
                                                return None
                                else:
                                    # remove line breaks (only important for Python expressions)
                                    com2[nest]=string.replace(com2[nest],'\n','')
    # old parsing style, being phased out
                                    if kw[nest][4]==parsing.ABORT:
                                        return None # SCRIPT ABORT EXIT POINT
                                    if kw[nest][4]==parsing.MOVIE: # copy literal single line, no breaks
                                        next[nest] = ()
                                        if not secure:
                                            input[nest] = string.split(com1[nest],' ',1)
                                        else:
                                            print 'Error: Movie commands disallowed in this file. '
                                            return None
                                    if len(input[nest])>1:
                                        args[nest] = parsing.split(input[nest][1],kw[nest][3])
                                        while 1:
                                            nArg = len(args[nest]) - 1
                                            c = 0
                                            while c < nArg:
                                                if (string.count(args[nest][c],'(')!=
                                                     string.count(args[nest][c],')')):
                                                    tmp=args[nest][c+1]
                                                    args[nest].remove(tmp)
                                                    args[nest][c]=string.strip(args[nest][c])+\
                                                                      ','+string.strip(tmp)
                                                    nArg = nArg-1
                                                    break;
                                                c = c + 1
                                            if c == nArg:
                                                break;
                                        if len(args[nest])==1 and len(args[nest][0])==0:
                                            args[nest] = []
                                    else:
                                        args[nest] = []
                                    if kw[nest][1]<= len(args[nest]) <= kw[nest][2]:
                                        args[nest] = map(string.strip,args[nest])
                                        if kw[nest][4]<parsing.RUN:
        #                           
        # this is where old-style commands are invoked
        #
                                            result=apply(kw[nest][0],args[nest])
        #                           
                                        elif kw[nest][4]==parsing.SPAWN:
                                            if not secure:
                                                # spawn command
                                                if len(args[nest])==1: # default: module
                                                    parsing.run_file_as_module(exp_path(args[nest][0]),spawn=1)
                                                elif args[nest][1]=='main':
                                                    parsing.spawn_file(exp_path(args[nest][0]),
                                                                                 __main__.__dict__,
                                                                                 __main__.__dict__)
                                                elif args[nest][1]=='private':
                                                    parsing.spawn_file(exp_path(args[nest][0]),
                                                                                 __main__.__dict__,
                                                                                 {})
                                                elif args[nest][1]=='local':
                                                    parsing.spawn_file(exp_path(args[nest][0]),
                                                                                 pymol_names,{})
                                                elif args[nest][1]=='global':
                                                    parsing.spawn_file(exp_path(args[nest][0]),
                                                                                 pymol_names,pymol_names)
                                                elif args[nest][1]=='module':
                                                    parsing.run_file_as_module(exp_path(args[nest][0]),spawn=1)
                                            else:
                                                next[nest] = ()                                    
                                                print 'Error: spawn disallowed in this file.'
                                                return None
                                        elif kw[nest][4]==parsing.RUN:
                                            if not secure:
                                                # run command
                                                if len(args[nest])==1: # default: global
                                                    parsing.run_file(exp_path(args[nest][0]),pymol_names,pymol_names)
                                                elif args[nest][1]=='main':
                                                    parsing.run_file(exp_path(args[nest][0]),__main__.__dict__,__main__.__dict__)
                                                elif args[nest][1]=='private':
                                                    parsing.run_file(exp_path(args[nest][0]),__main__.__dict__,{})
                                                elif args[nest][1]=='local':
                                                    parsing.run_file(exp_path(args[nest][0]),pymol_names,{})
                                                elif args[nest][1]=='global':
                                                    parsing.run_file(exp_path(args[nest][0]),pymol_names,pymol_names)
                                                elif args[nest][1]=='module':
                                                    parsing.run_file_as_module(exp_path(args[nest][0]),spawn=0)
                                            else:
                                                next[nest] = ()                                    
                                                print 'Error: run disallowed in this file.'
                                                return None                                    
                                        elif (kw[nest][4]==parsing.EMBED):
                                            next[nest] = ()
                                            if secure or nest==0: # only legal on top level and p1m files
                                                l = len(args[nest])
                                                if l>0:
                                                    key = args[nest][0]
                                                else:
                                                    key = os.path.splitext(os.path.basename(sc_path[nest]))[0]
                                                if l>1:
                                                    format = args[nest][1]
                                                else:
                                                    format = 'pdb'
                                                if l>2:
                                                    embed_sentinel[nest] = args[nest][2]
                                                else:
                                                    embed_sentinel[nest] = "embed end"
                                                list = []
                                                dict = embed_dict.get(nest,{})
                                                dict[key] = ( format, list )
                                                embed_dict[nest] = dict
                                                embed_list[nest] = list
                                                embed_type[nest] = 0 # not a python block
                                            else:
                                                print 'Error: embed only legal in p1m files'
                                                raise None
                                        elif (kw[nest][4]==parsing.SKIP):
                                            next[nest] = ()
                                            l = len(args[nest])
                                            if l>0:
                                                embed_sentinel[nest] = args[nest][0]
                                            else:
                                                embed_sentinel[nest] = "skip end"
                                            embed_type[nest] = 2 # skip block
                                            embed_line[nest] = 0                                            
                                        elif (kw[nest][4]==parsing.PYTHON_BLOCK):
                                            next[nest] = ()
                                            if not secure: 
                                                l = len(args[nest])
                                                if l>0:
                                                    embed_sentinel[nest] = args[nest][0]
                                                else:
                                                    embed_sentinel[nest] = "python end"
                                                list = []
                                                embed_list[nest] = list
                                                embed_type[nest] = 1 # python block
                                                embed_line[nest] = 0
                                            else:
                                                print 'Error: Python blocks disallowed in this file.'
                                                raise None
                                        else:
                                            print 'Error: unknown keyword mode: '+str(kw[nest][4])
                                            raise QuietException
                                    else:
                                        print 'Error: invalid arguments for %s command.' % com
    #
    # non-keyword command handling
    #
                            elif len(input[nest][0]):
                                if input[nest][0][0]=='@':
                                    path = exp_path(string.strip(com2[nest][1:]))
                                    if string.lower(path[-3:])=='p1m':
                                        nest_secure = 1
                                    else:
                                        nest_secure = secure
                                    script[nest] = open(path,'r')
                                    nest=nest+1
                                    cont[nest]=''
                                    sc_path[nest]=path
                                    embed_sentinel[nest]=None
                                    while 1:
                                        com0[nest]  = script[nest-1].readline()
                                        if not com0[nest]: break
                                        inp_cmd = com0[nest]
                                        tmp_cmd = string.strip(inp_cmd)
                                        if len(tmp_cmd):
                                            if tmp_cmd[0] not in ['#','_','/']: # suppress comments, internals, python
                                                if embed_sentinel[nest]==None:
                                                    print "PyMOL>"+tmp_cmd
                                            elif tmp_cmd[0]=='_' and \
                                                  tmp_cmd[1:2] in [' ','']: # "_ " remove echo suppression signal
                                                inp_cmd=inp_cmd[2:]
                                        pp_result = parse(inp_cmd,nest_secure)
                                        if pp_result==None: # RECURSION
                                            break # abort command gets us out
                                        elif pp_result==0: # QuietException
                                            if cmd.get_setting_legacy("stop_on_exceptions"):
                                                p_result = 0 # signal an error occurred
                                                print"PyMOL: stopped on exception."
                                                break;
                                    nest=nest-1
                                    script[nest].close()
                                else: # nothing found, try literal python
                                    com2[nest] = string.strip(com2[nest])
                                    if len(com2[nest])>0:
                                        if not secure:
                                            exec(com2[nest]+"\n",pymol_names,pymol_names)
                                        elif input[nest][0][0:1]!='#':
                                            print 'Error: unrecognized keyword: '+input[nest][0]
                    if (len(next[nest])>1) and p_result:
                              # continue parsing if no error or break has occurred
                        nest=nest+1
                        com0[nest] = next[nest-1][1]
                        next[nest-1]=()
                        cont[nest]=''
                        embed_sentinel[nest]=None
                        p_result = parse(com0[nest]) # RECURSION
                        nest=nest-1
        except QuietException:
            if cmd._feedback(fb_module.parser,fb_mask.blather):
                print "Parser: QuietException caught"
            p_result = 0 # notify caller that an error was encountered
        except CmdException:
            if cmd._feedback(fb_module.parser,fb_mask.blather):         
                print "Parser: CmdException caught."
            p_result = 0
        except:
            traceback.print_exc()
            if cmd._feedback(fb_module.parser,fb_mask.blather):
                print "PyMOL: Caught an unknown exception."
            p_result = 0 # notify caller that an error was encountered
        return p_result  # 0 = Exception, None = abort, 1 = ok

    # command and filename completion

    def _same_(a,b):
        if a==b:
            return a
        else:
            return None

    remove_lists_re = re.compile("\[[^\]]*\]")
    
    def complete_sc(st,sc,type_name,postfix):
        result = None
        sc=sc() # invoke lambda functions (if any)
        amb = sc.interpret(st)
        if amb==None:
            print " parser: no matching %s."%type_name
        elif type(amb)==types.StringType:
            result = amb+postfix
        else:
            amb.sort()
            print " parser: matching %s:"%type_name
            flist = filter(lambda x:x[0]!='_',amb)
            lst = parsing.list_to_str_list(flist)
            for a in lst:
                print a
            # now append up to point of ambiguity
            if not len(flist):
                css = []
            else:
                css = map(None,flist[0]) # common sub-string (css)
                for a in flist:
                    ac = map(None,a)
                    tmp = css
                    css = []
                    for c in range(len(tmp)):
                        if tmp[c]!=ac[c]:
                            break
                        css.append(tmp[c])
            css = filter(None,css)
            css = string.join(css,'')
            if len(css)>len(st):
                result = css
        return result

    def complete(st):
        result = None
        pre = ''
        flag = 0
        if (string.find(st,' ')<0) and (string.find(st,'@'))<0:
            try:
                result = complete_sc(st,cmd.kwhash,'commands',' ')
            except:
                traceback.print_exc()
        else:
            full = cmd.kwhash.interpret(re.sub(r" .*","",st))
            st_no_lists = remove_lists_re.sub("",st)
            count = string.count(st_no_lists,',') # which argument are we on
            if cmd.is_string(full):
                if count<len(cmd.auto_arg):
                    if cmd.auto_arg[count].has_key(full): # autocomplete arguments
                        flag = 1
                        try:
                            pre = re.sub(r"^[^ ]* ",' ',st,count=1) # trim command
                            if re.search(r",",pre)!=None:
                                pre = re.sub(r"[^\, ]*$","",pre,count=1) 
                                pre = re.sub(r",\s*[^\, ]*$",", ",pre,count=1) # trim 1 arg
                            else:
                                pre = re.sub("[^ ]*$","",pre,count=1) # trim 1 arg               
                            pre = re.sub(r"^ *",'',pre)
                            pre = full+' '+pre
                            pat = re.sub(r".*[\, ]",'',st)
        #               print ":"+pre+":"+pat+":"
                            result = apply(complete_sc,tuple([pat]+cmd.auto_arg[count][full]),{})
                        except:
                            traceback.print_exc()
            if not flag: # otherwise fallback onto filename completion
                if(st[:1]=='@'):
                    st=st[1:]
                    pre = '@'
                loc = reduce(max,[string.rfind(st,','),
                                        string.rfind(st,' '),
                                        string.rfind(st,']'),
                                        string.rfind(st,')')])+1
                st3 = st[loc:]
                flist = glob.glob(exp_path(st3)+"*")
                lst = map(None,st3)
                lst.reverse()
                st3 = string.join(lst,'')
                lf = len(flist)
                if lf == 0:
                    print " parser: no matching files."
                elif lf==1:
                    result = st[0:loc]+flist[0]
                    if os.path.isdir(flist[0]):
                        result = result + "/"
                else:
                    flist.sort()
                    print " parser: matching files:"
                    lst = parsing.list_to_str_list(flist)
                    for a in lst:
                        print a
                    # now append as much up to point of ambiguity
                    css = map(None,flist[0]) # common sub-string (css)
                    for a in flist:
                        ac = map(None,a)
                        tmp = css
                        css = []
                        for c in range(len(tmp)):
                            if tmp[c]!=ac[c]:
                                break
                            css.append(tmp[c])
                    tmp = css
                    css = []
                    for a in tmp:
                        if a==None:
                            break
                        css.append(a)
                    css = string.join(css,'')
                    if len(css)>len(st3):
                        result = st[0:loc]+css
        if result!=None:
            result = pre+result
        return result
