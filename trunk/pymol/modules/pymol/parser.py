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

import pymol
import traceback
import string
import cmd
import exceptions
import new
import re
import parsing
import types
import glob

QuietException = parsing.QuietException
      
pymol_names = pymol.__dict__

# parsing state implemented with dictionaries to enable safe recursion
# to arbitrary depths

com0 = {}    # verbose line, as read in
com1 = {}    # line w/o trailing whitespace
com2 = {}    # non-compound command
cont = {}    # continued characters from previous lines (i.e., before \ )
script = {}  # file handles
kw = {}      # row in the cmd.keyword table for the current command
input = {}   # list of length two - command and unparsed arguments string
next = {}    # characters for next command (i.e., after ; )
args = {}    # parsed non-keyword argument string
kw_args = {} # parser keyword argument string

# The resulting value from a pymol command (if any) is stored in the
# parser.result global variable.  However, script developers will
# geerally want to switch to the Python API for any of this kind of
# stuff.

result = None

# initialize parser

nest=0
com0[nest]=""
cont[nest]=""

# main parser routine

def parse(s):
   global com0,com1,com2,cont,script,kw,input
   global next,nest,args,kw_args,cmd,cont,result
   local_names = {}
   com0[nest] = s
   try:
      com1[nest] = string.rstrip(com0[nest]) # strips trailing whitespace
      if len(com1[nest]) > 0:
#         print com1
#         print com1[nest]
         if str(com1[nest][-1]) == "\\":
            cont[nest] = cont[nest] + com1[nest][:-1]
         else:
            if cont[nest] != '':
               com1[nest] = cont[nest] + com1[nest]
               cont[nest] = ''
# this routine splits up the line first based on semicolon 
            next[nest] = parsing.split(com1[nest],';',1)
# com2[nest] now a full non-compound command            
            com2[nest] = next[nest][0]
            input[nest] = string.split(com2[nest],' ',1)
            if len(input[nest]):
               input[nest][0] = string.strip(input[nest][0])
               com = input[nest][0]
               if com[0:1]=='/':
                  # explicit literal python 
                  com2[nest] = string.strip(com2[nest][1:])
                  if len(com2[nest])>0:
                      exec(com2[nest],pymol_names,pymol_names)
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
                        if kw[nest][4]>=parsing.LITERAL: # treat literally
                           next[nest] = ()
                           com2[nest]=com1[nest]
                        (args[nest],kw_args[nest]) = \
                           parsing.prepare_call(
                              kw[nest][0],
                              parsing.parse_arg(com2[nest],mode=kw[nest][4]),
                              kw[nest][4]) # will raise exception on failure
                        result=apply(kw[nest][0],args[nest],kw_args[nest])
                     else:
# old parsing style, being phased out
                        if kw[nest][4]==parsing.ABORT:
                           return None # SCRIPT ABORT EXIT POINT
                        if kw[nest][4]==parsing.SINGLE: # single line, no breaks
                           next[nest] = ()
                           input[nest] = string.split(com1[nest],' ',1)
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
                              # spawn command
                              if len(args[nest])==1: # default: module
                                 parsing.run_as_module(args[nest][0],spawn=1)
                              elif args[nest][1]=='local':
                                 parsing.run_as_thread(args[nest][0],
                                                       pymol_names,local_names)
                              elif args[nest][1]=='global':
                                 parsing.run_as_thread(args[nest][0],
                                                       pymol_names,pymol_names)
                              elif args[nest][1]=='module':
                                 parsing.run_as_module(args[nest][0],spawn=1)
                           elif kw[nest][4]==parsing.RUN:
                              # run command
                              if len(args[nest])==1: # default: global
                                 execfile(args[nest][0],pymol_names,pymol_names)
                              elif args[nest][1]=='local':
                                 execfile(args[nest][0],pymol_names,local_names)
                              elif args[nest][1]=='global':
                                 execfile(args[nest][0],pymol_names,pymol_names)
                              elif args[nest][1]=='module':
                                 parsing.run_as_module(args[nest][0],spawn=0)
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
                        script[nest] = open(string.strip(com2[nest][1:]),'r')
                        nest=nest+1
                        cont[nest]=''
                        while 1:
                           com0[nest]  = script[nest-1].readline()
                           if not com0[nest]: break
                           inp_cmd = com0[nest]
                           tmp_cmd = string.strip(inp_cmd)
                           if len(tmp_cmd):
                              if tmp_cmd[0] not in ['#','_','/']: # suppress comments, internals, python
                                 print "PyMOL>"+tmp_cmd
                           if(parse(inp_cmd)==None): # RECURSION
                              break # abort command gets us out
                        nest=nest-1
                        script[nest].close()
                     else: # nothing found, try literal python
                        com2[nest] = string.strip(com2[nest])
                        if len(com2[nest])>0:
                           exec(com2[nest],pymol_names,pymol_names)
            if len(next[nest])>1:
               nest=nest+1
               com0[nest] = next[nest-1][1]
               next[nest-1]=()
               cont[nest]=''
               parse(com0[nest]) # RECURSION
               nest=nest-1
   except QuietException:
      pass
   except:
      traceback.print_exc()
   return 1

# command and filename completion

def _same_(a,b):
   if a==b:
      return a
   else:
      return None
   
def complete(st):
   result = None
   if string.find(st,' ')<0:
      amb = cmd.kwhash.interpret(st)
      if amb==None:
         print " parser: no matching commands."
      elif type(amb)==types.StringType:
         result = amb+' '
      else:
         amb.sort()
         print " parser: possibilities:"
         flist = filter(lambda x:x[0]!='_',amb)
         lst = parsing.list_to_str_list(flist)
         for a in lst:
            print a
   else: # filename completion
      lst = map(None,st)
      lst.reverse()
      st2 = string.join(lst,'')
      mo = re.search("^[^,\[\( ]+",st2,1)
      if mo!=None:
         pat = mo.group(0)
         lst = map(None,pat)
         lst.reverse()
         st3 = string.join(lst,'')
         loc = string.find(st,st3)
      else:
         st3 = ''
         loc = len(st)
      flist = glob.glob(st3+"*")
      lf = len(flist)
      if lf == 0:
         print " parser: no matching files."
      elif lf==1:
         result = st[0:loc]+flist[0]
      else:
         flist.sort()
         print " parser: possibilities:"
         lst = parsing.list_to_str_list(flist)
         css = map(None,flist[0])
         for a in lst:
            print a
         for a in flist:
            ac = map(None,a)
            css = map(_same_,css,ac)
         css = filter(None,css)
         css = string.join(css,'')
         if len(css)>len(st3):
            result = st[0:loc]+css
   return result
