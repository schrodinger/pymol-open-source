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

pymol_names = pymol.__dict__

nest=0

com0 = {}
com1 = {}
com2 = {}
cont = {}
script = {}
kw = {}
input = {}
next = {}
args = {}

com0[nest]=""
cont[nest]=""

def parse(s):
   global com0,com1,com2,cont,script,kw,input,next,nest,args,cmd,cont
   local_names = {}
   com0[nest] = s
   try:
      com1[nest] = string.rstrip(com0[nest])
      if len(com1[nest]) > 0:
         if com1[nest][-1] == '\\':
            cont[nest] = cont[nest] + com1[nest][:-1]
         else:
            if cont[nest] != '':
               com1[nest] = cont[nest] + com1[nest]
               cont[nest] = ''
            next[nest] = cmd._split(com1[nest],';',1)
            com2[nest] = next[nest][0]
            input[nest] = string.split(com2[nest],' ',1)
            if len(input[nest]):
               input[nest][0] = string.strip(input[nest][0])
               com = input[nest][0]
               if com[0:1]=='/': # explicit literal python 
                  com2[nest] = string.strip(com2[nest][1:])
                  if len(com2[nest])>0:
                      exec(com2[nest],pymol_names,pymol_names)
               else:
                  if cmd.kwhash.has_key(com):
                     com = cmd.kwhash[com]
                     if not com:
                        print 'Error: ambiguous command: ',
                        com = input[nest][0]
                        lcm = len(com)
                        for a in cmd.keyword.keys():
                           if a[0:lcm] == com:
                              print a+' ',
                        raise RuntimeError
                  if cmd.keyword.has_key(com):
                     kw[nest] = cmd.keyword[com]
                     if kw[nest][4]==1:
                        next[nest] = ()
                        input[nest] = string.split(com1[nest],' ',1)   
                     if len(input[nest])>1:
                        args[nest] = cmd._split(input[nest][1],kw[nest][3])
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
                        if kw[nest][4]<2:
                           result=apply(kw[nest][0],args[nest])
                        elif kw[nest][4]==3:
                           thread.start_new_thread(execfile,
                                                   (args[nest][0],
                                                    pymol_names,local_names))
                        elif len(args[nest])==1:
                           execfile(args[nest][0],pymol_names,local_names)
                        elif args[nest][1]=='local':
                           execfile(args[nest][0],pymol_names,{})
                        elif args[nest][1]=='global':
                           execfile(args[nest][0],pymol_names,pymol_names)
                     else:
                        print 'Error: invalid arguments for %s command.' % com
                  elif len(input[nest][0]):
                     if input[nest][0][0]=='@':
                        script[nest] = open(input[nest][0][1:],'r')
                        nest=nest+1
                        cont[nest]=''
                        while 1:
                           com0[nest] = script[nest-1].readline()
                           if not com0[nest]: break
                           parse(com0[nest])
                           exec(pymol,pymol_names,pymol_names)
                        nest=nest-1
                        script[nest].close()
                     else: # nothing found, try literal python
                        com2[nest] = string.strip(com2[nest])
                        if len(com2[nest])>0:
                           exec(com2[nest],pymol_names,pymol_names)
            if len(next[nest])>1:
               nest=nest+1
               com0[nest] = next[nest-1][1]
               cont[nest]=''
               parse(com0[nest])
               exec(pymol,pymol_names,pymol_names)
               nest=nest-1
   except RuntimeError,e:
      pass
   except StandardError:
      traceback.print_exc()
