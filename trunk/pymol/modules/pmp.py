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

# pmp.py 
# Python parser module for PyMol
#
# This compiled code is run in the global namespace...
#
# **Only the PyMOL binary needs to import this file

pmp_nest=0

pmp_cmd = {}
pmp_cmd1 = {}
pmp_cmd2 = {}
pmp_cont = {}
pmp_script = {}
pmp_kw = {}
pmp_input = {}
pmp_next = {}
pmp_args = {}

pmp_cmd[pmp_nest]=""
pmp_cont[pmp_nest]=""

pymol_str="\
try:\n\
   pmp_cmd1[pmp_nest] = string.rstrip(pmp_cmd[pmp_nest])\n\
   if len(pmp_cmd1[pmp_nest]) > 0:\n\
      if pmp_cmd1[pmp_nest][-1] == '\\\\':\n\
         pmp_cont[pmp_nest] = pmp_cont[pmp_nest] + pmp_cmd1[pmp_nest][:-1]\n\
      else:\n\
         if pmp_cont[pmp_nest] != '':\n\
            pmp_cmd1[pmp_nest] = pmp_cont[pmp_nest] + pmp_cmd1[pmp_nest]\n\
            pmp_cont[pmp_nest] = ''\n\
         pmp_next[pmp_nest] = pm._split(pmp_cmd1[pmp_nest],';',1)\n\
         pmp_cmd2[pmp_nest] = pmp_next[pmp_nest][0]\n\
         pmp_input[pmp_nest] = string.split(pmp_cmd2[pmp_nest],' ',1)\n\
         if len(pmp_input[pmp_nest]):\n\
            pmp_input[pmp_nest][0] = string.strip(pmp_input[pmp_nest][0])\n\
            pmp_com = pmp_input[pmp_nest][0]\n\
            if pm.kwhash.has_key(pmp_com):\n\
               pmp_com = pm.kwhash[pmp_com]\n\
               if not pmp_com:\n\
                  print 'Error: ambiguous command.'\n\
                  raise SyntaxError\n\
            if pm.keyword.has_key(pmp_com):\n\
               pmp_kw[pmp_nest] = pm.keyword[pmp_com]\n\
               if pmp_kw[pmp_nest][4]==1:\n\
                  pmp_next[pmp_nest] = ()\n\
                  pmp_input[pmp_nest] = string.split(pmp_cmd1[pmp_nest],' ',1)   \n\
               if len(pmp_input[pmp_nest])>1:\n\
                  pmp_args[pmp_nest] = pm._split(pmp_input[pmp_nest][1],pmp_kw[pmp_nest][3])\n\
                  while 1:\n\
                     pmp_nArg = len(pmp_args[pmp_nest]) - 1\n\
                     pmp_c = 0\n\
                     while pmp_c < pmp_nArg:\n\
                        if (string.count(pmp_args[pmp_nest][pmp_c],'(')!=\n\
                            string.count(pmp_args[pmp_nest][pmp_c],')')):\n\
                           pmp_tmp=pmp_args[pmp_nest][pmp_c+1]\n\
                           pmp_args[pmp_nest].remove(pmp_tmp)\n\
                           pmp_args[pmp_nest][pmp_c]=string.strip(pmp_args[pmp_nest][pmp_c])+\
                                                      ','+string.strip(pmp_tmp)\n\
                           pmp_nArg = pmp_nArg-1\n\
                           break;\n\
                        pmp_c = pmp_c + 1\n\
                     if pmp_c == pmp_nArg:\n\
                        break;\n\
                  if len(pmp_args[pmp_nest])==1 and len(pmp_args[pmp_nest][0])==0:\n\
                     pmp_args[pmp_nest] = []\n\
               else:\n\
                  pmp_args[pmp_nest] = []\n\
               if pmp_kw[pmp_nest][1]<= len(pmp_args[pmp_nest]) <= pmp_kw[pmp_nest][2]:\n\
                  pmp_args[pmp_nest] = map(string.strip,pmp_args[pmp_nest])\n\
                  if pmp_kw[pmp_nest][4]<2:\n\
                     result=apply(pmp_kw[pmp_nest][0],pmp_args[pmp_nest])\n\
                  elif pmp_kw[pmp_nest][4]==3:\n\
                     thread.start_new_thread(execfile,(pmp_args[pmp_nest][0],globals(),locals()))\n\
                  elif len(pmp_args[pmp_nest])==1:\n\
                     execfile(pmp_args[pmp_nest][0],globals(),locals())\n\
                  elif pmp_args[pmp_nest][1]=='local':\n\
                     execfile(pmp_args[pmp_nest][0],globals(),{})\n\
                  elif pmp_args[pmp_nest][1]=='global':\n\
                     execfile(pmp_args[pmp_nest][0],globals(),locals())\n\
               else:\n\
                  print 'Error: invalid arguments for %s command.' % pmp_com\n\
            elif len(pmp_input[pmp_nest][0]):\n\
               if pmp_input[pmp_nest][0][0]=='@':\n\
                  pmp_script[pmp_nest] = open(pmp_input[pmp_nest][0][1:],'r')\n\
                  pmp_nest=pmp_nest+1\n\
                  pmp_cont[pmp_nest]=''\n\
                  while 1:\n\
                     pmp_cmd[pmp_nest] = pmp_script[pmp_nest-1].readline()\n\
                     if not pmp_cmd[pmp_nest]: break\n\
                     exec(pymol,globals(),globals())\n\
                  pmp_nest=pmp_nest-1\n\
                  pmp_script[pmp_nest].close()\n\
               else:   \n\
                  pmp_cmd2[pmp_nest] = string.strip(pmp_cmd2[pmp_nest])\n\
                  if len(pmp_cmd2[pmp_nest])>0:\n\
                     exec(pmp_cmd2[pmp_nest],globals(),globals())\n\
         if len(pmp_next[pmp_nest])>1:\n\
            pmp_nest=pmp_nest+1\n\
            pmp_cmd[pmp_nest] = pmp_next[pmp_nest-1][1]\n\
            pmp_cont[pmp_nest]=''\n\
            exec(pymol,globals(),globals())\n\
            pmp_nest=pmp_nest-1\n\
except SyntaxError:\n\
   pass\n\
except NameError,e:\n\
   print 'NameError: unrecognized command:',e\n"
      
   
pymol = compile(pymol_str,'PyMOL','exec')

