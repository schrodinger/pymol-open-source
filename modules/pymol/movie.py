#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information. 
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-* Peter Haebel
#-* 
#-*
#Z* -------------------------------------------------------------------

import cmd
import glob
import math

def load(*args):
   nam = "mov"
   if len(args)>1:
      nam = args[1]
   fils = glob.glob(args[0])
   fils.sort()
   if not len(fils):
      print "Error: no matching files"
   else:
      for a in fils:
         cmd.load(a,nam)

def rock(first,last,angle=30,phase=0,loop=1,axis='y'):
   first=int(first)
   last=int(last)
   angle=float(angle)
   phase=float(phase)
   loop=int(loop)
   nstep = (last-first)+1
   if nstep<0:
      nstep = 1
   if loop:
      subdiv = nstep
   else:
      subdiv = nstep+1
   ang_cur = math.pi*phase/180
   ang_inc = 2*math.pi/subdiv
   ang_cur = ang_cur - ang_inc
   a = 0
   while a<nstep:
      last = angle*math.sin(ang_cur)/2
      ang_cur = ang_cur + ang_inc
      disp = angle*math.sin(ang_cur)/2
      diff = disp-last
      # com = "mdo %d:turn %s,%8.3f" % (first+a,axis,diff)
      # cmd.do(com)
      cmd.mdo("%d"%(first+a),"turn %s,%8.3f"% (axis,diff))      
      a = a + 1

def roll(first,last,loop=1,axis='y'):
   first=int(first)
   last=int(last)
   loop=int(loop)
   n = last - first
   if loop:
      step = 2*math.pi/(n+1)
   else:
      step = 2*math.pi/n   
   a = 0
   deg = (180*step/math.pi)
   while a<=n:
      # com = "mdo %d:turn %s,%8.3f" % (first+a,axis,deg)
      # cmd.do(com)
      cmd.mdo("%d" % (first+a), "turn %s,%8.3f" % (axis,deg))
      a = a + 1


def zoom(first,last,step=1,loop=1,axis='z'):
   # Author: Peter Haebel
   first=int(first)
   last=int(last)
   step=int(step)
   loop=int(loop)
   n = last - first
   a = 0
   while a<=n:
      if (loop and a>n/2):
         s = -step
      else:
         s = step
      # com = "mdo %d:move %s,%8.3f" % (first+a,axis,s)
      # cmd.do(com)
      cmd.mdo("%d" % (first+a),"move %s,%8.3f" % (axis,s))
      a = a + 1

def nutate(first,last,angle=30,phase=0,loop=1,shift=math.pi/2.0,factor=0.01):
   first=int(first)
   last=int(last)
   angle=float(angle)
   phase=float(phase)
   loop=int(loop)
   nstep = (last-first)+1
   if nstep<0:
      nstep = 1
   if loop:
      subdiv = nstep
   else:
      subdiv = nstep+1
   ang_cur = math.pi*phase/180
   ang_inc = 2*math.pi/subdiv
   ang_cur = ang_cur - ang_inc
   a = 0
   while a<nstep:
      lastx = angle*math.sin(ang_cur)/2
      lasty = angle*math.sin(ang_cur+shift)/2
      ang_cur = ang_cur + ang_inc
      nextx = angle*math.sin(ang_cur)/2
      nexty = angle*math.sin(ang_cur+shift)/2      
      # com = "mdo %d:turn %s,%8.3f" % (first+a,axis,diff)
      # cmd.do(com)
      cmd.mdo("%d"%(first+a),"turn x,%8.3f;turn y,%8.3f;turn y,%8.3f;turn x,%8.3f"%
              (-lastx,-lasty,nexty,nextx))
      a = a + 1

def screw(first,last,step=1,angle=30,phase=0,loop=1,axis='y'):
   # Author: Peter Haebel
   first=int(first)
   last=int(last)
   step=int(step)
   angle=float(angle)
   phase=float(phase)
   loop=int(loop)
   nstep = (last-first)+1
   if nstep<0:
      nstep = 1
   if loop:
      subdiv = nstep
   else:
      subdiv = nstep+1
   ang_cur = math.pi*phase/180
   ang_inc = 2*math.pi/subdiv
   ang_cur = ang_cur - ang_inc
   a = 0
   while a<nstep:
      if (loop and a>=nstep/2):
         s = -step
      else:
         s = step
      last = angle*math.sin(ang_cur)/2
      ang_cur = ang_cur + ang_inc
      disp = angle*math.sin(ang_cur)/2
      diff = disp-last
      # com = "mdo %d:turn %s,%8.3f; move z,%8.3f" % (first+a,axis,diff,s)
      # cmd.do(com)
      cmd.mdo("%d" % (first+a), "turn %s,%8.3f; move z,%8.3f" % (axis,diff,s))
      a = a + 1

