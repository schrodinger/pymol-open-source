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

import cmd
import math
import string
import glob


def cbag(s):
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("carbon","(elem C and "+s+")")

def cbac(s):
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("cyan","(elem C and "+s+")")

def cbay(s):
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("yellow","(elem C and "+s+")")

def cbas(s):
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("salmon","(elem C and "+s+")")

def cbap(s):
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("purple","(elem C and "+s+")")

def cbaw(s):
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("hydrogen","(elem C and "+s+")")


def mrock(fir,las,dsp,pha,loop):
   n = las - fir
   ang = pha * math.pi
   if loop:
      step = 2*math.pi/(n+1)
      last = -(math.sin(ang+step*n)*dsp);
   else:
      last=0
      step = 2*math.pi/n   
   a = 0
   while a<=n:
      deg = (math.sin(ang)*dsp)
      com = "mdo %d:turn y,%8.3f;turn y,%8.3f" % (fir+a,last,deg)
      last = -deg;
      print com
      cmd.do(com)
      ang = ang + step
      a = a + 1

def mroll(fir,las,loop):
   n = las - fir
   if loop:
      step = 2*math.pi/(n+1)
   else:
      step = 2*math.pi/n   
   a = 0
   deg = (180*step/math.pi)
   while a<=n:
      com = "mdo %d:turn y,%8.3f" % (fir+a,deg)
      print com
      cmd.do(com)
      a = a + 1

def hbond(a,b,cutoff=3.3):
   st = "(%s and (%s around %4.2f) and elem N,O),(%s and (%s around %4.2f) and elem N,O),%4.2f" % (a,b,cutoff,b,a,cutoff,cutoff)
   cmd.dist("hbond",st)
        
def mload(*args):
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
   
