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
import pymol

def cbag(selection):
   s = str(selection)
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("carbon","(elem C and "+s+")")

def cbac(selection):
   s = str(selection)
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("cyan","(elem C and "+s+")")

def cbay(selection):
   s = str(selection)   
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("yellow","(elem C and "+s+")")

def cbas(selection):
   s = str(selection)   
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("salmon","(elem C and "+s+")")

def cbap(selection):
   s = str(selection)   
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("purple","(elem C and "+s+")")

def cbaw(selection):
   s = str(selection)   
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("hydrogen","(elem C and "+s+")")

def cbab(selection):
   s = str(selection)   
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("slate","(elem C and "+s+")")

def mrock(first,last,angle,phase,loop):
   fir=int(first)
   las=int(last)
   dsp=float(angle)
   pha=float(phase)
   loop=int(loop)
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
   fir=int(first)
   las=int(last)
   loop=int(loop)
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
   
def cbc(selection='(all)'):
   '''
   Color all chains a different color
   '''
   pymol.stored.chain = {}
   cmd.iterate("(%s)"%selection,"stored.chain[chain]=1")
   c = 7
   for a in pymol.stored.chain.keys():
      if len(a):
         print ("%d,(chain %s)"%(c,a))
         cmd.color("%d"%c,"(chain %s)"%a)
         c = c + 1

color_chains = cbc

def sum_charge(*arg):
   result = None
   try:
      obj = "all"
      if len(arg):
         obj = arg

      pymol.stored._sum_charge = 0.0
      cmd.iterate("(%s)"%obj,
                  "stored._sum_charge=stored._sum_charge+partial_charge")
      result = pymol.stored._sum_charge
      print " sum_charge: %6.4f"%result
   except:
      print " sum_charge: an error occurred."
   return result


def ff_copy(src,dst):
   pymol._rcopy = pymol.Scratch_Storage()
   pymol._rcopy.pc={}
   pymol._rcopy.tt={}
   cmd.iterate("(%s)"%src,"_rcopy.pc[name]=partial_charge")
   cmd.alter("(%s)"%dst,"partial_charge=_rcopy.pc[name]")
   cmd.iterate("(%s)"%src,"_rcopy.tt[name]=text_type")
   cmd.alter("(%s)"%dst,"text_type=_rcopy.tt[name]")
   del pymol._rcopy
   
def b2vdw(*arg):
   if not len(arg):
      sele = 'all'
   else:
      sele = arg[0]
   # use B values to create RMS VDW spheres
   # rms = sqrt(b/(8*(PI^2)))
   cmd.alter("(%s)"%sele,"vdw=math.sqrt(b/78.9568352087)")
   
   
