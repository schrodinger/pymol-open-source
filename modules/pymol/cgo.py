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

import string
from chempy import cpv
#import popen2
import os
from pymol import cmd

POINTS             = 0.0
LINES              = 1.0
LINE_LOOP          = 2.0
LINE_STRIP         = 3.0
TRIANGLES          = 4.0
TRIANGLE_STRIP     = 5.0
TRIANGLE_FAN       = 6.0
QUADS              = 7.0
QUAD_STRIP         = 8.0
POLYGON            = 9.0                                                            

STOP               =  0.0
NULL               =  1.0
BEGIN              =  2.0
END                =  3.0
VERTEX             =  4.0
NORMAL             =  5.0
COLOR              =  6.0
SPHERE             =  7.0
TRIANGLE           =  8.0
CYLINDER           =  9.0
LINEWIDTH          = 10.0
WIDTHSCALE         = 11.0


def molauto(*arg):
   name = "mols"
   sele = "(all)"
   marg = "-nice"
   la = len(arg)
   if la:
      name = arg[0]
   if la>1:
      sele = arg[1]
   if la>2:
      marg = arg[2]
   cmd.save("molauto.pdb",sele)

   os.system("molauto %s -nocentre molauto.pdb | molscript -r > molauto.r3d"%marg)
   f = open("molauto.r3d")
   rr = RenderReader(f)
   f.close()
   cmd.load_cgo(rr.obj,name)

# the following implementation causes full-blown system crashes on some machines.
#   (stdout,stdin) = popen2.popen2("molauto %s -nocentre molauto.pdb | molscript -r > molauto.r3d"%marg)
#
#   if stdin:
#      stdin.close()
#      rr = RenderReader(stdout)
#      cmd.load_cgo(rr.obj,name)

def from_r3d(fname):
   result = None
   input = open(fname)
   if input:
      rr = RenderReader(input)
      result = rr.obj
   return result

class RenderReader:

   def append_last(self):
      if self.app_fn:
         apply(self.app_fn)
         self.app_fn=None
      
   def append_tri(self):
      if self.l_vert and not self.l_norm:
         d0 = cpv.sub(self.l_vert[0],self.l_vert[1])
         d1 = cpv.sub(self.l_vert[0],self.l_vert[2])
         n0 = cpv.cross_product(d0,d1)
         n0 = cpv.normalize_failsafe(n0)
         n1 = [-n0[0],-n0[1],-n0[2]]
         ns = cpv.scale(n0,0.002)
         if not self.tri_flag:
            self.obj.append(BEGIN)
            self.obj.append(TRIANGLES)
            self.tri_flag = 1
         self.obj.append(COLOR)  # assuming unicolor
         self.obj.extend(self.t_colr[0])
         self.obj.append(NORMAL)
         self.obj.extend(n0)
         self.obj.append(VERTEX)
         self.obj.extend(cpv.add(self.l_vert[0],ns))
         self.obj.append(COLOR)  # assuming unicolor
         self.obj.extend(self.t_colr[1])
         self.obj.append(NORMAL)
         self.obj.extend(n0)
         self.obj.append(VERTEX)
         self.obj.extend(cpv.add(self.l_vert[1],ns))
         self.obj.append(COLOR)  # assuming unicolor
         self.obj.extend(self.t_colr[2])
         self.obj.append(NORMAL)
         self.obj.extend(n0)
         self.obj.append(VERTEX)
         self.obj.extend(cpv.add(self.l_vert[2],ns))
         self.obj.append(COLOR)  # assuming unicolor
         self.obj.extend(self.t_colr[0])
         self.obj.append(NORMAL)
         self.obj.extend(n1)
         self.obj.append(VERTEX)
         self.obj.extend(cpv.sub(self.l_vert[0],ns))
         self.obj.append(COLOR)  # assuming unicolor
         self.obj.extend(self.t_colr[1])
         self.obj.append(NORMAL)
         self.obj.extend(n1)
         self.obj.append(VERTEX)
         self.obj.extend(cpv.sub(self.l_vert[1],ns))
         self.obj.append(COLOR)  # assuming unicolor
         self.obj.extend(self.t_colr[2])
         self.obj.append(NORMAL)
         self.obj.extend(n1)
         self.obj.append(VERTEX)
         self.obj.extend(cpv.sub(self.l_vert[2],ns))

      elif self.l_vert and self.t_colr and self.l_norm:
         if not self.tri_flag:
            self.obj.append(BEGIN)
            self.obj.append(TRIANGLES)
            self.tri_flag = 1
         self.obj.append(COLOR) # assuming unicolor
         self.obj.extend(self.t_colr[0])
         self.obj.append(NORMAL)
         self.obj.extend(self.l_norm[0])
         self.obj.append(VERTEX)
         self.obj.extend(self.l_vert[0])
         self.obj.append(COLOR) # assuming unicolor
         self.obj.extend(self.t_colr[1])
         self.obj.append(NORMAL)
         self.obj.extend(self.l_norm[1])
         self.obj.append(VERTEX)
         self.obj.extend(self.l_vert[1])
         self.obj.append(COLOR) # assuming unicolor
         self.obj.extend(self.t_colr[2])
         self.obj.append(NORMAL)
         self.obj.extend(self.l_norm[2])
         self.obj.append(VERTEX)
         self.obj.extend(self.l_vert[2])
      self.l_vert=None
      self.t_colr=None
      self.l_norm=None

   def append_cyl(self):
      if self.l_vert and self.c_colr and self.l_radi:
         if self.tri_flag:
            self.tri_flag=0
            self.obj.append(END)
         self.obj.append(CYLINDER)
         d = cpv.sub(self.l_vert[1],self.l_vert[0])
         d = cpv.normalize_failsafe(d)
         d0 = cpv.scale(d,self.l_radi/4.0)
         self.obj.extend(cpv.add(self.l_vert[0],d0))
         self.obj.extend(cpv.sub(self.l_vert[1],d0))
         self.obj.append(self.l_radi)
         self.obj.extend(self.c_colr[0])
         self.obj.extend(self.c_colr[1])
      self.l_vert=None
      self.c_colr=None
      self.l_radi=None

   def tri(self,f):
      self.append_last()
      l = f.readline()
      if l:
         self.app_fn=self.append_tri
         s = string.split(l)
         self.l_vert = [[float(s[0]),float(s[1]),float(s[2])],
                  [float(s[3]),float(s[4]),float(s[5])],
                  [float(s[6]),float(s[7]),float(s[8])]]
         self.t_colr_t = [float(s[9]),float(s[10]),float(s[11])]
         self.t_colr = [self.t_colr_t,self.t_colr_t,self.t_colr_t]

   def tri_normal(self,f):
      l = f.readline()
      if l:
         s = string.split(l)
         self.l_norm = [[float(s[0]),float(s[1]),float(s[2])],
                  [float(s[3]),float(s[4]),float(s[5])],
                  [float(s[6]),float(s[7]),float(s[8])]]

   def cyl(self,f):
      self.append_last()
      l = f.readline()
      if l:
         self.app_fn = self.append_cyl
         s = string.split(l)
         self.l_vert = [[float(s[0]),float(s[1]),float(s[2])],
                  [float(s[4]),float(s[5]),float(s[6])]]
         self.l_radi = float(s[3])
         self.c_colr_t = [float(s[8]),float(s[9]),float(s[10])]
         self.c_colr = [self.c_colr_t,self.c_colr_t]

   def sphere(self,f):
      self.append_last()
      l = f.readline()
      if l:
         s = string.split(l)
         self.obj.append(COLOR)
         self.obj.extend([float(s[4]),float(s[5]),float(s[6])])
         self.obj.append(SPHERE)
         self.obj.extend([float(s[0]),float(s[1]),float(s[2]),float(s[3])])


   def __init__(self,input):
      self.app_fn = None
      self.l_vert = None
      self.t_colr = None
      self.c_colr = None
      self.l_radi = None
      self.l_norm = None
      self.o_vert = None
      self.tri_flag = 0
      self.cc = 0
      self.obj = []
      for a in range(20):
         input.readline()
      dispatch = [
         None,
         self.tri,
         self.sphere,
         self.cyl,
         None,
         self.cyl,
         None,
         self.tri_normal,
         ]
      ld = len(dispatch)
      while 1:
         l = input.readline()
         if not l:
            break
         v = string.split(l)
         n=int(v[0])
         if(n<ld):
            dd = dispatch[n]
            if dd:
               apply(dd,(input,))

      self.append_last()
      if self.tri_flag:
         self.obj.append(END)
      input.close()
