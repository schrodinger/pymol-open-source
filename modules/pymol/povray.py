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

import os
import traceback

povray_exe = "x-povray"

def render_from_string(pov_inp,prefix):
   r = None
   try:
      pov = prefix +".pov"
      png = prefix +".png"
      f=open(pov,'w')
      f.write(pov_inp)
      f.close()
      if os.path.exists(png):
         os.unlink(png)
      os.system("%s +I%s +O%s +W640 +H480 -A0.01"%(
                povray_exe,pov,png))
      if os.path.exists(png):
         r = 1
   except:
      traceback.print_exc()
   return r
