#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* Copyright (c) Schrodinger, LLC.
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

if True:

    import os
    import traceback

    povray_exe = "povray"

    def render_from_string(header,pov_inp,prefix,width,height,antialias):
        r = None
        try:
            pov = prefix +".pov"
            png = prefix +".png"
            f=open(pov,'w')
            f.write(header)
            f.write(pov_inp)
            f.close()
            if os.path.exists(png):
                os.unlink(png)
            if antialias:
                ant="+A"
            else:
                ant=""
            os.system("%s -D +I%s +O%s +W%d +H%d %s"%(
                         povray_exe,pov,png,width,height,ant))
            if os.path.exists(png):
                r = 1
        except:
            traceback.print_exc()
        return r
