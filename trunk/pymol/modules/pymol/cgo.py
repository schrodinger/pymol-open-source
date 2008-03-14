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
from cmd import DEFAULT_ERROR, DEFAULT_SUCCESS, _raising

POINTS             = 0.0
LINES              = 1.0
LINE_LOOP          = 2.0
LINE_STRIP         = 3.0
TRIANGLES          = 4.0
TRIANGLE_STRIP     = 5.0
TRIANGLE_FAN       = 6.0
#QUADS              = 7.0
#QUAD_STRIP         = 8.0
#POLYGON            = 9.0                                                            

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
ENABLE             = 12.0
DISABLE            = 13.0
SAUSAGE            = 14.0
CUSTOM_CYLINDER    = 15.0
DOTWIDTH           = 16.0
ALPHA_TRIANGLE     = 17.0
ELLIPSOID          = 18.0

#SHAPE_VERTEX       = 16.0
#SHAPE_COLOR        = 17.0
#SHAPE_NORMAL       = 18.0

FONT               = 19.0
FONT_SCALE         = 20.0
FONT_VERTEX        = 21.0
FONT_AXES          = 22.0

CHAR               = 23.0

ALPHA              = 25.0
QUADRIC            = 26.0 # NOTE: Only works with ellipsoids and disks
CONE               = 27.0 

LIGHTING           = float(0x0B50)

def molauto(*arg,**kw):
    _self = kw.get('_self',cmd)
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
    _self.save("molauto.pdb",sele)
    print "molauto %s -nocentre molauto.pdb | molscript -r > molauto.r3d"%marg
    os.system("molauto %s -nocentre molauto.pdb | molscript -r > molauto.r3d"%marg)
    f = open("molauto.r3d")
    rr = RenderReader(f)
    f.close()
    _self.load_cgo(rr.obj,name)

# the following implementation causes full-blown system crashes on some machines.
#   (stdout,stdin) = popen2.popen2("molauto %s -nocentre molauto.pdb | molscript -r > molauto.r3d"%marg)
#
#   if stdin:
#      stdin.close()
#      rr = RenderReader(stdout)
#      cmd.load_cgo(rr.obj,name)

def measure_text(font,text,
                      axes=[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]):
    w = 0
    x = axes[0]
    for char in text:
        if font.has_key(char):
            w = w + font[char][0]*x[0]
    return w

def wire_text(cgo,font,pos,text,
                  axes = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]): # modifies pos
    x = axes[0]
    y = axes[1]
    for char in text:
        if font.has_key(char):
            fc = font[char]
            stroke = 0
            w = fc[0]
            f = fc[1]
            c = 0
            l = len(f)-2
            while c<l:
                if not f[c]:
                    if stroke: cgo.append(END)
                    cgo.append(BEGIN)
                    cgo.append(LINE_STRIP)
                    stroke = 1
                ax = f[c+1]
                ay = f[c+2]
                cgo.append(VERTEX)
                cgo.append(pos[0]+x[0]*ax+y[0]*ay)
                cgo.append(pos[1]+x[1]*ax+y[1]*ay)
                cgo.append(pos[2]+x[2]*ax+y[2]*ay)         
                c = c + 3
            pos[0] = pos[0] + w*x[0]
            pos[1] = pos[1] + w*x[1]
            pos[2] = pos[2] + w*x[2]
            if stroke: cgo.append(END)

def cyl_text(cgo,font,pos,text,radius=0.1,color=[1.0,1.0,1.0],
                  axes = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]): # modifies pos
    x = axes[0]
    y = axes[1]
    for char in text:
        if font.has_key(char):
            fc = font[char]
            stroke = 0
            w = fc[0]
            f = fc[1]
            c = 0
            l = len(f)-2
            while c<l:
                ax = f[c+1]
                ay = f[c+2]
                next = [(pos[0]+x[0]*ax+y[0]*ay),
                          (pos[1]+x[1]*ax+y[1]*ay),
                          (pos[2]+x[2]*ax+y[2]*ay)]
                if f[c]:
                    if stroke:
                        cgo.append(SAUSAGE)
                        cgo.extend(last)
                        cgo.extend(next)
                        cgo.append(radius)
                        cgo.extend(color)
                        cgo.extend(color)
                else:
                    stroke = 1
                last = next
                c = c + 3
            pos[0] = pos[0] + w*x[0]
            pos[1] = pos[1] + w*x[1]
            pos[2] = pos[2] + w*x[2]


def from_r3d(fname):
    result = DEFAULT_ERROR
    input = None
    if string.find(fname,':')>1:
        import urllib
        input = urllib.urlopen(fname)
    elif os.path.exists(fname):
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
            self.obj.append(SAUSAGE)
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

    def quadric(self,f):
        self.append_last()
        l = f.readline()
        if l:
            s = string.split(l)
            self.obj.append(COLOR)
            self.obj.extend([float(s[4]),float(s[5]),float(s[6])])
            self.obj.append(QUADRIC)
            self.obj.extend([float(s[0]),float(s[1]),float(s[2]),float(s[3])])
        l = f.readline()
        if l:
            s = string.split(l)
            self.obj.extend([float(s[0]),float(s[1]),float(s[2]),
                             float(s[3]),float(s[4]),float(s[5]),
                             float(s[6]),float(s[7]),float(s[8]),float(s[9])])
        
    def mat_prop(self,f):
        self.append_last()
        l = f.readline()
        print "mat_prop"+l
        if l:
            s = string.split(l)
            (mphong, mspec, sr, sg, sb, clrity) = map(float,s[0:6])
            if clrity>0.999:
                clrity=0.999
            self.obj.extend([ALPHA, 1.0-clrity])
            opts1 = int(s[6])
            opts4 = int(s[9])
            for x in range(opts4):
                f.readline();

    def mat_reset(self,f):
        self.append_last()
        self.obj.extend([ALPHA, 1.0])

    def reset(self,f):
        pass
    
    def __init__(self,input):
        # Author: Warren DeLano
        # Modifications: Robert Campbell
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
            self.mat_prop,
            self.reset,
            None,
            None,
            None,
            None,
            self.quadric,
            ]
        ld = len(dispatch)
        while 1:
            l = input.readline()
            if not l:
                break
            if l[0] != '#':
                v = string.split(l)
                n=int(v[0])
                if(n<ld):
                    dd = dispatch[n]
                    if dd:
                        apply(dd,(input,))
                    else:
                        # skip over lines that don't match desired object type
                        input.readline()
                elif ( n != 9 ):
                    # don't read another line if render object type 0
                    input.readline()
        self.append_last()
        if self.tri_flag:
            self.obj.append(END)
        input.close()
