#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* Copyright (c) Schrodinger, LLC.
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

from chempy import cpv
#import popen2
import os
from pymol import cmd
from .cmd import DEFAULT_ERROR, DEFAULT_SUCCESS, _raising

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

PICK_COLOR         = 31.0 # 0x1F [PICK_COLOR, index, bond/cPickable_t]

LIGHTING           = float(0x0B50)

# enum cPickable_t:
cPickableAtom = -1.0
cPickableLabel = -2.0
cPickableGadget = -3.0
cPickableNoPick = -4.0
cPickableThrough = -5.0


def molauto(name="mols", sele="(all)", marg="-nice", _self=cmd):
    _self.save("molauto.pdb",sele)
    print("molauto %s -nocentre molauto.pdb | molscript -r > molauto.r3d"%marg)
    os.system("molauto %s -nocentre molauto.pdb | molscript -r > molauto.r3d"%marg)
    f = open("molauto.r3d")
    rr = RenderReader(f)
    f.close()
    _self.load_cgo(rr.obj,name)

def measure_text(font,text,
                      axes=[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]):
    w = 0
    x = axes[0]
    for char in text:
        if char in font:
            w = w + font[char][0]*x[0]
    return w

def wire_text(cgo,font,pos,text,
                  axes = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]): # modifies pos
    x = axes[0]
    y = axes[1]
    for char in text:
        if char in font:
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
        if char in font:
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
    if '://' in fname:
        from urllib.request import urlopen
        input = urlopen(fname)
    elif os.path.exists(fname):
        input = open(fname)
    if input:
        rr = RenderReader(input)
        result = rr.obj
    return result

class RenderReader:

    def append_last(self):
        if self.app_fn:
            self.app_fn()
            self.app_fn=None

    def append_tri(self):
        if self.l_vert:
            d0 = cpv.sub(self.l_vert[0],self.l_vert[1])
            d1 = cpv.sub(self.l_vert[0],self.l_vert[2])
            n0 = cpv.cross_product(d0,d1)
            n0 = cpv.normalize_failsafe(n0)

            if not self.tri_flag:
                self.obj.append(BEGIN)
                self.obj.append(TRIANGLES)
                self.tri_flag = 1

            indices = [0, 1, 2]

            if not self.l_norm:
                # TODO could simplify this if ray tracing would support
                # object-level two_sided_lighting. Duplicating the
                # face with an offset is a hack and produces visible
                # lines on edges.
                n1 = [-n0[0],-n0[1],-n0[2]]
                ns = cpv.scale(n0,0.002)
                indices = [0, 1, 2, 4, 3, 5]
                l_vert_offsetted =     [cpv.add(v, ns) for v in self.l_vert]
                l_vert_offsetted.extend(cpv.sub(v, ns) for v in self.l_vert)
                self.l_vert = l_vert_offsetted
                self.l_norm = [n0, n0, n0, n1, n1, n1]
            elif cpv.dot_product(self.l_norm[0], n0) < 0:
                indices = [0, 2, 1]

            for i in indices:
                self.obj.append(COLOR) # assuming unicolor
                self.obj.extend(self.t_colr[i % 3])
                self.obj.append(NORMAL)
                self.obj.extend(self.l_norm[i])
                self.obj.append(VERTEX)
                self.obj.extend(self.l_vert[i])

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
            s = l.split()
            self.l_vert = [[float(s[0]),float(s[1]),float(s[2])],
                        [float(s[3]),float(s[4]),float(s[5])],
                        [float(s[6]),float(s[7]),float(s[8])]]
            self.t_colr_t = [float(s[9]),float(s[10]),float(s[11])]
            self.t_colr = [self.t_colr_t,self.t_colr_t,self.t_colr_t]

    def tri_normal(self,f):
        l = f.readline()
        if l:
            s = l.split()
            self.l_norm = [[float(s[0]),float(s[1]),float(s[2])],
                        [float(s[3]),float(s[4]),float(s[5])],
                        [float(s[6]),float(s[7]),float(s[8])]]

    def cyl(self,f):
        self.append_last()
        l = f.readline()
        if l:
            self.app_fn = self.append_cyl
            s = l.split()
            self.l_vert = [[float(s[0]),float(s[1]),float(s[2])],
                        [float(s[4]),float(s[5]),float(s[6])]]
            self.l_radi = float(s[3])
            self.c_colr_t = [float(s[8]),float(s[9]),float(s[10])]
            self.c_colr = [self.c_colr_t,self.c_colr_t]

    def sphere(self,f):
        self.append_last()
        l = f.readline()
        if l:
            s = l.split()
            self.obj.append(COLOR)
            self.obj.extend([float(s[4]),float(s[5]),float(s[6])])
            self.obj.append(SPHERE)
            self.obj.extend([float(s[0]),float(s[1]),float(s[2]),float(s[3])])

    def quadric(self,f):
        self.append_last()
        l = f.readline()
        if l:
            s = l.split()
            self.obj.append(COLOR)
            self.obj.extend([float(s[4]),float(s[5]),float(s[6])])
            self.obj.append(QUADRIC)
            self.obj.extend([float(s[0]),float(s[1]),float(s[2]),float(s[3])])
        l = f.readline()
        if l:
            s = l.split()
            self.obj.extend([float(s[0]),float(s[1]),float(s[2]),
                             float(s[3]),float(s[4]),float(s[5]),
                             float(s[6]),float(s[7]),float(s[8]),float(s[9])])

    def mat_prop(self,f):
        self.append_last()
        l = f.readline()
        print("mat_prop"+l)
        if l:
            s = l.split()
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
                v = l.split()
                n=int(v[0])
                if(n<ld):
                    dd = dispatch[n]
                    if dd:
                        dd(input)
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


def torus(center=(0., 0., 0.), normal=(0., 0., 1.), radius=1., color='',
        cradius=.25, samples=20, csamples=20, _self=cmd):
    '''
    Generate and return a torus CGO with given center, normal
    and ring radius.
    '''
    from math import cos, sin, pi

    if color and isinstance(color, str):
        color = list(_self.get_color_tuple(color))
    obj = []

    axis = cpv.cross_product(normal, (0., 0., 1.))
    angle = -cpv.get_angle(normal, (0., 0., 1.))
    matrix = cpv.rotation_matrix(angle, cpv.normalize(axis))

    obj_vertex = lambda x, y, z: obj.extend([VERTEX] + cpv.add(center,
        cpv.transform(matrix, [x, y, z])))
    obj_normal = lambda x, y, z: obj.extend([NORMAL] +
        cpv.transform(matrix, [x, y, z]))

    r = radius
    cr = cradius
    rr = 1.5 * cr
    dv = 2 * pi / csamples
    dw = 2 * pi / samples
    v = 0.0
    w = 0.0

    while w < 2 * pi:
        v = 0.0
        c_w = cos(w)
        s_w = sin(w)
        c_wdw = cos(w + dw)
        s_wdw = sin(w + dw)

        obj.append(BEGIN)
        obj.append(TRIANGLE_STRIP)

        if color:
            obj.append(COLOR)
            obj.extend(color)

        while v < 2 * pi + dv:
            c_v = cos(v)
            s_v = sin(v)
            c_vdv = cos(v + dv)
            s_vdv = sin(v + dv)
            obj_normal(
                (r + rr * c_v) * c_w - (r + cr * c_v) * c_w,
                (r + rr * c_v) * s_w - (r + cr * c_v) * s_w,
                (rr * s_v - cr * s_v))
            obj_vertex(
                (r + cr * c_v) * c_w,
                (r + cr * c_v) * s_w,
                cr * s_v)
            obj_normal(
                (r + rr * c_vdv) * c_wdw - (r + cr * c_vdv) * c_wdw,
                (r + rr * c_vdv) * s_wdw - (r + cr * c_vdv) * s_wdw,
                rr * s_vdv - cr * s_vdv)
            obj_vertex(
                (r + cr * c_vdv) * c_wdw,
                (r + cr * c_vdv) * s_wdw,
                cr * s_vdv)
            v += dv

        obj.append(END)
        w += dw

    return obj

def from_plystr(contents, surfacenormals=True, alphaunit=1.):
    '''
    PLY - Polygon File Format
    '''
    import sys
    if isinstance(contents, bytes):
        contents = contents.decode(errors='ignore')

    lines_iter = iter(contents.splitlines())

    if next(lines_iter) != 'ply':
        raise ValueError('not a ply file')

    types = {'char': int, 'uchar': int, 'short': int, 'ushort': int,
            'int': int, 'uint': int, 'float': float, 'double': float,
            'int8': int, 'uint8': int, 'int16': int, 'uint16': int,
            'int32': int, 'uint32': int, 'float32': float, 'float64': float}

    elements = []

    for line in lines_iter:
        a = line.split()

        if not a:
            continue

        command = a[0]

        if command == 'end_header':
            break

        if command == 'format' or command == 'comment':
            continue

        if command == 'element':
            element = {'type': a[1], 'count': int(a[2]), 'properties': []}
            elements.append(element)
            continue

        if command == 'property':
            element['properties'].append(a[1:])
            continue

        if command == 'obj_info':
            # ignore
            continue

        print('unknown instruction: ' + command)

    table = {}

    for element in elements:
        table[element['type']] = records = []
        properties = element['properties']
        for i in range(element['count']):
            line = next(lines_iter)
            a_iter = iter(line.split())
            rec = {}
            for prop in properties:
                if prop[0] == 'list':
                    rec[prop[-1]] = [types[prop[2]](next(a_iter))
                        for j in range(types[prop[2]](next(a_iter)))]
                else:
                    rec[prop[-1]] = types[prop[0]](next(a_iter))
                    # auto-detect alpha divider
                    if alphaunit == 1. and prop[-1] == 'alpha' and rec['alpha'] > 1:
                        alphaunit = 255.
                    # detect if normals need to be calculated
                    if surfacenormals and prop[-1] == 'nx':
                        surfacenormals = False
            records.append(rec)

    vertices = table['vertex']

    obj = []

    enum_quad = [0, 1, 2, 2, 3, 0]

    def colorfromelem(elem):
        if 'red' in elem:
            obj.append(COLOR)
            obj.append(elem['red'] / 255.0)
            obj.append(elem['green'] / 255.0)
            obj.append(elem['blue'] / 255.0)
        if 'alpha' in elem:
            obj.append(ALPHA)
            obj.append(elem['alpha'] / alphaunit)

    def compute_surface_normals():
        '''
        Compute average normals from all adjacent triangles
        on each vertex
        '''
        from functools import reduce

        # don't use cpv.normalize which has an RSMALL4 limit
        normalize = lambda v: cpv.scale(v, 1. / cpv.length(v))

        for face in table['face']:
            if 'vertex_index' in face:
                indices = face['vertex_index']
            elif 'vertex_indices' in face:
                indices = face['vertex_indices']
            else:
                return

            f_vert = [vertices[i] for i in indices]
            f_xyz = [(v['x'], v['y'], v['z']) for v in f_vert]

            try:
                normal = normalize(cpv.cross_product(
                    cpv.sub(f_xyz[1], f_xyz[0]),
                    cpv.sub(f_xyz[2], f_xyz[1])))
            except ZeroDivisionError:
                continue

            for v in f_vert:
                v.setdefault('normals', []).append(normal)

        for v in vertices:
            try:
                v['nx'], v['ny'], v['nz'] = normalize(
                        reduce(cpv.add, v.pop('normals')))
            except (KeyError, ZeroDivisionError):
                continue

    if surfacenormals and 'face' in table:
        compute_surface_normals()

    for ptype, op in [
            ('face', TRIANGLES),
            ('edge', LINES),
            ('range_grid', POINTS),
            ]:
        if ptype in table:
            obj.append(BEGIN)
            obj.append(op)

            for face in table[ptype]:
                colorfromelem(face)

                if 'vertex_index' in face:
                    indices = face['vertex_index']
                elif 'vertex_indices' in face:
                    indices = face['vertex_indices']
                else:
                    indices = [face['vertex1'], face['vertex2']]

                N = len(indices)

                for i in (enum_quad if N == 4 else range(N)):
                    vertex = vertices[indices[i]]
                    colorfromelem(vertex)

                    if 'nx' in vertex:
                        obj.append(NORMAL)
                        obj.append(vertex['nx'])
                        obj.append(vertex['ny'])
                        obj.append(vertex['nz'])

                    obj.append(VERTEX)
                    obj.append(vertex['x'] * 100)
                    obj.append(vertex['y'] * 100)
                    obj.append(vertex['z'] * 100)
            obj.append(END)

    obj.append(STOP)
    return obj
