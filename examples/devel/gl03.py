'''
This example creates two callback objects which render a moving wave.
One is a static multi-state object, the other is a dynamic single-state
object which updates coordinates every frame.
'''

import numpy
from OpenGL.GL import *
from pymol.callback import Callback
from pymol import cmd

# define a callback object

class myCallback(Callback):
    '''
    Static callback object. Can be loaded into a multistate object
    '''
    def __init__(self, b=0):
        N = 126

        self.vert = numpy.zeros([N, 3], float)
        self.norm = numpy.zeros([N, 3], float)

        # setup x,z for vertices
        self.vert[0::2, 0] = \
        self.vert[1::2, 0] = numpy.arange(N / 2, dtype=float) / 10.0
        self.vert[1::2, 2] = 1.0

        self.update(b)

    def update(self, b):
        '''
        update y for vertices and x,y for normals
        '''
        a = self.vert[0::2, 0]
        x = a + b / 10.0

        # vertices y
        self.vert[0::2, 1] = \
        self.vert[1::2, 1] = numpy.sin(x)

        # normals x,y
        dydx = -numpy.cos(x) * (numpy.sqrt(2.0) / 2.0)
        self.norm[0::2, 0] = dydx
        self.norm[0::2, 1] = numpy.sqrt(1.0 - dydx**2)
        self.norm[1::2] = self.norm[0::2]

    def renderDrawArrays(self):
        '''
        Render with glDrawArrays.
        This should be faster than renderImmediate.
        '''
        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_NORMAL_ARRAY)
        glDisableClientState(GL_COLOR_ARRAY)

        glColor3f(1.0, 1.0, 1.0)
        glVertexPointer(3, GL_FLOAT, 0, self.vert)
        glNormalPointer(GL_FLOAT, 0, self.norm)
        glDrawArrays(GL_TRIANGLE_STRIP, 0, len(self.vert))

        glDisableClientState(GL_VERTEX_ARRAY)
        glDisableClientState(GL_NORMAL_ARRAY)

    def renderImmediate(self):
        '''
        Render in immediate mode
        '''
        glBegin(GL_TRIANGLE_STRIP)
        glColor3f(1.0, 1.0, 1.0)
        for v, n in zip(self.vert, self.norm):
            glNormal3fv(n)
            glVertex3fv(v)
        glEnd()

    __call__ = renderImmediate

    def get_extent(self):
        return [
            self.vert.min(0).tolist(),
            self.vert.max(0).tolist(),
        ]

class myCallbackDynamic(myCallback):
    '''
    Dynamic object which updates vertex coordinates every frame
    '''
    def __call__(self):
        self.update(cmd.get_frame() * 2)
        myCallback.__call__(self)

# load it into PyMOL

for b in range(0, 63):
    cmd.load_callback(myCallback(b), 'gl03', b + 1)

cmd.load_callback(myCallbackDynamic(), 'gl03dyn', 1)

# give us a nice view

cmd.turn('z',20)
cmd.turn('y',20)
cmd.turn('x',20)

cmd.mset('1-62')
cmd.mplay()

