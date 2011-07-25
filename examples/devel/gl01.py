from OpenGL.GL import *
from OpenGL.GL.ARB.vertex_buffer_object import *
from pymol.callback import Callback
from pymol import cmd

# this is a trivial example of creating a callback object which
# makes OpenGL calls using the PyOpenGL bindings

# define a callback object

class myCallback(Callback):

   def __call__(self):
      
      glBegin(GL_LINES)
      
      glColor3f(1.0,1.0,1.0)
      
      glVertex3f(0.0,0.0,0.0)
      glVertex3f(1.0,0.0,0.0)

      glVertex3f(0.0,0.0,0.0)
      glVertex3f(0.0,2.0,0.0)

      glVertex3f(0.0,0.0,0.0)
      glVertex3f(0.0,0.0,3.0)

      glEnd()

   def get_extent(self):
      return [[0.0,0.0,0.0],[1.0,2.0,3.0]]

# load it into PyMOL

cmd.load_callback(myCallback(),'gl01')




