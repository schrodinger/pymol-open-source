from pymol.opengl.gl import *
from pymol.callback import Callback
from pymol import cmd

# this is an example of creating a callback object which
# uses GL arrays.  Requires Numeric Python:

# NOTE: NUMERIC NOT SUPPORTED BY PYMOL 1.1 and up

from Numeric import *

# define a callback object

class myCallback(Callback):

   def __init__(self):
      self.vert = zeros([126,3],Float)
      self.norm = zeros([126,3],Float)

      for a in range(0,63):
         x = a/10.0
         self.vert[2*a][0]=x
         self.vert[2*a][1]=math.sin(x)
         self.vert[2*a][2]=0.0

         self.vert[2*a+1][0]=self.vert[2*a][0]
         self.vert[2*a+1][1]=self.vert[2*a][1]
         self.vert[2*a+1][2]=1.0

         dydx = -math.cos(x)*sqrt(2)/2.0# is this right? I can't remember...
         self.norm[2*a][0]=dydx
         self.norm[2*a][1]=math.sqrt(1.0-dydx*dydx)
         self.norm[2*a][2]=0.0

         self.norm[2*a+1][0]=self.norm[2*a][0]
         self.norm[2*a+1][1]=self.norm[2*a][1]
         self.norm[2*a+1][2]=self.norm[2*a][2]

   def __call__(self):
      glColor3f(1.0,1.0,1.0)
      glVertexPointer(3,0,self.vert)
      glNormalPointer(3,0,self.norm)
      
      glEnable(GL_VERTEX_ARRAY)
      glEnable(GL_NORMAL_ARRAY)
      glDisable(GL_COLOR_ARRAY)
      glDrawArrays(GL_TRIANGLE_STRIP,0,126)

   def get_extent(self):
      return [[0.0,0.0,-1.0],
              [6.3,1.0,1.0]]


# load it into PyMOL

cmd.load_callback(myCallback(),'gl02')


# give us a nice view

cmd.turn('z',20)
cmd.turn('y',20)
cmd.turn('x',20)



