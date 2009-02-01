# This demonstration requires Python 2.2 and
# PyMOL compiled with Python 2.2

# Author: Greg Landum
# Commentary by Warren DeLano

# To Try:
# (1) in one process launch PyMOL with "-R" option "./pymol.com -R"
# (2) in another process "python xmlrpc01.py"

molData="""foo


  5  4  0  0  0                 1 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0
    0.5000    0.5000    0.5000 C   0  0  0  0  0
   -0.5000   -0.5000    0.5000 C   0  0  0  0  0
   -0.5000    0.5000   -0.5000 C   0  0  0  0  0
    0.5000   -0.5000   -0.5000 C   0  0  0  0  0
  1  2  1  0  0  0
  1  3  1  0  0  0
  1  4  1  0  0  0
  1  5  1  0  0  0
M  END
"""
import sys,xmlrpclib,tempfile,os

fName = tempfile.mktemp('.mol')
open(fName,'w+').write(molData)
port = 9123
if len(sys.argv)>1:
  port = int(sys.argv[1])

s = xmlrpclib.Server('http://localhost:%d'%(port))
s.do('delete fooby')
s.do('delete fooby-o')
s.do('load %s, fooby'%(fName))
s.sphere((.5,.5,.5),.5,(1,.5,1),'fooby-o',0)
s.sphere((-.5,-.5,.5),.5,(1,1,.5),'fooby-o')
s.sphere((-.5,.5,-.5),.5,(1,.5,5),'fooby-o')
s.sphere((.5,-.5,-.5),.5,(.5,1,.5),'fooby-o')

s.cylinder((.5,.5,.5),(-.5,-.5,.5),.1,(.5,.5,.5),'fooby-o')
s.cylinder((.5,.5,.5),(-.5,.5,-.5),.1,(.5,.5,.5),'fooby-o')
s.cylinder((.5,.5,.5),(.5,-.5,-.5),.1,(.5,.5,.5),'fooby-o')
s.cylinder((-.5,-.5,.5),(.5,-.5,-.5),.1,(.5,.5,.5),'fooby-o')
s.cylinder((-.5,-.5,.5),(-.5,.5,-.5),.1,(.5,.5,.5),'fooby-o')
s.cylinder((.5,-.5,-.5),(-.5,.5,-.5),.1,(.5,.5,.5),'fooby-o')

#s.label((1,0,0),"foobar")

os.unlink(fName)
