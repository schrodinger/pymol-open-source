""" simple demo for using the PyMol RPC server

  Author: Greg Landrum (Landrum@RationalDiscovery.com)
  Created:       February 2004
  $LastChangedDate$
  License:  PyMol
  Requires:
            - a python xmlrpclib distribution containing the SimpleXMLRPCServer
              module (1.0 or greater should be fine)
  RD Version: $Rev$            
"""
try:
    import xmlrpclib
except ImportError:
    import xmlrpc.client as xmlrpclib

def startServer(host='localhost',startPort=9123,nToTry=5):
  done = 0
  offset = 0
  while offset < nToTry:
    c = xmlrpclib.Server('http://%s:%d'%(host,startPort+offset))
    try:
      c.ping()
    except:
      print('Failed on port %d, trying another'%(startPort+offset))
      offset = offset + 1
    else:
      done = 1
      break
  if done:
    return c,startPort+offset
  else:
    return None,-1


molBlock="""3d.mol


 27 28  0  0  0                 1 V2000
    1.7032    0.2061   -1.4783 C   0  0  0  0  0
   -1.1754    1.1362    0.5252 C   0  0  0  0  0
   -0.8291   -0.4052   -1.6110 C   0  0  0  0  0
   -1.2900   -0.2398   -0.1445 C   0  0  0  0  0
    1.7764    0.1474    0.0609 C   0  0  0  0  0
    1.2054    1.3597    0.8068 C   0  0  0  0  0
   -0.9224   -1.4043    0.7840 C   0  0  0  0  0
    1.4274   -1.2288    0.6468 C   0  0  0  0  0
    0.2875   -1.2046    1.4805 O   0  0  0  0  0
    0.3400    0.4288   -2.1474 C   0  0  0  0  0
    0.0057    1.8686    0.2581 O   0  0  0  0  0
    2.1441   -0.7361   -1.8853 H   0  0  0  0  0
    2.3779    1.0240   -1.8348 H   0  0  0  0  0
   -1.3300    1.0780    1.6257 H   0  0  0  0  0
   -1.9986    1.7818    0.1317 H   0  0  0  0  0
   -1.6979   -0.1528   -2.2717 H   0  0  0  0  0
   -0.6158   -1.4800   -1.8185 H   0  0  0  0  0
   -2.4054   -0.3977   -0.2240 H   0  0  0  0  0
    2.8858    0.2153    0.2548 H   0  0  0  0  0
    1.9363    2.1991    0.7014 H   0  0  0  0  0
    1.0950    1.2014    1.9010 H   0  0  0  0  0
   -0.9058   -2.3887    0.2658 H   0  0  0  0  0
   -1.6992   -1.5000    1.5819 H   0  0  0  0  0
    2.2587   -1.5596    1.3170 H   0  0  0  0  0
    1.3356   -2.0264   -0.1211 H   0  0  0  0  0
    0.4626    0.1900   -3.2341 H   0  0  0  0  0
    0.0668    1.5081   -2.1636 H   0  0  0  0  0
  1  5  1  1  0  0
  1 10  1  0  0  0
  1 12  1  0  0  0
  1 13  1  0  0  0
  2  4  1  0  0  0
  2 11  1  0  0  0
  2 14  1  1  0  0
  2 15  1  0  0  0
  3  4  1  1  0  0
  3 10  1  0  0  0
  3 16  1  6  0  0
  3 17  1  0  0  0
  4  7  1  1  0  0
  4 18  1  0  0  0
  5  6  1  1  0  0
  5  8  1  0  0  0
  5 19  1  0  0  0
  6 11  1  0  0  0
  6 20  1  0  0  0
  6 21  1  1  0  0
  7  9  1  1  0  0
  7 22  1  6  0  0
  7 23  1  1  0  0
  8  9  1  1  0  0
  8 24  1  1  0  0
  8 25  1  6  0  0
 10 26  1  6  0  0
 10 27  1  0  0  0
M  END
"""

if __name__=='__main__':
  import sys
  serv,port = startServer()
  if serv is not None:
    print('connected to PyMol rpc-server on port %d'%port)
  else:
    print('unable to connect to PyMol')
    sys.exit(-1)
  serv.loadMolBlock(molBlock,'sample-mol')
  serv.set('sphere_scale',0.25,'sample-mol')
  serv.do('show sticks;show spheres')
  serv.sphere((.28,-1.2,1.48),.5,(1,0,1),'demo')
  serv.sphere((0,1.87,.26),.5,(1,0,1),'demo')
  serv.cylinder((.28,-1.2,1.48),(0,1.87,.26),.1,(.5,0,.5),'demo')
  serv.zoom()
  
