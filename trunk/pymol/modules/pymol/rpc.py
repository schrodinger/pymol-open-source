""" an XML-RPC server to allow remote control of PyMol

  Author: Greg Landrum (Landrum@RationalDiscovery.com)
  Created:  January 2002
  License:  PyMol
  Requires:
            - a python xmlrpclib distribution containing the SimpleXMLRPCServer
              module (1.0 or greater should be fine)
            - python with threading enabled  
             
            
"""
import SimpleXMLRPCServer
import threading
from pymol import cmd,cgo

# initial port to try for the server
_xmlPort=9123
# number of alternate ports to try if the first fails
_nPortsToTry=5
def rpcCmd(cmdText):
  """ executes a PyMol API command

   return value is either the result of the command or the empty string

  """
  res = cmd.do(cmdText)
  if res is not None:
    return res
  else:
    return ''

def rpcLabel(pos,labelText,color=(1,1,1),id='lab1'):
  """ create a text label

    Arguments:
      pos: a 3 tuple with the position of the label
      text: a string with the label
      color: a 3 tuple with the color of the label. (1,1,1) is white
      id: (OPTIONAL) the name of the object to be created

    NOTE:
      at the moment this is, how you say, a hack
      Also at the moment the color argument is ignored

  """
  # FIX: this don't work
  r,g,b = color
  x,y,z = pos

  text = "HETATM%5d  C   UNK     1    %8.3f%8.3f%8.3f  1.00 10.00\n"%(1,x,y,z),
  cmd.read_pdbstr(text,id)
  cmd.label("(%s)"%(id),labelText)
  cmd.hide("nonbonded",id)
  return 1


def rpcSphere(pos,rad,color,id='cgo',extend=1):
  """ create a sphere

    Arguments:
      pos: a 3 tuple with the position of the sphere
      rad: a float with the radius
      color: a 3 tuple with the color of the sphere. (1,1,1) is white
      id: (OPTIONAL) the name of the object to be created
      extend: (OPTIONAL) if this is nonzero, the object will be cleared
        before adding the new sphere.  Otherwise the sphere is appended
        to the ojbect

  """
  r,g,b = color
  x,y,z = pos
  if extend:
    obj = cgoDict.get(id,[])
  else:
    obj = []
  obj.extend([cgo.COLOR,r,g,b,cgo.SPHERE,x,y,z,rad])
  cgoDict[id] = obj
  cmd.load_cgo(obj,id,1)
  return 1


def rpcCylinder(end1,end2,rad,color1,id='cgo',color2=None,extend=1):
  """ create a cylinder

    Arguments:
      end1: a 3 tuple with the position of end1 of the sphere
      end2: a 3 tuple with the position of end1 of the sphere
      rad: a float with the radius
      color1: a 3 tuple with the color of end1 of the sphere. (1,1,1) is white
      id: (OPTIONAL) the name of the object to be created
      color2: (OPTIONAL) a 3 tuple with the color of end2 of the sphere. (1,1,1) is white
      extend: (OPTIONAL) if this is nonzero, the object will be cleared
        before adding the new sphere.  Otherwise the sphere is appended
        to the ojbect

    NOTE: the reason that color2 follows id is that I think clients are
    going to be interested in setting the id more often than they are going
    to care about the second color.
    
  """
  global cgoDict
  
  if color2 is None: color2 = color1
  r1,g1,b1 = color1
  r2,g2,b2 = color2
  x1,y1,z1 = end1
  x2,y2,z2 = end2
  if extend:
    obj = cgoDict.get(id,[])
  else:
    obj = []
  obj.extend([cgo.CYLINDER,x1,y1,z1,x2,y2,z2,rad,r1,g1,b1,r2,g2,b2,])
  cgoDict[id] = obj
  cmd.load_cgo(obj,id,1)
  return 1

def launch_XMLRPC(hostname='localhost',port=_xmlPort,nToTry=_nPortsToTry):
  """ launches the xmlrpc server into a separate thread

    Arguments:
      hostname: name of the host for the server
      port: (OPTIONAL) the first port to try for the server
      nToTry: (OPTIONAL) the number of possible ports to try
        (in case the first can't be opened)

  """
  global cgoDict
  cgoDict = {}
  for i in range(nToTry):
    try:
      serv = SimpleXMLRPCServer.SimpleXMLRPCServer((hostname,port+i))
    except:
      serv = None
    else:
      break
  if serv:
    print 'xml-rpc server running on port %d'%(port+i)
    serv.register_function(rpcCmd,'do')
    serv.register_function(rpcSphere,'sphere')
    serv.register_function(rpcCylinder,'cylinder')
    #serv.register_function(rpcLabel,'label')
    t = threading.Thread(target=serv.serve_forever)
    t.setDaemon(1)
    t.start()
  else:
    print 'xml-rpc server could not be started'
