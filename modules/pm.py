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

# pm.py 
# Python interface module for PyMol
#
# **This is the only module which should be/need be imported by 
# **PyMol Programs

import re
import _pm
import string
import traceback
import thread
import threading
import __main__

def ready():
   return _pm.ready()

def setup_global_locks():
   __main__.lock_api = _pm.get_globals()['lock_api']
   
def lock():
   __main__.lock_api.acquire(1)

def lock_attempt():
   res = __main__.lock_api.acquire(blocking=0)
   if res:
      _pm.get_globals()['lock_state'] = 1;
   else:
      _pm.get_globals()['lock_state'] = None;

def unlock():
   __main__.lock_api.release()

def export_dots(a,b):
   lock()
   r = _pm.export_dots(a,int(b))
   unlock()
   return r

def do(a):
   lock()
   _pm.do(a);
   unlock()
   
def turn(a,b):
   lock()
   _pm.turn(a,float(b))
   unlock()
   
def render():
   lock()   
   _pm.render()
   unlock()
   
def ray():
   lock()   
   _pm.render()
   unlock()

def real_system(a):
   _pm.system(a)
   
def system(a):
   threading.Thread(target=real_system,args=([a])).start()
   lock()
   unlock()

def fit(a,b):
   if a[0]!='(': a="(%"+a+")"
   if b[0]!='(': b="(%"+b+")"
   lock()   
   _pm.fit("(%s in %s)" % (a,b),
         "(%s in %s)" % (b,a))
   unlock()

def expfit(a,b):
   lock()   
   _pm.fit(a,b)
   unlock()
   
def zoom(a):
   lock()   
   _pm.zoom(a)
   unlock()
   
def frame(a):
   lock()   
   _pm.frame(int(a))
   unlock()
   
def move(a,b):
   lock()   
   _pm.move(a,float(b))
   unlock()
   

def clip(a,b):
   lock()   
   _pm.clip(a,float(b))
   unlock()
   

def origin(a):
   lock()   
   _pm.origin(a)
   unlock()

def orient(*arg):
   lock()
   if len(arg)<1:
      a = "(all)"
   else:
      a = arg[0]
   _pm.orient(a)
   unlock()
   

def refresh():
   lock()
   if thread.get_ident() ==__main__.glutThread:
      _pm.refresh_now()
   else:
      _pm.refresh()
   unlock()
   

def dirty():
   lock()   
   _pm.dirty()
   unlock()

def set(a,b):
   lock()   
   _pm.set(a,b)
   unlock()
   

def reset():
   lock()   
   _pm.reset(0)
   unlock()

def reset_rate():
   lock()   
   _pm.reset_rate()
   unlock()
      
def delete(a):
   lock()   
   _pm.delete(a)
   unlock()

def quit_():
   lock()
   _pm.quit()
   unlock()
   
def quit():
   if thread.get_ident() ==__main__.glutThread:
      lock()
      _pm.do("quit_")
   else:
      _pm.do("quit_")
      thread.exit()
   
def png(a):
   lock()   
   fname = a
   if not re.search("\.png$",fname):
      fname = fname +".png"
   _pm.png(fname)
   unlock()
   

def mclear():
   lock()   
   _pm.mclear()
   unlock()
   

def mstop():
   lock()   
   _pm.mplay(0)
   unlock()
   

def mplay():
   lock()   
   _pm.mplay(1)
   unlock()
   

def mray():
   lock()   
   _pm.mplay(2)
   unlock()
   

def viewport(a,b):
   lock()   
   _pm.viewport(int(a),int(b))
   unlock()
   

def mdo(a,b):
   lock()   
   _pm.mdo(int(a)-1,b)
   unlock()
   

def dummy(*args):
   lock()   
   pass
   unlock()
   

def rock():
   lock()   
   _pm.rock()
   unlock()
   

def forward():
   lock()   
   _pm.setframe(5,1)
   unlock()

def backward():
   lock()   
   _pm.setframe(5,-1)
   unlock()

def beginning():
   lock()   
   _pm.setframe(0,0)
   unlock()

def ending():
   lock()   
   _pm.setframe(2,0)
   unlock()
         
def middle():
   lock()   
   _pm.setframe(3,0)
   unlock()
   
def save(*arg):
   lock()
   fname = 'save.pdb'
   sele = '( all )'
   state = -1
   format = 'pdb'
   if len(arg)==1:
      fname = arg[0]
   elif len(arg)==2:
      fname = arg[0]
      sele = arg[1]
   elif len(arg)==3:
      fname = arg[0]
      sele = arg[1]
      state = arg[2]
   elif len(arg)==4:
      fname = arg[0]
      sele = arg[1]
      state = arg[2]
      format = arg[3]
   if (len(arg)>0) and (len(arg)<4):
      if re.search("\.pdb$",fname):
         format = 'pdb'
      elif re.search("\.mol$",fname):
         formet = 'mol'
      elif re.search("\.sdf$",fname):
         formet = 'sdf'
   if format=='pdb':
      f=open(fname,"w")
      if f:
         f.write(_pm.get_pdb(sele,int(state)))
         f.close()
      
   unlock()

def get_feedback():
   l = []
   lock()
   unlock()
   r = _pm.get_feedback()
   while r:
      l.append(r)
      r = _pm.get_feedback()
   return l

def load(*args):
   lock()   
   ftype = 0
   if re.search("\.pdb$",args[0]):
      ftype = 0
   elif re.search("\.mol$",args[0]):
      ftype = 1
   if len(args)==1:
      oname = re.sub("[^/]*\/","",args[0])
      oname = re.sub("\.pdb|\.mol","",oname)
      _pm.load(oname,args[0],-1,ftype)
   elif len(args)==2:
      oname = string.strip(args[1])
      _pm.load(oname,args[0],-1,ftype)
   elif len(args)==3:
      oname = string.strip(args[1])
      _pm.load(oname,args[0],int(args[2])-1,ftype)
   else:
      print "argument error."
   unlock()
   

def read_mol(*args):
   lock()   
   ftype = 3
   if len(args)==2:
      oname = string.strip(args[1])
      _pm.load(oname,args[0],-1,ftype)
   elif len(args)==3:
      oname = string.strip(args[1])
      _pm.load(oname,args[0],int(args[2])-1,ftype)
   else:
      print "argument error."
   unlock()
   

def select(*args):
   lock()   
   if len(args)==1:
      sel_cnt = _pm.get("sel_counter") + 1.0
      _pm.set("sel_counter","%1.0f" % sel_cnt)
      sel_name = "sel%02.0f" % sel_cnt
      sel = args[0]
   else:
      sel_name = args[0]
      sel = args[1]
   _pm.select(sel_name,sel)
   unlock()
   

def color(*args):
   lock()   
   if len(args)==2:
      _pm.color(args[0],args[1],0)
   else:
      _pm.color(args[0],args[1],1)   
   unlock()
   

def mpng(a):
   _pm.do("pm.mpng_('"+a+"')")

def mpng_(*args):
   lock()   
   fname = args[0]
   if re.search("\.png$",fname):
      fname = re.sub("\.png$","",fname)
   _pm.mpng_(fname)
   unlock()
   
def show(*args):
   lock()   
   if len(args)==2:
      if repres.has_key(args[0]):      
         repn = repres[args[0]];
         _pm.showhide(args[1],repn,1);
   elif args[0]=='all':
         _pm.showhide("!",0,1);
   unlock()
   

def hide(*args):
   lock()   
   if len(args)==2:
      if repres.has_key(args[0]):      
         repn = repres[args[0]];
         _pm.showhide(args[1],repn,0);
   elif args[0]=='all':
         _pm.showhide("!",0,0);
   unlock()
   

def mmatrix(a):
   lock()   
   if a=="clear":
      _pm.mmatrix(0)
   elif a=="store":
      _pm.mmatrix(1)
   elif a=="recall":
      _pm.mmatrix(2)
   unlock()
   

def enable(nam):
   lock()   
   _pm.onoff(nam,1);
   unlock()
   

def disable(nam):
   lock()   
   _pm.onoff(nam,0);
   unlock()
   

def mset(seq):
   lock()   
   output=[]
   input = string.split(seq," ")
   last = 0
   for x in input:
      if x[0]>"9" or x[0]<"0":
         if x[0]=="x":
            cnt = int(x[1:])-1
            while cnt>0:
               output.append(str(last))
               cnt=cnt-1
         elif x[0]=="-":
            dir=1
            cnt=last
            last = int(x[1:])-1
            if last<cnt:
               dir=-1
            while cnt!=last:
               cnt=cnt+dir
               output.append(str(cnt))
      else:
         val = int(x) - 1
         output.append(str(val))
         last=val
   _pm.mset(string.join(output," "))
   unlock()
   
   
keyword = { 
   'backward'    : [backward     , 0 , 0 , ',' , 0 ],
   'beginning'   : [beginning    , 0 , 0 , ',' , 0 ],
   'clip'        : [clip         , 2 , 2 , ',' , 0 ],
   'color'       : [color        , 2 , 3 , ',' , 0 ],
   'disable'     : [disable      , 1 , 1 , ',' , 0 ],
   'delete'      : [delete       , 1 , 1 , ',' , 0 ],
   'enable'      : [enable       , 1 , 1 , ',' , 0 ],
   'export_dots' : [export_dots  , 2 , 2 , ',' , 0 ],
   'fit'         : [fit          , 2 , 2 , ',' , 0 ],
   'fork'        : [dummy        , 1 , 1 , ',' , 3 ],
   'forward'     : [forward      , 0 , 0 , ',' , 0 ],
   'frame'       : [frame        , 1 , 1 , ',' , 0 ],
   'hide'        : [hide         , 1 , 2 , ',' , 0 ],
   'load'        : [load         , 1 , 3 , ',' , 0 ],
   'move'        : [move         , 2 , 2 , ',' , 0 ],
   'mset'        : [mset         , 1 , 1 , ',' , 0 ],
   'mdo'         : [mdo          , 2 , 2 , ':' , 1 ],
   'mpng'        : [mpng         , 1 , 2 , ',' , 0 ],
   'mplay'       : [mplay        , 0 , 0 , ',' , 0 ],
   'mray'        : [mray         , 0 , 0 , ',' , 0 ],
   'mstop'       : [mstop        , 0 , 0 , ',' , 0 ],
   'mclear'      : [mclear       , 0 , 0 , ',' , 0 ],
   'middle'      : [middle       , 0 , 0 , ',' , 0 ],
   'mmatrix'     : [mmatrix      , 1 , 1 , ',' , 0 ],
   'origin'      : [origin       , 1 , 1 , ',' , 0 ],
   'orient'      : [orient       , 0 , 1 , ',' , 0 ],
   'ray'         : [render       , 0 , 0 , ',' , 0 ],
   'refresh'     : [refresh      , 0 , 0 , ',' , 0 ],
   'render'      : [render       , 0 , 0 , ',' , 0 ],
   'reset'       : [reset        , 0 , 0 , ',' , 0 ],
   'reset_rate'  : [reset_rate   , 0 , 0 , ',' , 0 ],
   'rewind'      : [beginning    , 0 , 0 , ',' , 0 ],
   'rock'        : [rock         , 0 , 0 , ',' , 0 ],
   'run'         : [dummy        , 1 , 2 , ',' , 2 ],
   'select'      : [select       , 1 , 2 , '=' , 0 ],
   'set'         : [set          , 2 , 2 , '=' , 0 ],
   'show'        : [show         , 1 , 2 , ',' , 0 ],
   'system'      : [system       , 1 , 1 , ',' , 0 ],
   'turn'        : [turn         , 2 , 2 , ',' , 0 ],
   'quit'        : [quit         , 0 , 0 , ',' , 0 ],
   'quit_'       : [quit_        , 0 , 0 , ',' , 0 ],
   'png'         : [png          , 1 , 1 , ',' , 0 ],
   'viewport'    : [viewport     , 2 , 2 , ',' , 0 ],
   'save'        : [save         , 0 , 4 , ',' , 0 ],
   'zoom'        : [zoom         , 1 , 1 , ',' , 0 ]
   }

repres = {
   'lines'       : 0,
   'sticks'      : 1,
   'dots'        : 2,
   'mesh'        : 3,
   'spheres'     : 4,
   'ribbon'      : 5,
   'surface'     : 6
}



