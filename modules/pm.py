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
import __main__
import os


def sort(*args):
   lock()
   if len(args)==0:
      r = _pm.sort("")
   else:
      r = _pm.sort(args[0])
   unlock()
   return r

def mem():
   lock()
   r = _pm.mem()
   unlock()
   return r

def ready():
   return _pm.ready()

def copy(src,dst):
   lock()
   r = _pm.copy(src,dst)
   unlock()
   return r

def alter(sele,expr):
   lock()
   r = _pm.alter(sele,expr)
   unlock()   
   return r

def _stereo(flag):
   if flag:
      os.system("/usr/gfx/setmon -n 1024x768_96s")
   else:
      os.system("/usr/gfx/setmon -n 72hz")

def stereo(a):
   r = None
   if a=="on":
      lock()
      if _pm.stereo(1):
         r = _stereo(1)
      else:
         print " error: stereo not available"
      unlock();
   else:
      lock()
      if _pm.stereo(0):
         r = _stereo(0)
      unlock();
   return r
   
def overlap(*args):
   state = [0,0]
   if len(args)==3:
      state[0]=int(args[2][0])
      if state[0]<1: state[0]=1;
      state[1]=int(args[2][1])
      if state[1]<1: state[1]=1
   lock()
   r = _pm.overlap(args[0],args[1],state[0]-1,state[1]-1)
   unlock()
   return r

def distance(*args):
   la = len(args)
   if la==0:
      a="pk1"
      b="pk3"
   elif la==1:
      a=args[0]
      b="pk1"
   elif la==2:
      a=args[0]
      b=args[1]
   if a[0]!='(': a="(%"+a+")"
   if b[0]!='(': b="(%"+b+")"
   lock()   
   r = _pm.distance(a,b)
   unlock()
   return r

def _alter_do(at):
   ns = {'type': at[1],
         'name': at[2],
         'resn': at[3],
         'chain': at[4],
         'resi': at[5],
         'x': at[6],
         'y': at[7],
         'z': at[8],
         'q': at[9],
         'b': at[10],
         'segi': at[11]}
   exec at[0] in _pm.get_globals(),ns
   type = string.upper(string.strip(ns['type']))
   type = type[:6]
   name = string.upper(string.strip(ns['name']))
   name = name[:4]
   resn = string.upper(string.strip(ns['resn']))
   resn = resn[:3]
   chain = string.upper(string.strip(ns['chain']))
   chain = chain[:1]
   resi = string.upper(string.strip(str(ns['resi'])))
   resi = resi[:4]
   x = float(ns['x'])
   y = float(ns['y'])
   z = float(ns['z'])
   b = float(ns['b'])
   q = float(ns['q'])
   segi = string.strip(ns['segi'])
   segi = segi[:4]
   return [type,name,resn,chain,resi,x,y,z,q,b,segi]

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

def count_states(*args):
   lock()
   if not len(args):
      a = "(all)"
   else:
      a=args[0]
   r = _pm.count_states(a)
   unlock()
   return r

def do(a):
   lock()
   r = _pm.do(a);
   unlock()
   return r

def turn(a,b):
   lock()
   r = _pm.turn(a,float(b))
   unlock()
   return r

def render():
   lock()   
   r = _pm.render()
   unlock()
   return r

def ray():
   lock()   
   r = _pm.render()
   unlock()
   return r

def real_system(a):
   r = _pm.system(a)
   return r

def system(a):
   real_system(a)

def intra_fit(*args):
   lock()
   if len(args)<2:
      b=-1
   else:
      b=int(args[1])-1
   r = _pm.intrafit(args[0],b,2)
   unlock()
   return r

def intra_rms(*args):
   lock()
   if len(args)<2:
      b=-1
   else:
      b=int(args[1])-1
   r = _pm.intrafit(args[0],b,1)
   unlock()
   return r

def intra_rms_cur(*args):
   lock()
   if len(args)<2:
      b=-1
   else:
      b=int(args[1])-1
   r = _pm.intrafit(args[0],b,0)
   unlock()
   return r

def fit(a,b):
   if a[0]!='(': a="(%"+a+")"
   if b[0]!='(': b="(%"+b+")"
   lock()   
   r = _pm.fit("(%s in %s)" % (a,b),
         "(%s in %s)" % (b,a),2)
   unlock()
   return r

def rms(a,b):
   if a[0]!='(': a="(%"+a+")"
   if b[0]!='(': b="(%"+b+")"
   lock()   
   r = _pm.fit("(%s in %s)" % (a,b),
         "(%s in %s)" % (b,a),0)
   unlock()
   return r

def rms_cur(a,b):
   if a[0]!='(': a="(%"+a+")"
   if b[0]!='(': b="(%"+b+")"
   lock()   
   r = _pm.fit("(%s in %s)" % (a,b),
         "(%s in %s)" % (b,a),1)
   unlock()
   return r

def pairfit(*args):
   lock()   
   r = _pm.fit_pairs(args)
   unlock()
   return r

def expfit(a,b):
   lock()   
   r = _pm.fit(a,b,2)
   unlock()
   return r
   
def zoom(a):
   lock()   
   r = _pm.zoom(a)
   unlock()
   return r
   
def frame(a):
   lock()   
   r = _pm.frame(int(a))
   unlock()
   return r

def move(a,b):
   lock()   
   r = _pm.move(a,float(b))
   unlock()
   return r

def clip(a,b):
   lock()   
   r = _pm.clip(a,float(b))
   unlock()
   return r

def origin(a):
   lock()   
   r = _pm.origin(a)
   unlock()
   return r

def orient(*arg):
   lock()
   if len(arg)<1:
      a = "(all)"
   else:
      a = arg[0]
   r = _pm.orient(a)
   unlock()
   return r

def refresh():
   lock()
   if thread.get_ident() ==__main__.glutThread:
      r = _pm.refresh_now()
   else:
      r = _pm.refresh()
   unlock()
   return r

def dirty():
   lock()   
   r = _pm.dirty()
   unlock()
   return r

def set(a,b):
   lock()   
   r = _pm.set(a,b)
   unlock()
   return r

def reset():
   lock()   
   r = _pm.reset(0)
   unlock()
   return r

def reset_rate():
   lock()   
   r = _pm.reset_rate()
   unlock()
   return r

def delete(a):
   lock()   
   r = _pm.delete(a)
   unlock()
   return r

def _quit():
   lock()
   r = _pm.quit()
   unlock()
   return r

def quit():
   if thread.get_ident() ==__main__.glutThread:
      lock()
      r = _pm.do("_quit")
   else:
      r = _pm.do("_quit")
      thread.exit()
   return r

def png(a):
   lock()   
   fname = a
   if not re.search("\.png$",fname):
      fname = fname +".png"
   r = _pm.png(fname)
   unlock()
   return r

def mclear():
   lock()   
   r = _pm.mclear()
   unlock()
   return r

def _special(k,x,y):
   k=int(k)
   if special.has_key(k):
      if special[k][1]:
         if special[k][2]:
            apply(special[k][1],(x,y))
         else:
            apply(special[k][1],())
   return None

def mstop():
   lock()   
   r = _pm.mplay(0)
   unlock()
   return r

def mplay():
   lock()   
   r = _pm.mplay(1)
   unlock()
   return r

def mray():
   lock()   
   r = _pm.mplay(2)
   unlock()
   return r

def viewport(a,b):
   r = _pm.viewport(int(a),int(b))
   
def mdo(a,b):
   lock()   
   r = _pm.mdo(int(a)-1,b)
   unlock()
   return r

def dummy(*args):
   lock()   
   pass
   unlock()
   return None

def rock():
   lock()   
   r = _pm.rock()
   unlock()
   return r

def forward():
   lock()   
   r = _pm.setframe(5,1)
   unlock()
   return r

def backward():
   lock()   
   r = _pm.setframe(5,-1)
   unlock()
   return r

def beginning():
   lock()   
   r = _pm.setframe(0,0)
   unlock()
   return r

def ending():
   lock()   
   r=_pm.setframe(2,0)
   unlock()
   return r

def middle():
   lock()   
   r = _pm.setframe(3,0)
   unlock()
   return r

def save(*arg):
   r = 1
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
         f.write(_pm.get_pdb(sele,int(state)-1))
         f.close()
         r = None
         print " Save: wrote \""+fname+"\"."
   unlock()
   return r

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
   r = 1
   lock()   
   ftype = 0
   if re.search("\.pdb$",args[0]):
      ftype = 0
   elif re.search("\.mol$",args[0]):
      ftype = 1
   elif re.search("\.mmod$",args[0]):
      ftype = 4
   if len(args)==1:
      oname = re.sub("[^/]*\/","",args[0])
      oname = re.sub("\.pdb|\.mol|\.mmod","",oname)
      r = _pm.load(oname,args[0],-1,ftype)
   elif len(args)==2:
      oname = string.strip(args[1])
      r = _pm.load(oname,args[0],-1,ftype)
   elif len(args)==3:
      oname = string.strip(args[1])
      r = _pm.load(oname,args[0],int(args[2])-1,ftype)
   elif len(args)==4:
      oname = string.strip(args[1])
      r = _pm.load(oname,args[0],int(args[2])-1,int(args[3]))
   else:
      print "argument error."
   unlock()
   return r

def read_mol(*args):
   r = 1
   lock()   
   ftype = 3
   if len(args)==2:
      oname = string.strip(args[1])
      r = _pm.load(oname,args[0],-1,ftype)
   elif len(args)==3:
      oname = string.strip(args[1])
      r = _pm.load(oname,args[0],int(args[2])-1,ftype)
   else:
      print "argument error."
   unlock()
   return r

def load_mmolstr(*args):
   r = 1
   lock()   
   ftype = 6
   if len(args)==2:
      oname = string.strip(args[1])
      r = _pm.load(oname,args[0],-1,ftype)
   elif len(args)==3:
      oname = string.strip(args[1])
      r = _pm.load(oname,args[0],int(args[2])-1,ftype)
   else:
      print "argument error."
   unlock()
   return r
   
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
   r = _pm.select(sel_name,sel)
   unlock()
   return r

def color(*args):
   lock()   
   if len(args)==2:
      r = _pm.color(args[0],args[1],0)
   else:
      r = _pm.color(args[0],"(all)",0)   
   unlock()
   return r

def colordef(nam,col):
   r = 1
   lock()
   c = string.split(col)
   if len(c)==3:
      r = _pm.colordef(nam,float(c[0]),float(c[1]),float(c[2]))
   else:
      print "invalid color vector"
   unlock()
   return r

def mpng(a):
   if thread.get_ident() ==__main__.glutThread:
      r = _mpng(a)
   else:
      r = _pm.do("pm._mpng('"+a+"')")
   return r

def _mpng(*args):
   lock()   
   fname = args[0]
   if re.search("\.png$",fname):
      fname = re.sub("\.png$","",fname)
   r = _pm.mpng_(fname)
   unlock()
   return r

def show(*args):
   r=1
   lock()   
   if len(args)==2:
      if repres.has_key(args[0]):      
         repn = repres[args[0]];
         r = _pm.showhide(args[1],repn,1);
   elif args[0]=='all':
         r = _pm.showhide("!",0,1);
   unlock()
   return r

def hide(*args):
   r = 1
   lock()   
   if len(args)==2:
      if repres.has_key(args[0]):      
         repn = repres[args[0]];
         r = _pm.showhide(args[1],repn,0);
   elif args[0]=='all':
         r = _pm.showhide("!",0,0);
   unlock()
   return r

def mmatrix(a):
   r = 1
   lock()   
   if a=="clear":
      r = _pm.mmatrix(0)
   elif a=="store":
      r = _pm.mmatrix(1)
   elif a=="recall":
      r = _pm.mmatrix(2)
   unlock()
   return r

def enable(nam):
   lock()   
   r = _pm.onoff(nam,1);
   unlock()
   return r

def disable(nam):
   lock()   
   r = _pm.onoff(nam,0);
   unlock()
   return r

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
   r = _pm.mset(string.join(output," "))
   unlock()
   return r


keyword = { 
   'alter'         : [alter        , 2 , 2 , ',' , 0 ],
   'backward'      : [backward     , 0 , 0 , ',' , 0 ],
   'beginning'     : [beginning    , 0 , 0 , ',' , 0 ],
   'clip'          : [clip         , 2 , 2 , ',' , 0 ],
   'color'         : [color        , 1 , 2 , ',' , 0 ],
   'colordef'      : [colordef     , 2 , 2 , ',' , 0 ],
   'copy'          : [copy         , 2 , 2 , ',' , 0 ],
   'count_states'  : [count_states , 0 , 1 , ',' , 0 ],   
   'delete'        : [delete       , 1 , 1 , ',' , 0 ],
   'disable'       : [disable      , 1 , 1 , ',' , 0 ],
   'dist'          : [distance     , 0 , 2 , ',' , 0 ],
   'enable'        : [enable       , 1 , 1 , ',' , 0 ],
   'export_dots'   : [export_dots  , 2 , 2 , ',' , 0 ],
   'fit'           : [fit          , 2 , 2 , ',' , 0 ],
   'fork'          : [dummy        , 1 , 1 , ',' , 3 ],
   'forward'       : [forward      , 0 , 0 , ',' , 0 ],
   'frame'         : [frame        , 1 , 1 , ',' , 0 ],
   'hide'          : [hide         , 1 , 2 , ',' , 0 ],
   'intra_fit'     : [intra_fit    , 1 , 2 , ',' , 0 ],
   'intra_rms'     : [intra_rms    , 1 , 2 , ',' , 0 ],
   'intra_rms_cur' : [intra_rms_cur, 1 , 2 , ',' , 0 ],
   'load'          : [load         , 1 , 4 , ',' , 0 ],
   'mem'           : [mem          , 0 , 0 , ',' , 0 ],
   'move'          : [move         , 2 , 2 , ',' , 0 ],
   'mset'          : [mset         , 1 , 1 , ',' , 0 ],
   'mdo'           : [mdo          , 2 , 2 , ':' , 1 ],
   'mpng'          : [mpng         , 1 , 2 , ',' , 0 ],
   'mplay'         : [mplay        , 0 , 0 , ',' , 0 ],
   'mray'          : [mray         , 0 , 0 , ',' , 0 ],
   'mstop'         : [mstop        , 0 , 0 , ',' , 0 ],
   'mclear'        : [mclear       , 0 , 0 , ',' , 0 ],
   'middle'        : [middle       , 0 , 0 , ',' , 0 ],
   'mmatrix'       : [mmatrix      , 1 , 1 , ',' , 0 ],
   'origin'        : [origin       , 1 , 1 , ',' , 0 ],
   'orient'        : [orient       , 0 , 1 , ',' , 0 ],
   'overlap'       : [overlap      , 2 , 3 , ',' , 0 ],
   'pairfit'       : [pairfit      , 2 ,98 , ',' , 0 ],
   'ray'           : [render       , 0 , 0 , ',' , 0 ],
   'refresh'       : [refresh      , 0 , 0 , ',' , 0 ],
   'render'        : [render       , 0 , 0 , ',' , 0 ],
   'reset'         : [reset        , 0 , 0 , ',' , 0 ],
   'reset_rate'    : [reset_rate   , 0 , 0 , ',' , 0 ],
   'rewind'        : [beginning    , 0 , 0 , ',' , 0 ],
   'rock'          : [rock         , 0 , 0 , ',' , 0 ],
   'run'           : [dummy        , 1 , 2 , ',' , 2 ],
   'rms'           : [rms          , 2 , 2 , ',' , 0 ],
   'rms_cur'       : [rms_cur      , 2 , 2 , ',' , 0 ],
   'save'          : [save         , 0 , 4 , ',' , 0 ],
   'select'        : [select       , 1 , 2 , '=' , 0 ],
   'sel'           : [select       , 1 , 2 , '=' , 0 ],
   'set'           : [set          , 2 , 2 , '=' , 0 ],
   'show'          : [show         , 1 , 2 , ',' , 0 ],
   'sort'          : [sort         , 0 , 1 , ',' , 0 ],
   '_special'      : [_special     , 3 , 3 , ',' , 0 ],
   'stereo'        : [stereo       , 1 , 1 , ',' , 0 ],
   'system'        : [system       , 1 , 1 , ',' , 0 ],
   'turn'          : [turn         , 2 , 2 , ',' , 0 ],
   'quit'          : [quit         , 0 , 0 , ',' , 0 ],
   '_quit'         : [_quit        , 0 , 0 , ',' , 0 ],
   'png'           : [png          , 1 , 1 , ',' , 0 ],
   'viewport'      : [viewport     , 2 , 2 , ',' , 0 ],
   'zoom'          : [zoom         , 1 , 1 , ',' , 0 ]
   }

repres = {
   'lines'         : 0,
   'sticks'        : 1,
   'dots'          : 2,
   'mesh'          : 3,
   'spheres'       : 4,
   'ribbon'        : 5,
   'surface'       : 6
}

special = {
   1        : [ 'F1'        , None                   , 0 , None ],
   2        : [ 'F2'        , None                   , 0 , None ],
   3        : [ 'F3'        , None                   , 0 , None ],
   4        : [ 'F4'        , None                   , 0 , None ],
   5        : [ 'F5'        , None                   , 0 , None ],
   6        : [ 'F6'        , None                   , 0 , None ],
   7        : [ 'F7'        , None                   , 0 , None ],
   8        : [ 'F8'        , None                   , 0 , None ],
   9        : [ 'F9'        , None                   , 0 , None ],
   10       : [ 'F10'       , None                   , 0 , None ],
   11       : [ 'F11'       , None                   , 0 , None ],
   12       : [ 'F12'       , None                   , 0 , None ],
   100      : [ 'left'      , backward               , 0 , None ],
   101      : [ 'up'        , None                   , 0 , None ],
   102      : [ 'right'     , forward                , 0 , None ],
   103      : [ 'down'      , None                   , 0 , None ],
   104      : [ 'pgup'      , None                   , 0 , None ],
   105      : [ 'pgdown'    , None                   , 0 , None ],
   106      : [ 'home'      , beginning              , 0 , None ],
   107      : [ 'end'       , ending                 , 0 , None ],
   108      : [ 'insert'    , rock                   , 0 , None ]   
}
