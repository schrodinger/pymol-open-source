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
# ** PyMol API Based Programs

import re
import _pm
import string
import traceback
import thread
import __main__
import os

def commands():
   '''
COMMANDS
 
   INPUT/OUTPUT  load     save     delete   quit
   VIEW          turn     move     clip     rock
                 show     hide     enable   disable
                 reset    refresh
                 zoom     origin   orient
   MOVIES        mplay    mstop    mset     mdo
                 mpng     mmatrix  frame
                 rewind   middle   ending
                 forward  backward
   IMAGING       png      mpng
   RAYTRACING    ray      
   MAPS          isomesh  isodot
   DISPLAY       cls      viewport splash
   SELECTIONS    select
   SETTINGS      set
   ATOMS         alter
   FITTING       fit      rms      rms_cur  intra_fit    
   COLORS        color    colordef
   HELP          help     commands
   DISTANCES     dist      
   STEREO        stereo
   SYMMETRY      symexp
 
Try "help <command-name>" for more information on a given command.
Additional help is also available for "selections", for
special "api" commands, and "keyboard" for keyboard shortcuts.
   '''
   help('commands')

def api():
   '''
API COMMANDS
 
   KEY BINDING   set_key

   '''
   help('api')

   
def keyboard():
   '''
KEYBOARD COMMANDS and MODIFIERS
 
   TAB          Toggle onscreen text.
 
   INSERT       Toggle rocking.
 
   LEFT ARROW   Go backwards one frame.
   RIGHT ARROW  Go forwards one frame.
   END          Go to end of movie.
   HOME         Go to beginning of movie.
 
 ATOM SELECTIONS (These only work on the "lines" representation!)
  
   CTRL/left mouse click    Pick atom and store as selection (pk1).
   CTRL/middle mouse click  Pick atom and store as selection (pk2).
   CTRL/right mouse click   Pick atom and store as selection (pk3).
 
   CTRL-SHIFT/left mouse click    Pick atom and add to selection (pk1).
   CTRL-SHIFT/middle mouse click  Pick atom and add to selection (pk2).
   CTRL-SHIFT/right mouse click   Pick atom and add to selection (pk3).
   '''
   help('keyboard')

def mouse():
   '''
MOUSE CONTROLS
 
   Transformations can be changed by setting the buttom_mode variable.
   The current configuration is visible on screen with the following
   abbreviations:
 
      R-XYZ    = Rotates about X, Y, and Z axes
      R-Z      = Rotates about the Z axis
      Trans-XY = Translates along the X and Y axes
      Trans-Z  = Translates along Z axis
      Clip-NF  = Y motion moves the near clipping plane while
                 X motion moves the far one
      Clip-N   = Motion affects only the near clipping plane
      Clip-F   = Motion affects only the far clipping plane
 
ATOM SELECTIONS (These only work on the "lines" representation!)
 
   CTRL-left mouse click    Pick atom and store as selection (pk1).
   CTRL-middle mouse click  Pick atom and store as selection (pk2).
   CTRL-right mouse click   Pick atom and store as selection (pk3).
 
   CTRL-SHIFT-left mouse click    Pick atom and add to selection (pk1).
   CTRL-SHIFT-middle mouse click  Pick atom and add to selection (pk2).
   CTRL-SHIFT-right mouse click   Pick atom and add to selection (pk3).
   '''
   help('mouse')

def examples():
   '''
EXAMPLE ATOM SELECTIONS
 
   select bk = ( name ca or name c or name n )
      * can be abbreviated as *
   sel bk = (n;ca,c,n)
 
   select hev = ( not hydro )
      * can be abbreviated as *
   sel hev = (!h;)
 
   select site = ( byres ( resi 45:52 expand 5 ))
      * can be abbreviated as *
   sel site = (b;(i;45:52 x;5))
 
   select combi = ( hev and not site )
      * can be abbreviated as *
   sel combi = (hev&!site)
   '''
### -------------------------------------------------------------------
def selections():
   '''
DESCRIPTION
 
   Selections are surrounded by parenthesis and contain an
   expression consisting of predicates, logical operations, objects,
   named selections and additional parenthesis:
       ( [... [(...) ... ]] )
 
   Properties  
      name <atom names>           n;<atom names>          
      resn <residue names>        r; 
      resi <residue identifiers>  i;<residue identifiers>
      chain <chain ID>            c;<chain identifiers>
      segi <segment identifiers>  s;<segment identifiers>
      elem <element symbol>       e;<element symbols>
      b <operator> <value>        -
   Generic 
      hydro                       h;
      all                         *
   Logical
      and                         &
      or                          |
      not                         !
   Modifiers                        
      byres <selection>           b;<selection>
      around <distance>           a;<distance>
      expand <distance>           e;<distance>
      in <selection>              -
 
   Objects and selections can be referred to by name from within
   subsequent selections.
 
   "help examples" for some example selections
   '''
   help('selections')

def _split(*arg): # custom split-and-trim
   '''
split(string,token[,count]) -> list of strings
 
UTILITY FUNCTION, NOT PART OF THE API
Breaks strings up by tokens but preserves quoted strings and
parenthetical groups (such as atom selections).
'''
   str = arg[0]
   tok = arg[1]
   if len(arg)>2:
      mx=arg[2]
   else:
      mx=0
   pair = { '(':')','[':']','{':'}',"'":"'",'"':'"' }
   plst = pair.keys()
   stack = []
   lst = []
   c = 0
   nf = 0
   l = len(str)
   wd = ""
   while str[c]==tok:
      c = c + 1
   while c<l:
      ch = str[c]
      if (ch in tok) and (len(stack)==0):
         lst.append(string.strip(wd))
         nf = nf + 1
         if mx:
            if nf==mx:
               wd = string.strip(str[c+1:])
               break;
         wd = ''
         w = 0
      else:
         if len(stack):
            if ch==stack[0]:
               stack = stack[1:]
            elif (ch in plst):
               stack[:0]=[pair[ch]]
         elif (ch in plst):
            stack[:0]=[pair[ch]]
         wd = wd + ch
      c = c + 1
   if len(wd):
      lst.append(string.strip(wd))
   return lst
   

def sort(*arg):
   '''
UNSUPPORTED - LIKELY TO BE REMOVED
'''
   lock()
   if len(arg)==0:
      r = _pm.sort("")
   else:
      r = _pm.sort(arg[0])
   unlock()
   return r

def cls():
   lock()
   r = _pm.cls()
   unlock()
   return r
   
def mem():
   '''
DESCRIPTION

"mem" Dumps current memory state to standard output. This is a
debugging feature, not an official part of the API.

'''
   lock()
   r = _pm.mem()
   unlock()
   return r

def dist(*arg):
   '''
DESCRIPTION
 
"dist" creates a new distance object between two
selections.  It will display all distances within a cutoff.
 
USAGE
 
   dist 
   dist (selection1), (selection2)
   dist name = (selection1), (selection1) [,dist [,mode] ]
 
   name = name of distance object 
   selection1,selection2 = atom selections
   cutoff = maximum distance to display
   mode = 0 (default)
 
PYMOL API
 
   pm.dist( string name, string selection1, string selection2,
          string cutoff, string mode )
   returns None
 
NOTES
 
   "dist" alone will show distances between selections (pk1) and (pk3)
   created by left and right button atom picks (hold down CTRL)
'''
   la = len(arg)
   if la<2:
      lock()
      cnt = _pm.get("dist_counter") + 1.0
      _pm.set("dist_counter","%1.0f" % cnt)
      nam = "dist%02.0f" % cnt
      unlock()
      if la==0:
         argst = "(pk1),(pk3)"
      else:
         argst = arg[0]
   else:
      nam = arg[0]
      argst = string.join(arg[1:],',')
   arg = _split(argst,',')
   la = len(arg)
   if la<2:
      print "Error: invalid arguments for dist command."
      raise RuntimeError
   else:
      sel1 = arg[0]
      sel2 = arg[1]
      optarg1=-1.0
      optarg2=0
      if(la>2):
         optarg1 = float(arg[2])
      lock()
      r = _pm.dist(nam,sel1,sel2,optarg2,optarg1)
      unlock()
   return r

def help(*arg):
   '''
USAGE
 
   help <command>
   '''
   if len(arg):
      cmd = arg[0]
   else:
      cmd = 'commands'
   if kwhash.has_key(cmd):
      cc = kwhash[cmd]
      if cc:
         cmd=cc
   if keyword.has_key(cmd):
      doc = keyword[cmd][0].__doc__
      if doc:
         print " \n",string.strip(doc),"\n \n"
      else:
         print "Error: sorry no help available on that command."
   elif help_only.has_key(cmd):
      doc = help_only[cmd][0].__doc__
      if doc:
         print " \n",string.strip(doc),"\n \n"
      else:
         print "Error: sorry no help available on that command."      
   else:
      print "Error: unrecognized command"

def symexp(*arg):
   '''
DESCRIPTION
 
"symexp" creates all symmetry related objects for the specified object
that occurs within a cutoff about an atom selection.  The new objects
are labeled using the prefix provided
 
USAGE
 
   symexp prefix = object, (selection), cutoff
   symexp 
 
PYMOL API
 
   pm.symexp( string prefix, string object, string selection,
           string cutoff) 

   '''
   if len(arg)<2:
      print "Error: invalid arguments for symexp command."
      raise RuntimeError
   elif len(arg)<3:
      nam= arg[0]
      argst = arg[1]
   else:
      nam= arg[0]
      argst = string.join(arg[1:],',')
   arg = _split(argst,',')
   la = len(arg)
   if la<3:
      print "Error: invalid arguments for symexp command."
      raise RuntimeError
   else:
      obj=arg[0]
      sele=arg[1]
      dist=arg[2]
      lock()
      r = _pm.symexp(nam,obj,sele,float(dist))
      unlock()
   return r

def isomesh(nam,argst):
   '''
DESCRIPTION
 
"isomesh" creates a mesh isosurface object from a map object.
 
USAGE
 
   isomesh name = map-object, level [,(selection) [,buffer] ] 
   '''
   arg = _split(argst,',')
   la = len(arg)
   if la<1:
      print "Error: invalid arguments for isomesh command."
      raise RuntimeError
   else:
      map=arg[0]
      mopt=0
      optarg1=''
      optarg2=''
      lvl = 1.0
      if la>1:
         lvl = float(arg[1])
      if la>2:
         if arg[2][0] == '(':
            mopt = 1
            optarg1=arg[2]
      if la>3:
         optarg2 = arg[3]
      lock()
      r = _pm.isomesh(nam,0,map,mopt,optarg1,optarg2,lvl,0)
      unlock()
   return r

def isodot(nam,argst):
   '''
DESCRIPTION
 
"isomesh" creates a dot isosurface object from a map object.
 
USAGE
 
   isodot name = map-object, level [,(selection) [,buffer] ] 
   '''
   arg = _split(argst,',')
   la = len(arg)
   if la<1:
      print "Error: invalid arguments for isodot command."
      raise RuntimeError
   else:
      map=arg[0]
      mopt=0
      optarg1=''
      optarg2=''
      lvl = 1.0
      if la>1:
         lvl = float(arg[1])
      if la>2:
         if arg[2][0] == '(':
            mopt = 1
            optarg1=arg[2]
      if la>3:
         optarg2 = arg[3]
      lock()
      r = _pm.isomesh(nam,0,map,mopt,optarg1,optarg2,lvl,1)
      unlock()
   return r

def ready():
   '''
INTERNAL USAGE
   '''
   return _pm.ready()

def splash():
   '''
DESCRIPTION
 
"splash" shows the splash screen information.
   '''
   lock()
   r = _pm.splash()
   unlock()
   return r

def copy(dst,src):
   '''
DESCRIPTION
 
"copy" creates a new object that is an identical copy of an
existing object
 
USAGE
 
   copy name = object
 
PYMOL API
 
   pm.copy(new-object-name,object)
   '''
   lock()
   r = _pm.copy(src,dst)
   unlock()
   return r

def alter(sele,expr):
   '''
DESCRIPTION
 
"alter" changes one or more atomic properties over a selection using
the python evaluator with a separate name space for each atom.  The
symbols defined in the name space are:
 
   name, resn, resi, chain, x, y, z,
   q, b, segi, and type (0=ATOM,1=HETATM) 
 
All strings in the expression must be explicitly quoted.  This
operation typically takes several seconds per thousand atoms altered.
 
USAGE
 
   alter (selection),expression
 
EXAMPLES
 
   alter (chain A),chain='B'
   alter (all),resi=str(int(resi)+100)
   '''
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
   '''
DESCRIPTION
 
"stereo" activates or deactives stereo mode.  Currently only high-end
stereo graphics are supported on the SGI (stereo in a window), and it
is necessary to launch the program with a "-s" option to activate
this feature.
 
USAGE
 
   stereo on
   stereo off
   '''
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
   
def overlap(*arg):
   '''
UNSUPPORTED FEATURE - LIKELY TO CHANGE
   '''
   if len(arg)==3:
      state[0]=int(arg[2][0])
      if state[0]<1: state[0]=1;
      state[1]=int(arg[2][1])
      if state[1]<1: state[1]=1
   lock()
   r = _pm.overlap(arg[0],arg[1],state[0]-1,state[1]-1)
   unlock()
   return r

def distance(*arg):
   '''
OBSOLETE - TO BE REMOVED
   '''
   la = len(arg)
   if la==0:
      a="pk1"
      b="pk3"
   elif la==1:
      a=arg[0]
      b="pk1"
   elif la==2:
      a=arg[0]
      b=arg[1]
   if a[0]!='(': a="(%"+a+")"
   if b[0]!='(': b="(%"+b+")"
   lock()   
   r = _pm.distance(a,b)
   unlock()
   return r

def _alter_do(at):
   '''
INTERNAL
   '''
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
   '''
INTERNAL
   '''
   __main__.lock_api = _pm.get_globals()['lock_api']
   
def lock():
   '''
INTERNAL
   '''
   __main__.lock_api.acquire(1)

def lock_attempt():
   '''
INTERNAL
   '''
   res = __main__.lock_api.acquire(blocking=0)
   if res:
      _pm.get_globals()['lock_state'] = 1;
   else:
      _pm.get_globals()['lock_state'] = None;

def unlock():
   '''
INTERNAL
   '''
   __main__.lock_api.release()

def export_dots(a,b):
   '''
UNSUPPORTED - WILL BE REMOVED
   '''
   lock()
   r = _pm.export_dots(a,int(b))
   unlock()
   return r

def count_states(*arg):
   '''
UNDOCUMENTED
   '''
   lock()
   if not len(arg):
      a = "(all)"
   else:
      a=arg[0]
   r = _pm.count_states(a)
   unlock()
   return r

def do(a):
   '''
DESCRIPTION
 
   "do" makes it possible for python programs to issue simple PyMOL
   commands as if they were entered on the command line.
    
PYMOL API
 
   pm.do( command )
 
USAGE (PYTHON)
 
   import pm
   pm.do("load file.pdb")
   '''
   lock()
   r = _pm.do(a);
   unlock()
   return r

def turn(a,b):
   '''
DESCRIPTION
  
   "turn" rotates the world about one of the three primary axes
      
USAGE
  
   turn axis, angle
 
EXAMPLES
 
   turn x,90
   turn y,45
    
PYMOL API
 
   pm.turn( string axis, float angle )
   '''
   lock()
   r = _pm.turn(a,float(b))
   unlock()
   return r

def ray():
   '''
DESCRIPTION
  
   "ray" creates a ray traced image of the current frame. This
   can take some time (up to several minutes).
      
USAGE
 
   ray
 
PYMOL API
  
   pm.ray()
   
   '''
   lock()   
   r = _pm.render()
   unlock()
   return r

def real_system(a):
   '''
UNSUPPORTED
   '''
   r = _pm.system(a)
   return r

def system(a):
   '''
UNSUPPORTED
   '''
   real_system(a)

def intra_fit(*arg):
   '''
DESCRIPTION
  
   "intra_fit" fits all states of an object to an atom selection
   in the specified state.  It returns the rms values to python
   as an array.
      
USAGE
 
   intra_fit (selection),state
 
PYMOL API
  
   pm.intra_fit( string selection, int state )
    
EXAMPLES
 
   intra_fit ( name ca )
   
PYTHON EXAMPLE
 
   import pm
   rms = pm.intra_fit("(name ca)",1)
   '''
   lock()
   if len(arg)<2:
      b=-1
   else:
      b=int(arg[1])-1
   r = _pm.intrafit(arg[0],b,2)
   unlock()
   return r

def intra_rms(*arg):
   '''
DESCRIPTION
  
   "intra_rms" calculates rms fit values for all states of an object
   over an atom selection relative to the indicated state.  
   Coordinates are left unchanged.  The rms values are returned
   as a python array.
      
PYMOL API
 
   pm.intra_rms( string selection, int state)
 
PYTHON EXAMPLE
 
   import pm
   rms = pm.intra_rms("(name ca)",1)
   '''
   lock()
   if len(arg)<2:
      b=-1
   else:
      b=int(arg[1])-1
   r = _pm.intrafit(arg[0],b,1)
   unlock()
   return r

def intra_rms_cur(*arg):
   '''
DESCRIPTION
  
   "intra_rms_cur" calculates rms values for all states of an object
   over an atom selection relative to the indicated state without
   performing any fitting.  The rms values are returned
   as a python array.
      
PYMOL API
 
   pm.intra_rms_cur( string selection, int state)
 
PYTHON EXAMPLE
 
   import pm
   rms = pm.intra_rms_cur("(name ca)",1)
   '''
   lock()
   if len(arg)<2:
      b=-1
   else:
      b=int(arg[1])-1
   r = _pm.intrafit(arg[0],b,0)
   unlock()
   return r

def fit(a,b):
   '''
DESCRIPTION
  
   "fit" superimposes the model in the first selection on to the model
   in the second selection.  Only matching atoms in both selections
   will be used for the fit.
   
USAGE
 
   fit (selection), (selection)
 
EXAMPLES
 
   fit ( mutant and name ca ), ( wildtype and name ca )
   '''
   if a[0]!='(': a="(%"+a+")"
   if b[0]!='(': b="(%"+b+")"
   lock()   
   r = _pm.fit("(%s in %s)" % (a,b),
         "(%s in %s)" % (b,a),2)
   unlock()
   return r

def rms(a,b):
   '''
DESCRIPTION
  
   "rms" computes a RMS fit between two atom selections, but does not
   tranform the models after performing the fit.
   
USAGE
 
   rms (selection), (selection)
 
EXAMPLES
 
   fit ( mutant and name ca ), ( wildtype and name ca )
   '''
   if a[0]!='(': a="(%"+a+")"
   if b[0]!='(': b="(%"+b+")"
   lock()   
   r = _pm.fit("(%s in %s)" % (a,b),
         "(%s in %s)" % (b,a),0)
   unlock()
   return r

def rms_cur(a,b):
   '''
DESCRIPTION
  
   "rms_cur" computes the RMS difference between two atom
   selections without performing any fitting.
   
USAGE
 
   rms_cur (selection), (selection)
   '''
   if a[0]!='(': a="(%"+a+")"
   if b[0]!='(': b="(%"+b+")"
   lock()   
   r = _pm.fit("(%s in %s)" % (a,b),
         "(%s in %s)" % (b,a),1)
   unlock()
   return r

def pairfit(*arg):
   '''
DESCRIPTION
  
   "pair_fit" fits a set of atom pairs between two models.  Each atom
   in each pair must be specified individually, which can be tedious
   to enter manually.  Script files are recommended when using this
   command.
   
USAGE
 
   fit_pairs (selection), (selection) [ (selection), (selection) [ ...] ]
   '''
   lock()   
   r = _pm.fit_pairs(arg)
   unlock()
   return r

def expfit(a,b):
   '''
   ??? OBSOLETE
   '''
   lock()   
   r = _pm.fit(a,b,2)
   unlock()
   return r
   
def zoom(a):
   '''
DESCRIPTION
  
   "zoom" scales and translates the window and the origin to cover the
   atom selection.
      
USAGE
 
   zoom object-or-selection
   zoom (selection)
 
PYMOL API
 
   pm.orient( string object-or-selection )
   '''
   lock()   
   r = _pm.zoom(a)
   unlock()
   return r
   
def frame(a):
   '''
DESCRIPTION
  
   "frame" sets the viewer to the indicated movie frame.
   
USAGE
 
   frame frame-number
 
PYMOL API
 
   pm.frame( int frame_number )
 
NOTES
 
   Frame numbers are 1-based 
   '''
   lock()   
   r = _pm.frame(int(a))
   unlock()
   return r

def move(a,b):
   '''
DESCRIPTION
  
   "move" translates the world about one of the three primary axes.
      
USAGE
  
   move axis,angle
 
EXAMPLES
 
   move x,3
   move y,-1
    
PYMOL API
 
   pm.move( string axis, float distance )
   '''
   lock()   
   r = _pm.move(a,float(b))
   unlock()
   return r

def clip(a,b):
   '''
DESCRIPTION
  
   "clip" translates the near and far clipping planes
      
USAGE
  
   clip {near|far}, distance
 
EXAMPLES
 
   clip near, -5
   clip far, 10
    
PYMOL API
 
   pm.clip( string plane, float distance )
   '''
   lock()   
   r = _pm.clip(a,float(b))
   unlock()
   return r

def origin(a):
   '''
DESCRIPTION
  
   "origin" sets the center of rotation about a selection
      
USAGE
 
   origin object-or-selection
   origin (selection)
 
PYMOL API
 
   pm.origin( string object-or-selection )
   '''
   lock()   
   r = _pm.origin(a)
   unlock()
   return r

def orient(*arg):
   '''
DESCRIPTION
  
   "orient" aligns the principal components of the atoms in the
   selection with the XYZ axes.  The function is similar to the
   orient command in X-PLOR.
      
USAGE
 
   orient object-or-selection
   orient (selection)
 
PYMOL API
 
   pm.orient( string object-or-selection )
   '''
   lock()
   if len(arg)<1:
      a = "(all)"
   else:
      a = arg[0]
   r = _pm.orient(a)
   unlock()
   return r

def is_glut_thread():
   '''
INTERNAL
   '''
   if thread.get_ident() == __main__.glutThread:
      return 1
   else:
      return 0

def refresh():
   '''
INTERNAL
   '''
   lock()
   if thread.get_ident() == __main__.glutThread:
      r = _pm.refresh_now()
   else:
      r = _pm.refresh()
   unlock()
   return r

def dirty():
   '''
INTERNAL
   '''
   lock()   
   r = _pm.dirty()
   unlock()
   return r

def set(a,b):
   '''
DESCRIPTION
  
   "set" changes one of the PyMOL state variables
      
USAGE
 
   set variable = value
 
PYMOL API
 
   pm.set ( string variable, string value )
   '''
   lock()   
   r = _pm.set(a,b)
   unlock()
   return r

def reset():
   '''
DESCRIPTION
  
   "reset" restores the rotation matrix to identity, sets the origin
   to the center of mass (approx.) and zooms the window and clipping
   planes to cover all objects.
   
USAGE
 
   reset 
 
PYMOL API
 
   pm.reset ( )
   '''
   lock()   
   r = _pm.reset(0)
   unlock()
   return r

def meter_reset():
   '''
UNDOCUMENTED
   '''
   lock()   
   r = _pm.reset_rate()
   unlock()
   return r

def delete(a):
   '''
DESCRIPTION
  
   "delete" removes an object or a selection. 
   
USAGE
 
   delete object-or-selection-name
   delete all
 
PYMOL API
 
   pm.delete ( string object-or-selection-name )
   '''
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
   '''
DESCRIPTION
  
   "quit" terminates the program. 
   
USAGE
 
   quit
 
PYMOL API
 
   pm.quit()
   '''
   if thread.get_ident() ==__main__.glutThread:
      lock()
      r = _pm.do("_quit")
   else:
      r = _pm.do("_quit")
      thread.exit()
   return r

def png(a):
   '''
DESCRIPTION
  
   "png" writes a png format image file of the current image to disk.
   
USAGE
 
   png filename
 
PYMOL API
 
   pm.png( string filename )
   '''
   if thread.get_ident() ==__main__.glutThread:
      r = _png(a)
   else:
      r = _pm.do("pm._png('"+a+"')")
   return r

def _png(a):
   lock()   
   fname = a
   if not re.search("\.png$",fname):
      fname = fname +".png"
   r = _pm.png(fname)
   unlock()
   return r

def mclear():
   '''
DESCRIPTION
  
   "mclear" clears the movie frame image cache.
   
USAGE
 
   mclear
 
PYMOL API
 
   pm.mclear()
   '''
   lock()   
   r = _pm.mclear()
   unlock()
   return r

def _special(k,x,y):
   k=int(k)
   if special.has_key(k):
      if special[k][1]:
         if special[k][2]:
            apply(special[k][1],special[k][3])
         else:
            apply(special[k][1],())
   return None

def set_key(*arg):  
   '''
DESCRIPTION
  
   "set_key" binds a specific python function to a key press.
   
PYMOL API
 
   pm.set_key( string key, function fn, tuple arguments)
 
PYTHON EXAMPLE
 
   import pm
 
   def color_blue(object):
      pm.color("blue",object)
    
   pm.set_key( 'F1' , make_it_blue, ( "object1" ) )
   pm.set_key( 'F2' , make_it_blue, ( "object2" ) )
 
   // would turn object1 blue when the F1 key is pressed and
   // would turn object2 blue when the F2 key is pressed.
   '''
   key=arg[0]
   cmd=arg[1]
   if len(arg)>2:
      cmd_arg=arg[2]
   else:
      cmd_arg=None 
   for a in special.keys():
      if special[a][0]==key:
         special[a][1]=cmd
         if cmd_arg:
            special[a][2]=1
            special[a][3]=cmd_arg
         else:
            special[a][2]=0
            special[a][3]=None

def mstop():
   '''
DESCRIPTION
  
   "mstop" stops the movie.
   
USAGE
 
   mstop
 
PYMOL API
 
   pm.mstop()
   '''
   lock()   
   r = _pm.mplay(0)
   unlock()
   return r

def mplay():
   '''
DESCRIPTION
  
   "mplay" starts the movie.
   
USAGE
 
   mplay
 
PYMOL API
 
   pm.mplay()
   '''
   lock()   
   r = _pm.mplay(1)
   unlock()
   return r

def mray():
   '''
DEPRECATED
   '''
   lock()   
   r = _pm.mplay(2)
   unlock()
   return r

def viewport(a,b):
   '''
DESCRIPTION
  
   "viewport" changes the size of the viewing port (and thus the size
   of all png files subsequently output)
      
USAGE
 
   viewport width, height
 
PYMOL API
  
   pm.viewport(int width, int height)
   '''
   r = _pm.viewport(int(a),int(b))
   
def mdo(a,b):
   '''
DESCRIPTION
  
   "mdo" sets up a command to be executed upon entry into the
   specified frame of the movie.  These commands are usually created
   by a PyMOL utility program (such as pmu.mrock).  Command can
   actually contain several commands separated by semicolons ';'
 
USAGE
 
   mdo frame : command
 
PYMOL API
  
   pm.mdo( int frame, string command )
 
EXAMPLE
 
   // Creates a single frame movie involving a rotation about X and Y
   
   load test.pdb
   mset 1
   mdo 1: turn x,5; turn y,5;
   mplay
   
NOTES
 
   The "mset" command must first be used to define the movie before
   "mdo" statements will have any effect.  Redefinition of the movie
   clears any existing mdo statements.
   '''
   lock()   
   r = _pm.mdo(int(a)-1,b)
   unlock()
   return r

def dummy(*arg):
   '''
DEBUGGING
   '''
   lock()   
   pass
   unlock()
   return None

def rock():
   '''
DESCRIPTION
  
   "rock" toggles Y axis rocking.
 
USAGE
 
   rock
 
PYMOL API
  
   pm.rock()
   '''
   lock()   
   r = _pm.rock()
   unlock()
   return r

def forward():
   '''
DESCRIPTION
  
   "forward" moves the movie one frame forward.
 
USAGE
 
   forward
 
PYMOL API
  
   pm.forward()
   '''
   lock()   
   r = _pm.setframe(5,1)
   unlock()
   return r

def backward():
   '''
DESCRIPTION
  
   "backward" moves the movie back one frame.
 
USAGE
 
   backward
 
PYMOL API
  
   pm.backward()
   '''
   lock()   
   r = _pm.setframe(5,-1)
   unlock()
   return r


def rewind():
   '''
DESCRIPTION
  
   "rewind" goes to the beginning of the movie.
 
USAGE
 
   rewind
 
PYMOL API
  
   pm.rewind()
   '''
   lock()   
   r = _pm.setframe(0,0)
   unlock()
   return r

def ending():
   '''
DESCRIPTION
  
   "ending" goes to the end of the movie.
 
USAGE
 
   ending
 
PYMOL API
  
   pm.ending()
   '''
   lock()   
   r=_pm.setframe(2,0)
   unlock()
   return r

def middle():
   '''
DESCRIPTION
  
   "middle" goes to the middle of the movie.
 
USAGE
 
   middle
 
PYMOL API
  
   pm.middle()
   '''
   lock()   
   r = _pm.setframe(3,0)
   unlock()
   return r

def test(): # generic test routine for development
   '''
DEBUGGING
   '''
   lock()   
   r=_pm.test()
   unlock()
   return r

def dump(fnam,obj):
   '''
DEBUGGING
   '''
   lock()
   r = _pm.dump(fnam,obj)
   unlock()
   return r

def save(*arg):
   '''
DESCRIPTION
  
   "save" writes selected atoms to a PDB file
 
USAGE
 
   save filename [,(selection) [,state] ]
 
PYMOL API
  
   pm.save(filename, selection, state)
   '''
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

def load(*arg):
   '''
DESCRIPTION
  
   "load" reads several file formats.  The file extension is used to
   determine the format.  PDB files must end in ".pdb", MOL files must
   end in ".mol", Macromodel files must end in ".mmod".  and XPLOR
   maps must end in ".xplor".
 
   If an object is specified, then the file is load into that object.
   Otherwise, an object is created with the same name as the file
   prefix.
 
USAGE
 
   load filname [,object ,[state]]
 
PYMOL API
  
   pm.load( filename [,object [,state]] )
   '''
   r = 1
   lock()   
   ftype = 0
   state = -1
   if re.search("\.pdb$",arg[0]):
      ftype = 0
   elif re.search("\.mol$",arg[0]):
      ftype = 1
   elif re.search("\.mmod$",arg[0]):
      ftype = 4
   elif re.search("\.xplor$",arg[0]):
      ftype = 7
   if len(arg)==1:
      oname = re.sub("[^/]*\/","",arg[0])
      oname = re.sub("\.pdb|\.mol|\.mmod|\.xplor","",oname)
      r = _pm.load(oname,arg[0],state,ftype)
   elif len(arg)==2:
      oname = string.strip(arg[1])
      r = _pm.load(oname,arg[0],state,ftype)
   elif len(arg)==3:
      oname = string.strip(arg[1])
      state = int(arg[2])-1
      r = _pm.load(oname,arg[0],state,ftype)
   elif len(arg)==4:
      if loadable.has_key(arg[3]):
         ftype = loadable[arg[3]]
      else:
         ftype = int(arg[3])
      state = int(arg[2])-1
      oname = string.strip(arg[1])
      r = _pm.load(oname,arg[0],state,ftype)
   else:
      print "argument error."
   unlock()
   return r

def read_molstr(*arg):
   r = 1
   lock()   
   ftype = 3
   if len(arg)==2:
      oname = string.strip(arg[1])
      r = _pm.load(oname,arg[0],-1,ftype)
   elif len(arg)==3:
      oname = string.strip(arg[1])
      r = _pm.load(oname,arg[0],int(arg[2])-1,ftype)
   else:
      print "argument error."
   unlock()
   return r

def read_mmodstr(*arg):
   r = 1
   lock()   
   ftype = 6
   if len(arg)==2:
      oname = string.strip(arg[1])
      r = _pm.load(oname,arg[0],-1,ftype)
   elif len(arg)==3:
      oname = string.strip(arg[1])
      r = _pm.load(oname,arg[0],int(arg[2])-1,ftype)
   else:
      print "argument error."
   unlock()
   return r
   
def select(*arg):
   '''
DESCRIPTION
  
   "select" creates a named selection from an atom selection.
 
USAGE
 
   select (selection)
   select selection-name = 
 
PYMOL API
  
   pm.load( filename [,object [,state]] )
 
EXAMPLES 
 
   select near = (pk1 expand 8)
   select bb = (name ca,n,c,o )
   '''
   lock()   
   if len(arg)==1:
      sel_cnt = _pm.get("sel_counter") + 1.0
      _pm.set("sel_counter","%1.0f" % sel_cnt)
      sel_name = "sel%02.0f" % sel_cnt
      sel = arg[0]
   else:
      sel_name = arg[0]
      sel = arg[1]
   r = _pm.select(sel_name,sel)
   unlock()
   return r

def color(*arg):
   '''
DESCRIPTION
  
   "color" changes the color of an object or an atom selection.
   
USAGE
 
   color color-name
   color color-name, object-name
   color color-name, (selection)
 
PYMOL API
  
   pm.color( string color, string color-name )
 
EXAMPLES 
 
   color yellow, (name C*)
   '''
   lock()   
   if len(arg)==2:
      r = _pm.color(arg[0],arg[1],0)
   else:
      r = _pm.color(arg[0],"(all)",0)   
   unlock()
   return r

def colordef(nam,col):
   '''
WARNING - THIS FUNCTION WILL BE SOON OBSOLETED

DESCRIPTION
  
   "colordef" defines a new color with color indices (0.0-1.0)
   
USAGE
 
   colordef color, red-float green-float blue-float
 
PYMOL API
  
   pm.colordef( string color, string color-components )
 
EXAMPLES 
 
   colordef red, 1.0 0.0 0.0 
   '''
   r = 1
   lock()
   c = string.split(col)
   if len(c)==3:
      r = _pm.colordef(nam,float(c[0]),float(c[1]),float(c[2]))
   else:
      print "Error: invalid color vector."
   unlock()
   return r

def mpng(a):
   '''
DESCRIPTION
  
   "mpng" writes a series of numbered movie frames to png files with
   the specified prefix.  If the "ray_trace_frames" variable is
   non-zero, these frames will be ray-traced.  This operation can take
   several hours for a long movie.
 
   Be sure to disable "cache_frames" when issuing this operation on a
   long movie (typically >100 frames to avoid running out of memory).
   
USAGE
 
   mpng prefix
 
PYMOL API
 
   pm.mpng( string prefix )
   '''
   if thread.get_ident() ==__main__.glutThread:
      r = _mpng(a)
   else:
      r = _pm.do("pm._mpng('"+a+"')")
   return r

def _mpng(*arg):
   lock()   
   fname = arg[0]
   if re.search("\.png$",fname):
      fname = re.sub("\.png$","",fname)
   r = _pm.mpng_(fname)
   unlock()
   return r

def show(*arg):
   '''
DESCRIPTION
  
   "show" turns on atom and bond representations.
 
   The available representation are:
    
      lines     spheres   mesh      ribbon
      sticks    dots      surface     
   
USAGE
 
   show
   show reprentation [,object]
   show reprentation [,(selection)]
 
PYMOL API
 
   pm.show( string representation, string object-or-selection )
 
EXAMPLES
 
   show lines,(name ca or name c or name n)
   show ribbon
 
NOTES
 
   "show" alone will turn on lines for all bonds.
   '''
   r=1
   lock()
   l = len(arg)
   if not l:
      r = _pm.showhide("(all)",0,1); # show lines by default       
   elif l==2:
      rep = arg[0]
      if rephash.has_key(rep):
         rep = rephash[rep]
      if repres.has_key(rep):      
         repn = repres[rep];
         r = _pm.showhide(arg[1],repn,1);
      else:
         print "Error: unrecognized or ambiguous representation"
   elif arg[0]=='all':
      r = _pm.showhide("(all)",0,1); # show lines by default 
   elif arg[0][0]=='(':
      r = _pm.showhide(arg[0],0,1);
   else:
      rep = arg[0]
      if rephash.has_key(rep):
         rep = rephash[rep]
      if repres.has_key(rep):      
         repn = repres[rep];
         r = _pm.showhide("(all)",repn,1);
      else:
         print "Error: unrecognized or ambiguous representation"
   unlock()
   return r

def hide(*arg):
   '''
DESCRIPTION
  
   "hide" turns of atom and bond representations.
 
   The available representation are:
    
      lines     spheres   mesh      ribbon
      sticks    dots      surface     
   
USAGE
 
   hide reprentation [,object]
   hide reprentation [,(selection)]
 
PYMOL API
 
   pm.hide( string representation, string object-or-selection )
 
EXAMPLES
 
   hide lines,all
   hide ribbon
   '''
   r = 1
   l = len(arg)
   lock()
   if not l:
      r = _pm.showhide("!",0,0);      
   elif l==2:
      rep = arg[0]
      if rephash.has_key(rep):
         rep = rephash[rep]
      if repres.has_key(rep):      
         repn = repres[rep];
         r = _pm.showhide(arg[1],repn,0);
      else:
         print "Error: unrecognized or ambiguous representation"
   elif arg[0]=='all':
      r = _pm.showhide("!",0,0);
   elif arg[0][0]=='(':
      r = _pm.showhide(arg[0],-1,0);
   else:
      rep = arg[0]
      if rephash.has_key(rep):
         rep = rephash[rep]
      if repres.has_key(rep):
         repn = repres[arg[0]];
         r = _pm.showhide("(all)",repn,0);
      else:
         print "Error: unrecognized or ambiguous representation"
   unlock()
   return r

def mmatrix(a):
   '''
DESCRIPTION
  
   "mmatrix" sets up a matrix to be used for the first frame of the movie.
   
USAGE
 
   mmatrix {clear|store|recall}
 
PYMOL API
 
   pm.mmatrix( string action )
 
EXAMPLES
 
   mmatrix store
   '''
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
   '''
DESCRIPTION
  
   "enable" enable display of an object and all currently visible representations.
   
USAGE
 
   enable object
   enable all
   
PYMOL API
 
   pm.enable( string object-name )
 
EXAMPLE
 
   enable my_object
   '''
   lock()   
   r = _pm.onoff(nam,1);
   unlock()
   return r

def disable(nam):
   '''
DESCRIPTION
  
   "disable" disables display of an object and all currently visible representations.
   
USAGE
 
   disable object
   disable all 
 
PYMOL API
 
   pm.disable( string object-name )
 
EXAMPLE
 
   disable my_object
   '''
   lock()   
   r = _pm.onoff(nam,0);
   unlock()
   return r

def mset(seq):
   '''
DESCRIPTION
  
   "mset" sets up a relationship between molecular states and movie
   frames.  This makes it possible to control which states are shown
   in which frame.
   
USAGE

   mset specification
 
PYMOL API
 
   pm.mset( string specification )
 
EXAMPLES

   mset 1         // simplest case, one state -> one frame
   mset 1 x10     // ten frames, all corresponding to state 1
   mset 1 x30 1 -15 15 x30 15 -1
     // more realistic example:
     // the first thirty frames are state 1
     // the next 15 frames pass through states 1-15
     // the next 30 frames are of state 15
     // the next 15 frames iterate back to state 1
   '''
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

def null():
   pass

keyword = { 
   'alter'         : [alter        , 2 , 2 , ',' , 0 ],
   'api'           : [api          , 0 , 0 , ',' , 0 ],
   'backward'      : [backward     , 0 , 0 , ',' , 0 ],
   'clip'          : [clip         , 2 , 2 , ',' , 0 ],
   'cls'           : [cls          , 0 , 0 , ',' , 0 ],
   'color'         : [color        , 1 , 2 , ',' , 0 ],
   'colordef'      : [colordef     , 2 , 2 , ',' , 0 ],
   'commands'      : [commands     , 0 , 0 , ',' , 0 ],
   'copy'          : [copy         , 2 , 2 , '=' , 0 ],
   'count_states'  : [count_states , 0 , 1 , ',' , 0 ],   
   'delete'        : [delete       , 1 , 1 , ',' , 0 ],
   'disable'       : [disable      , 1 , 1 , ',' , 0 ],
   'dist'          : [dist         , 0 , 2 , '=' , 0 ],
   'distance'      : [distance     , 0 , 2 , '=' , 0 ],
   'dump'          : [dump         , 2 , 2 , ',' , 0 ],
   'enable'        : [enable       , 1 , 1 , ',' , 0 ],
   'ending'        : [ending       , 0 , 0 , ',' , 0 ],
   'export_dots'   : [export_dots  , 2 , 2 , ',' , 0 ],
   'fit'           : [fit          , 2 , 2 , ',' , 0 ],
   'fork'          : [dummy        , 1 , 1 , ',' , 3 ],
   'forward'       : [forward      , 0 , 0 , ',' , 0 ],
   'frame'         : [frame        , 1 , 1 , ',' , 0 ],
   'help'          : [help         , 0 , 1 , ',' , 0 ],
   'hide'          : [hide         , 0 , 2 , ',' , 0 ],
   'intra_fit'     : [intra_fit    , 1 , 2 , ',' , 0 ],
   'intra_rms'     : [intra_rms    , 1 , 2 , ',' , 0 ],
   'intra_rms_cur' : [intra_rms_cur, 1 , 2 , ',' , 0 ],
   'isodot'        : [isodot       , 2 , 2 , '=' , 0 ],   
   'isomesh'       : [isomesh      , 2 , 2 , '=' , 0 ],
   'load'          : [load         , 1 , 4 , ',' , 0 ],
   'mem'           : [mem          , 0 , 0 , ',' , 0 ],
   'meter_reset'   : [meter_reset  , 0 , 0 , ',' , 0 ],
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
   'ray'           : [ray          , 0 , 0 , ',' , 0 ],
   'refresh'       : [refresh      , 0 , 0 , ',' , 0 ],
   'reset'         : [reset        , 0 , 0 , ',' , 0 ],
   'rewind'        : [rewind       , 0 , 0 , ',' , 0 ],
   'rock'          : [rock         , 0 , 0 , ',' , 0 ],
   'run'           : [dummy        , 1 , 2 , ',' , 2 ],
   'rms'           : [rms          , 2 , 2 , ',' , 0 ],
   'rms_cur'       : [rms_cur      , 2 , 2 , ',' , 0 ],
   'save'          : [save         , 0 , 4 , ',' , 0 ],
   'select'        : [select       , 1 , 2 , '=' , 0 ],
   'set'           : [set          , 2 , 2 , '=' , 0 ],
   'set_key'       : [set_key      , 2 , 1 , ',' , 0 ], # API only
   'show'          : [show         , 0 , 2 , ',' , 0 ],
   'sort'          : [sort         , 0 , 1 , ',' , 0 ],
   'splash'       : [splash       , 0 , 0 , ',' , 0 ],
   '_special'      : [_special     , 3 , 3 , ',' , 0 ],
   'stereo'        : [stereo       , 1 , 1 , ',' , 0 ],
   'symexp'        : [symexp       , 2 , 2 , '=' , 0 ],
   'system'        : [system       , 1 , 1 , ',' , 0 ],
   'test'          : [test         , 0 , 0 , ',' , 0 ],
   'turn'          : [turn         , 2 , 2 , ',' , 0 ],
   'quit'          : [quit         , 0 , 0 , ',' , 0 ],
   '_quit'         : [_quit        , 0 , 0 , ',' , 0 ],
   'png'           : [png          , 1 , 1 , ',' , 0 ],
   'viewport'      : [viewport     , 2 , 2 , ',' , 0 ],
   'zoom'          : [zoom         , 1 , 1 , ',' , 0 ]
   }

help_only = {
   'selections'    : [selections   , 0 , 0 , ',' , 0 ],
   'keyboard'      : [keyboard     , 0 , 0 , ',' , 0 ],
   'mouse'         : [mouse        , 0 , 0 , ',' , 0 ],
   'examples'      : [examples     , 0 , 0 , ',' , 0 ],
}

repres = {
   'lines'         : 0,
   'sticks'        : 1,
   'dots'          : 2,
   'mesh'          : 3,
   'spheres'       : 4,
   'ribbon'        : 5,
   'surface'       : 6,
   'dashes'        : 7,
   'labels'        : 8
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
   106      : [ 'home'      , rewind                 , 0 , None ],
   107      : [ 'end'       , ending                 , 0 , None ],
   108      : [ 'insert'    , rock                   , 0 , None ]   
}

loadable = {
   'pdb'   : 0,
   'mol'   : 1,
   'mmod'  : 4,
   'xplor' : 7
}

# build shortcuts list

kwhash = {}

for a in keyword.keys():
   for b in range(1,len(a)):
      sub = a[0:b]
      if kwhash.has_key(sub):
         kwhash[sub]=0
      else:
         kwhash[sub]=a

for a in keyword.keys():
   kwhash[a]=a

rephash = {}

for a in repres.keys():
   for b in range(1,len(a)):
      sub = a[0:b]
      if rephash.has_key(sub):
         rephash[sub]=0
      else:
         rephash[sub]=a

for a in repres.keys():
   rephash[a]=a

