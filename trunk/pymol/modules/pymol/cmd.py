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

# cmd.py 
# Python interface module for PyMol
#
# **This is the only module which should be/need be imported by 
# ** PyMol API Based Programs

import re
import _cmd
import string
import traceback
import thread
import types
import pymol
import os

lock_api = pymol.lock_api

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
   COLORS        color    set_color
   HELP          help     commands
   DISTANCES     dist      
   STEREO        stereo
   SYMMETRY      symexp
 
Try "help <command-name>" for more information on a given command.
 
Additional help topics include:
   "movies", "keyboard", "mouse", "selections",
   "examples", "launching", and "api".
   '''
   help('commands')

def api():
   '''
DESCRIPTION
 
The PyMOL Python Application Programming Interface (API) should be
accessed exclusively through the "pm" module.  All command-line
functions have an equivalent API method with the same name.
 
USAGE
 
   from pymol import cmd
   <result> = cmd.<methods>( <args> ) 
    
API-ONLY METHODS
 
   KEY BINDING   set_key
   (more documentation to come)
 
NOTES
 
   Although the PyMOL core is not multi-threaded, the API is
   threads-safe and can be called asynchronously by external python
   programs.  PyMOL handles the necessary locking to insure that
   internal states do not get corrupted.
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
 
   The configuration can be changed by setting the "button_mode"
   variable.  The current configuration is described on screen with
   the following abbreviations:
 
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
   help('examples')
   
def launching():
   '''
PyMOL COMMAND LINE OPTIONS 
 
   pymol.com [-cistwx] <file> ...
 
   -c   Command line mode, no GUI at all.
  
   -i   Disable the internal OpenGL GUI (object list, menus, etc.)
 
   -s   Enable stereo mode (not currently autodetected).
  
   -t   Use Tcl/Tk based external GUI module (pmg_tk).
   -w   Use wxPython based external GUI module (pmg_wx).
 
   -x   Disable the external GUI module.
 
   <file> can have extension:
    
      .pml    PyMOL command script to be run on startup
      .py     Python program to be run on startup
       
      .pdb    Protein Data Bank format file to be loaded on startup
      .mmod   Macromodel format to be loaded on startup
      .mol    MDL MOL file to be loaded on startup
      .xplor  X-PLOR Map file to be loaded on startup
   '''
   help('launching')

def movies():
   '''
MOVIES
 
   To create a movie, simply load multiple coordinate files
   into the same object.  This can be accomplish at the command line,
   using script files, or by writing PyMOL API-based programs.
 
   The commands:
   
load frame001.pdb,mov
load frame002.pdb,mov
 
   will create a two frame movie.  So will the following program:
   
from pymol import cmd
   
for a in ( "frame001.pdb","frame002.pdb" ):
   cmd.load(a,"mov")
 
   which can be executed at the command line using the "run" command.
   
   Python built-it glob module can be useful for loading movies.

from pymol import cmd
import glob
for a in ( glob.glob("frame*.pdb") ):
   cmd.load(a,"mov")
 
NOTE
 
   Because PyMOL stores all movie frames in memory, there is a
   a practical limit to the number of atoms in all coordinate files. 
   160 MB free RAM enables 500,000 atoms with line representations.
   Complex representations require significantly more memory.
   '''
   help('movies')
   
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
      flag <number>               f;
      alt <code>                  - 
   Generic 
      hydro                       h;
      all                         *
      visible                     v;
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
   try:
      lock()
      if len(arg)==0:
         r = _cmd.sort("")
      else:
         r = _cmd.sort(arg[0])
   finally:
      unlock()
   return r

def cls():
   try:
      lock()
      r = _cmd.cls()
   finally:
      unlock()
   return r
   
def mem():
   '''
DESCRIPTION

"mem" Dumps current memory state to standard output. This is a
debugging feature, not an official part of the API.

'''
   try:
      lock()
      r = _cmd.mem()
   finally:
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
 
   cmd.dist( string name, string selection1, string selection2,
          string cutoff, string mode )
   returns None
 
NOTES
 
   "dist" alone will show distances between selections (pk1) and (pk3)
   created by left and right button atom picks (hold down CTRL)
'''
   la = len(arg)
   if la<2:
      try:
         lock()
         cnt = _cmd.get("dist_counter") + 1.0
         _cmd.set("dist_counter","%1.0f" % cnt)
         nam = "dist%02.0f" % cnt
      finally:
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
      try:
         lock()
         r = _cmd.dist(nam,sel1,sel2,optarg2,optarg1)
      finally:
         unlock()
   return r

def show_help(cmd):
   set("text","1")
   print "PyMOL>help %s\n" % cmd
   help(cmd)
   print "(Hit TAB to hide)"
      
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
 
   cmd.symexp( string prefix, string object, string selection,
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
      try:
         lock()
         r = _cmd.symexp(nam,obj,sele,float(dist))
      finally:
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
      try:
         lock()
         r = _cmd.isomesh(nam,0,map,mopt,optarg1,optarg2,lvl,0)
      finally:
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
      try:
         lock()
         r = _cmd.isomesh(nam,0,map,mopt,optarg1,optarg2,lvl,1)
      finally:
         unlock()
   return r

def ready():
   '''
INTERNAL USAGE
   '''
   return _cmd.ready()

def splash():
   '''
DESCRIPTION
 
"splash" shows the splash screen information.
   '''
   set("text","1")
   try:
      lock()
      r = _cmd.splash()
   finally:
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
 
   cmd.copy(new-object-name,object)
   '''
   try:
      lock()
      r = _cmd.copy(src,dst)
   finally:
      unlock()
   return r

def alter(sele,expr):
   '''
DESCRIPTION
 
"alter" changes one or more atomic properties over a selection using
the python evaluator with a separate name space for each atom.  The
symbols defined in the name space are:
 
   name, resn, resi, chain,
   q, b, segi, and type (ATOM,HETATM) 
 
All strings in the expression must be explicitly quoted.  This
operation typically takes several seconds per thousand atoms altered.
 
USAGE
 
   alter (selection),expression
 
EXAMPLES
 
   alter (chain A),chain='B'
   alter (all),resi=str(int(resi)+100)
   '''
   try:
      lock()
      r = _cmd.alter(sele,expr)
   finally:
      unlock()   
   return r

def alter_state(state,sele,expr):
   '''
DESCRIPTION
 
"alter_state" changes the atomic coordinates of a particular state
using the python evaluator with a separate name space for each atom.
The symbols defined in the name space are:
 
   x,y,z
 
USAGE
 
   alter_state state,(selection),expression
 
EXAMPLES
 
   alter 1,(all),x=x+5
   '''
   try:
      lock()
      r = _cmd.alter_state(int(state)-1,sele,expr)
   finally:
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
is necessary to launching the program with a "-s" option to activate
this feature.
 
USAGE
 
   stereo on
   stereo off
   '''
   r = None
   if a=="on":
      try:
         lock()
         if _cmd.stereo(1):
            r = _stereo(1)
         else:
            print " error: stereo not available"
      finally:
         unlock();
   else:
      try:
         lock()
         if _cmd.stereo(0):
            r = _stereo(0)
      finally:
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
   try:
      lock()
      r = _cmd.overlap(arg[0],arg[1],state[0]-1,state[1]-1)
   finally:
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
   try:
      lock()   
      r = _cmd.distance(a,b)
   finally:
      unlock()
   return r

def setup_global_locks():
   '''
INTERNAL
   '''
   pass
#   pymol.lock_api = _cmd.get_globals()['lock_api']

   
def lock():
   '''
INTERNAL
   '''
   lock_api.acquire(1)

def lock_attempt():
   '''
INTERNAL
   '''
   res = lock_api.acquire(blocking=0)
   if res:
      pymol.lock_state = 1;
   else:
      pymol.lock_state = None;

def unlock():
   '''
INTERNAL
   '''
   lock_api.release()

def export_dots(a,b):
   '''
UNSUPPORTED - WILL BE REMOVED
   '''
   try:
      lock()
      r = _cmd.export_dots(a,int(b))
   finally:
      unlock()
   return r

def count_states(*arg):
   '''
UNDOCUMENTED
   '''
   try:
      lock()
      if not len(arg):
         a = "(all)"
      else:
         a=arg[0]
      r = _cmd.count_states(a)
   finally:
      unlock()
   return r

def do(a):
   '''
DESCRIPTION
 
   "do" makes it possible for python programs to issue simple PyMOL
   commands as if they were entered on the command line.
    
PYMOL API
 
   cmd.do( command )
 
USAGE (PYTHON)
 
   from pymol import cmd
   cmd.do("load file.pdb")
   '''
   try:
      lock()
      r = _cmd.do(a);
   finally:
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
 
   cmd.turn( string axis, float angle )
   '''
   try:
      lock()
      r = _cmd.turn(a,float(b))
   finally:
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
  
   cmd.ray()
   
   '''
   try:
      lock()   
      r = _cmd.render()
   finally:
      unlock()
   return r

def real_system(a):
   '''
UNSUPPORTED
   '''
   r = _cmd.system(a)
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
  
   cmd.intra_fit( string selection, int state )
    
EXAMPLES
 
   intra_fit ( name ca )
   
PYTHON EXAMPLE
 
   from pymol import cmd
   rms = cmd.intra_fit("(name ca)",1)
   '''
   if len(arg)<2:
      b=-1
   else:
      b=int(arg[1])-1
   try:
      lock()
      r = _cmd.intrafit(arg[0],b,2)
   finally:
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
 
   cmd.intra_rms( string selection, int state)
 
PYTHON EXAMPLE
 
   from pymol import cmd
   rms = cmd.intra_rms("(name ca)",1)
   '''
   if len(arg)<2:
      b=-1
   else:
      b=int(arg[1])-1
   try:
      lock()
      r = _cmd.intrafit(arg[0],b,1)
   finally:
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
 
   cmd.intra_rms_cur( string selection, int state)
 
PYTHON EXAMPLE
 
   from pymol import cmd
   rms = cmd.intra_rms_cur("(name ca)",1)
   '''
   if len(arg)<2:
      b=-1
   else:
      b=int(arg[1])-1
   try:
      lock()
      r = _cmd.intrafit(arg[0],b,0)
   finally:
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
   try:
      lock()   
      r = _cmd.fit("(%s in %s)" % (a,b),
                  "(%s in %s)" % (b,a),2)
   finally:
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
   try:
      lock()   
      r = _cmd.fit("(%s in %s)" % (a,b),
                  "(%s in %s)" % (b,a),0)
   finally:
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
   try:
      lock()   
      r = _cmd.fit("(%s in %s)" % (a,b),
                  "(%s in %s)" % (b,a),1)
   finally:
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
   try:
      lock()   
      r = _cmd.fit_pairs(arg)
   finally:
      unlock()
   return r

def expfit(a,b):
   '''
   ??? OBSOLETE
   '''
   try:
      lock()   
      r = _cmd.fit(a,b,2)
   finally:
      unlock()
   return r
   
def zoom(*arg):
   '''
DESCRIPTION
  
   "zoom" scales and translates the window and the origin to cover the
   atom selection.
      
USAGE
 
   zoom object-or-selection
   zoom (selection)
 
PYMOL API

   cmd.orient( string object-or-selection )
   '''
   if len(arg):
      a=arg[0]
   else:
      a="all"
   try:
      lock()   
      r = _cmd.zoom(a)
   finally:
      unlock()
   return r
   
def frame(a):
   '''
DESCRIPTION
  
   "frame" sets the viewer to the indicated movie frame.
   
USAGE
 
   frame frame-number
 
PYMOL API
 
   cmd.frame( int frame_number )
 
NOTES
 
   Frame numbers are 1-based 
   '''
   try:
      lock()   
      r = _cmd.frame(int(a))
   finally:
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
 
   cmd.move( string axis, float distance )
   '''
   try:
      lock()   
      r = _cmd.move(a,float(b))
   finally:
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
 
   cmd.clip( string plane, float distance )
   '''
   try:
      lock()   
      r = _cmd.clip(a,float(b))
   finally:
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
 
   cmd.origin( string object-or-selection )
   '''
   try:
      lock()   
      r = _cmd.origin(a)
   finally:
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
 
   cmd.orient( string object-or-selection )
   '''
   try:
      lock()
      if len(arg)<1:
         a = "(all)"
      else:
         a = arg[0]
      r = _cmd.orient(a)
   finally:
      unlock()
   return r

def is_glut_thread():
   '''
INTERNAL
   '''
   if thread.get_ident() == pymol.glutThread:
      return 1
   else:
      return 0

def refresh():
   '''
INTERNAL
   '''
   try:
      lock()
      if thread.get_ident() == pymol.glutThread:
         r = _cmd.refresh_now()
      else:
         r = _cmd.refresh()
   finally:
      unlock()
   return r

def dirty():
   '''
INTERNAL
   '''
   try:
      lock()   
      r = _cmd.dirty()
   finally:
      unlock()
   return r

def set(a,b):
   '''
DESCRIPTION
  
   "set" changes one of the PyMOL state variables
      
USAGE
 
   set variable = value
 
PYMOL API
 
   cmd.set ( string variable, string value )
   '''
   try:
      lock()   
      r = _cmd.set(a,b)
   finally:
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
 
   cmd.reset ( )
   '''
   try:
      lock()   
      r = _cmd.reset(0)
   finally:
      unlock()
   return r

def meter_reset():
   '''
UNDOCUMENTED
   '''
   try:
      lock()   
      r = _cmd.reset_rate()
   finally:
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
 
   cmd.delete ( string object-or-selection-name )
   '''
   try:
      lock()   
      r = _cmd.delete(a)
   finally:
      unlock()
   return r

def _quit():
   try:
      lock()
      r = _cmd.quit()
   finally:
      unlock()
   return r

def quit():
   '''
DESCRIPTION
  
   "quit" terminates the program. 
   
USAGE
 
   quit
 
PYMOL API
 
   cmd.quit()
   '''
   if thread.get_ident() ==pymol.glutThread:
      lock()
      r = _cmd.do("_quit")
   else:
      r = _cmd.do("_quit")
      thread.exit()
   return r

def png(a):
   '''
DESCRIPTION
  
   "png" writes a png format image file of the current image to disk.
   
USAGE
 
   png filename
 
PYMOL API
 
   cmd.png( string filename )
   '''
   if thread.get_ident() ==pymol.glutThread:
      r = _png(a)
   else:
      r = _cmd.do("cmd._png('"+a+"')")
   return r

def _png(a):
   try:
      lock()   
      fname = a
      if not re.search("\.png$",fname):
         fname = fname +".png"
      fname = os.path.expanduser(fname)
      fname = os.path.expandvars(fname)         
      r = _cmd.png(fname)
   finally:
      unlock()
   return r

def mclear():
   '''
DESCRIPTION
  
   "mclear" clears the movie frame image cache.
   
USAGE
 
   mclear
 
PYMOL API
 
   cmd.mclear()
   '''
   try:
      lock()   
      r = _cmd.mclear()
   finally:
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
 
   cmd.set_key( string key, function fn, tuple arguments)
 
PYTHON EXAMPLE
 
   from pymol import cmd
 
   def color_blue(object):
      cmd.color("blue",object)
    
   cmd.set_key( 'F1' , make_it_blue, ( "object1" ) )
   cmd.set_key( 'F2' , make_it_blue, ( "object2" ) )
 
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
 
   cmd.mstop()
   '''
   try:
      lock()   
      r = _cmd.mplay(0)
   finally:
      unlock()
   return r

def mplay():
   '''
DESCRIPTION
  
   "mplay" starts the movie.
   
USAGE
 
   mplay
 
PYMOL API
 
   cmd.mplay()
   '''
   try:
      lock()   
      r = _cmd.mplay(1)
   finally:
      unlock()
   return r

def mray():
   '''
DEPRECATED
   '''
   try:
      lock()   
      r = _cmd.mplay(2)
   finally:
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
  
   cmd.viewport(int width, int height)
   '''
   r = _cmd.viewport(int(a),int(b))
   
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
  
   cmd.mdo( int frame, string command )
 
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
   try:
      lock()   
      r = _cmd.mdo(int(a)-1,b)
   finally:
      unlock()
   return r

def dummy(*arg):
   '''
DEBUGGING
   '''
   try:
      lock()   
   finally:
      unlock()
   return None

def rock():
   '''
DESCRIPTION
  
   "rock" toggles Y axis rocking.
 
USAGE
 
   rock
 
PYMOL API
  
   cmd.rock()
   '''
   try:
      lock()   
      r = _cmd.rock()
   finally:
      unlock()
   return r

def forward():
   '''
DESCRIPTION
  
   "forward" moves the movie one frame forward.
 
USAGE
 
   forward
 
PYMOL API
  
   cmd.forward()
   '''
   try:
      lock()   
      r = _cmd.setframe(5,1)
   finally:
      unlock()
   return r

def backward():
   '''
DESCRIPTION
  
   "backward" moves the movie back one frame.
 
USAGE
 
   backward
 
PYMOL API
  
   cmd.backward()
   '''
   try:
      lock()   
      r = _cmd.setframe(5,-1)
   finally:
      unlock()
   return r


def rewind():
   '''
DESCRIPTION
  
   "rewind" goes to the beginning of the movie.
 
USAGE
 
   rewind
 
PYMOL API
  
   cmd.rewind()
   '''
   try:
      lock()   
      r = _cmd.setframe(0,0)
   finally:
      unlock()
   return r

def ending():
   '''
DESCRIPTION
  
   "ending" goes to the end of the movie.
 
USAGE
 
   ending
 
PYMOL API
  
   cmd.ending()
   '''
   try:
      lock()   
      r=_cmd.setframe(2,0)
   finally:
      unlock()
   return r

def middle():
   '''
DESCRIPTION
  
   "middle" goes to the middle of the movie.
 
USAGE
 
   middle
 
PYMOL API
  
   cmd.middle()
   '''
   try:
      lock()   
      r = _cmd.setframe(3,0)
   finally:
      unlock()
   return r

def test(): # generic test routine for development
   '''
DEBUGGING
   '''
   try:
      lock()   
      r=_cmd.test()
   finally:
      unlock()
   return r

def dump(fnam,obj):
   '''
DEBUGGING
   '''
   try:
      lock()
      r = _cmd.dump(fnam,obj)
   finally:
      unlock()
   return r

def save(*arg):
   '''
DESCRIPTION
  
   "save" writes selected atoms to a PDB file
 
USAGE
 
   save filename [,(selection) [,state] ]
 
PYMOL API
  
   cmd.save(filename, selection, state)
   '''
   r = 1
   try:
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
         fname = os.path.expanduser(fname)
         fname = os.path.expandvars(fname)
         f=open(fname,"w")
         if f:
            f.write(_cmd.get_pdb(sele,int(state)-1))
            f.close()
            r = None
            print " Save: wrote \""+fname+"\"."
   finally:
      unlock()
   return r


def get_model(*arg):
   '''
DESCRIPTION
  
   "get_model" returns a Chempy "Indexed" format model from a selection.
 
PYMOL API
 
   cmd.get_model( selection [,state] )
 
   '''
   r = 1
   try:
      lock()
      sele = "(all)"
      state = -1
      if len(arg)==1:
         sele = arg[0]
      elif len(arg)==2:
         sele = arg[0]
         state = arg[1]
      r = _cmd.get_model(sele,int(state)-1)
   finally:
      unlock()
   return r

def create(*arg):
   '''
DESCRIPTION
  
   "create" creates a new object from a selection
 
   NOTE: this command has not yet been throughly tested
 
USAGE
 
   create name = (selection) [,source_state [,target_state] ]
 
PYMOL API
  
   cmd.create(string name,string selection,int state,int target_state)
   '''
   name = arg[0]
   argst = string.join(arg[1:],',')
   arg = _split(argst,',')
   source = -1
   target = -1
   la = len(arg)
   sele = arg[0]
   if la>1:
      source = int(arg[1])-1
   if la>2:
      target = int(arg[2])-1
   try:
      lock()
      _cmd.create(name,sele,source,target)
   finally:
      unlock()
   return None

def get_feedback():
   l = []
   try:
      lock()
      r = _cmd.get_feedback()
      while r:
         l.append(r)
         r = _cmd.get_feedback()
   finally:
      unlock()
   return l

def load_model(*arg):
   '''
   '''
   r = 1
   try:
      lock()   
      ftype = 8
      state = -1
      model = arg[0];
      if len(arg)==2:
         oname = string.strip(arg[1])
         r = _cmd.load_object(oname,model,state,ftype)
      elif len(arg)==3:
         oname = string.strip(arg[1])
         state = int(arg[2])-1
         r = _cmd.load(oname,model,state,ftype)
      else:
         print "Error: invalid arguments."
   finally:
      unlock()
   return r
   
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
  
   cmd.load( filename [,object [,state]] )
   '''
   r = 1
   try:
      lock()   
      ftype = 0
      state = -1
      fname = arg[0];
      fname = os.path.expanduser(fname)
      fname = os.path.expandvars(fname)
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
         r = _cmd.load(oname,fname,state,ftype)
      elif len(arg)==2:
         oname = string.strip(arg[1])
         r = _cmd.load(oname,fname,state,ftype)
      elif len(arg)==3:
         oname = string.strip(arg[1])
         state = int(arg[2])-1
         r = _cmd.load(oname,fname,state,ftype)
      elif len(arg)==4:
         if loadable.has_key(arg[3]):
            ftype = loadable[arg[3]]
         else:
            ftype = int(arg[3])
         state = int(arg[2])-1
         oname = string.strip(arg[1])
         r = _cmd.load(oname,fname,state,ftype)
      else:
         print "argument error."
   finally:
      unlock()
   return r

def read_molstr(*arg):
   r = 1
   try:
      lock()
      ftype = 3
      if len(arg)==2:
         oname = string.strip(arg[1])
         r = _cmd.load(oname,arg[0],-1,ftype)
      elif len(arg)==3:
         oname = string.strip(arg[1])
         r = _cmd.load(oname,arg[0],int(arg[2])-1,ftype)
      else:
         print "argument error."
   finally:
      unlock()
   return r

def read_mmodstr(*arg):
   r = 1
   try:
      lock()   
      ftype = 6
      if len(arg)==2:
         oname = string.strip(arg[1])
         r = _cmd.load(oname,arg[0],-1,ftype)
      elif len(arg)==3:
         oname = string.strip(arg[1])
         r = _cmd.load(oname,arg[0],int(arg[2])-1,ftype)
      else:
         print "argument error."
   finally:
      unlock()
   return r
   
def select(*arg):
   '''
DESCRIPTION
  
   "select" creates a named selection from an atom selection.
 
USAGE
 
   select (selection)
   select selection-name = (selection)
 
PYMOL API
  
   cmd.select(string selection-name, string selection)
 
EXAMPLES 
 
   select near = (pk1 expand 8)
   select bb = (name ca,n,c,o )
   '''
   try:
      lock()   
      if len(arg)==1:
         sel_cnt = _cmd.get("sel_counter") + 1.0
         _cmd.set("sel_counter","%1.0f" % sel_cnt)
         sel_name = "sel%02.0f" % sel_cnt
         sel = arg[0]
      else:
         sel_name = arg[0]
         sel = arg[1]
      r = _cmd.select(sel_name,sel)
   finally:
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
  
   cmd.color( string color, string color-name )
 
EXAMPLES 
 
   color yellow, (name C*)
   '''
   try:
      lock()   
      if len(arg)==2:
         r = _cmd.color(arg[0],arg[1],0)
      else:
         r = _cmd.color(arg[0],"(all)",0)   
   finally:
      unlock()
   return r

def flag(*arg):
   '''
DESCRIPTION
  
   "flag" sets the indicated flag for atoms in the selection and
    clears the indicated flag for atoms not in the selection.  This
    is primarily useful for passing selection information into
    Chempy models.
   
USAGE

   flag flag_number = selection 
    
PYMOL API
  
   cmd.flag( int flag, string selection )
 
EXAMPLES  
 
   flag 0 = (name ca)
   flag 1 = (resi 45 x; 6)
 
   '''
   try:
      lock()   
      if len(arg)==2:
         r = _cmd.flag(int(arg[0]),arg[1])
   finally:
      unlock()
   return r

def set_color(nam,col):
   '''
DESCRIPTION
  
   "set_color" defines a new color with color indices (0.0-1.0)
   
USAGE
 
   set_color color-name = [ red-float, green-float, blue-float ]
 
PYMOL API
  
   cmd.set_color( string color, float-list color-components )
 
EXAMPLES 
 
   set_color red = [ 1.0, 0.0, 0.0 ]
   '''
   r = 1
   if isinstance(col,types.StringType):
      col = eval(col)
   if not (isinstance(col,types.ListType) or isinstance(col,types.TupleType)):
      print "Error: color specification must be a list such as [ 1.0, 0.0, 0.0 ]"
   elif len(col)!=3:
      print "Error: color specification must be a list such as [ 1.0, 0.0, 0.0 ]"
   else:
      try:
         lock()

         if len(col)==3:
            r = _cmd.colordef(nam,float(col[0]),float(col[1]),float(col[2]))
         else:
            print "Error: invalid color."
      finally:
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
 
   cmd.mpng( string prefix )
   '''
   if thread.get_ident() ==pymol.glutThread:
      r = _mpng(a)
   else:
      r = _cmd.do("cmd._mpng('"+a+"')")
   return r

def _mpng(*arg):
   try:
      lock()   
      fname = arg[0]
      if re.search("\.png$",fname):
         fname = re.sub("\.png$","",fname)
      fname = os.path.expanduser(fname)
      fname = os.path.expandvars(fname)
      r = _cmd.mpng_(fname)
   finally:
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
 
   cmd.show( string representation, string object-or-selection )
 
EXAMPLES
 
   show lines,(name ca or name c or name n)
   show ribbon
 
NOTES
 
   "show" alone will turn on lines for all bonds.
   '''
   r=1
   try:
      lock()
      l = len(arg)
      if not l:
         r = _cmd.showhide("(all)",0,1); # show lines by default       
      elif l==2:
         rep = arg[0]
         if rephash.has_key(rep):
            rep = rephash[rep]
         if repres.has_key(rep):      
            repn = repres[rep];
            r = _cmd.showhide(arg[1],repn,1);
         else:
            print "Error: unrecognized or ambiguous representation"
      elif arg[0]=='all':
         r = _cmd.showhide("(all)",0,1); # show lines by default 
      elif arg[0][0]=='(':
         r = _cmd.showhide(arg[0],0,1);
      else:
         rep = arg[0]
         if rephash.has_key(rep):
            rep = rephash[rep]
         if repres.has_key(rep):      
            repn = repres[rep];
            r = _cmd.showhide("(all)",repn,1);
         else:
            print "Error: unrecognized or ambiguous representation"
   finally:
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
 
   cmd.hide( string representation, string object-or-selection )
 
EXAMPLES
 
   hide lines,all
   hide ribbon
   '''
   r = 1
   l = len(arg)
   try:
      lock()
      if not l:
         r = _cmd.showhide("!",0,0);      
      elif l==2:
         rep = arg[0]
         if rephash.has_key(rep):
            rep = rephash[rep]
         if repres.has_key(rep):      
            repn = repres[rep];
            r = _cmd.showhide(arg[1],repn,0);
         else:
            print "Error: unrecognized or ambiguous representation"
      elif arg[0]=='all':
         r = _cmd.showhide("!",0,0);
      elif arg[0][0]=='(':
         r = _cmd.showhide(arg[0],-1,0);
      else:
         rep = arg[0]
         if rephash.has_key(rep):
            rep = rephash[rep]
         if repres.has_key(rep):
            repn = repres[rep];
            r = _cmd.showhide("(all)",repn,0);
         else:
            print "Error: unrecognized or ambiguous representation"
   finally:
      unlock()
   return r

def mmatrix(a):
   '''
DESCRIPTION
  
   "mmatrix" sets up a matrix to be used for the first frame of the movie.
   
USAGE
 
   mmatrix {clear|store|recall}
 
PYMOL API
 
   cmd.mmatrix( string action )
 
EXAMPLES
 
   mmatrix store
   '''
   r = 1
   try:
      lock()   
      if a=="clear":
         r = _cmd.mmatrix(0)
      elif a=="store":
         r = _cmd.mmatrix(1)
      elif a=="recall":
         r = _cmd.mmatrix(2)
   finally:
      unlock()
   return r

def enable(*arg):
   '''
DESCRIPTION
  
   "enable" enable display of an object and all currently visible representations.
   
USAGE
 
   enable object
   enable all
   
PYMOL API
 
   cmd.enable( string object-name )
 
EXAMPLE
 
   enable my_object
   '''
   if len(arg):
      nam = arg[0]
   else:
      nam = 'all'
   try:
      lock()   
      r = _cmd.onoff(nam,1);
   finally:
      unlock()
   return r

def disable(*arg):
   '''
DESCRIPTION
  
   "disable" disables display of an object and all currently visible representations.
   
USAGE
 
   disable object
   disable all 
 
PYMOL API
 
   cmd.disable( string object-name )
 
EXAMPLE
 
   disable my_object
   '''
   if len(arg):
      nam = arg[0]
   else:
      nam = 'all'
   try:
      lock()   
      r = _cmd.onoff(nam,0);
   finally:
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
 
   cmd.mset( string specification )
 
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
   try:
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
      r = _cmd.mset(string.join(output," "))
   finally:
      unlock()
   return r

def null():
   pass

keyword = { 
   'alter'         : [alter        , 2 , 2 , ',' , 0 ],
   'alter_state'   : [alter_state  , 3 , 3 , ',' , 0 ],
   'api'           : [api          , 0 , 0 , ',' , 0 ],
   'backward'      : [backward     , 0 , 0 , ',' , 0 ],
   'clip'          : [clip         , 2 , 2 , ',' , 0 ],
   'cls'           : [cls          , 0 , 0 , ',' , 0 ],
   'color'         : [color        , 1 , 2 , ',' , 0 ],
   'commands'      : [commands     , 0 , 0 , ',' , 0 ],
   'copy'          : [copy         , 2 , 2 , '=' , 0 ],
   'count_states'  : [count_states , 0 , 1 , ',' , 0 ],
   'create'        : [create       , 2 , 2 , '=' , 0 ],   
   'delete'        : [delete       , 1 , 1 , ',' , 0 ],
   'disable'       : [disable      , 0 , 1 , ',' , 0 ],
   'dist'          : [dist         , 0 , 2 , '=' , 0 ],
   'distance'      : [distance     , 0 , 2 , '=' , 0 ],
   'dump'          : [dump         , 2 , 2 , ',' , 0 ],
   'enable'        : [enable       , 0 , 1 , ',' , 0 ],
   'ending'        : [ending       , 0 , 0 , ',' , 0 ],
   'export_dots'   : [export_dots  , 2 , 2 , ',' , 0 ],
   'fit'           : [fit          , 2 , 2 , ',' , 0 ],
   'flag'          : [flag         , 2 , 2 , '=' , 0 ],
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
   'set_color'     : [set_color    , 2 , 2 , '=' , 0 ],
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
   'zoom'          : [zoom         , 0 , 1 , ',' , 0 ]
   }

help_only = {
   'selections'    : [selections   , 0 , 0 , ',' , 0 ],
   'keyboard'      : [keyboard     , 0 , 0 , ',' , 0 ],
   'mouse'         : [mouse        , 0 , 0 , ',' , 0 ],
   'examples'      : [examples     , 0 , 0 , ',' , 0 ],
   'launching'     : [launching    , 0 , 0 , ',' , 0 ],
   'movies'        : [movies       , 0 , 0 , ',' , 0 ],
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
   'xplor' : 7,
   'model' : 8,
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

