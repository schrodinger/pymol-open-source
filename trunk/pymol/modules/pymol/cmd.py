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

#
# Return conventions:
#   In general, functions should return "None" for error conditions
#   or something else if they succeed.  As of 2/07/2001, return codes
#   are a bit of a mess...
#

import re
import _cmd
import string
import traceback
import thread
import threading
import types
import pymol
import os
from glob import glob

from chempy import io
from chempy.sdf import SDF,SDFRec
from chempy import fragments

class QuietException:
   def __init__(self,args=None):
      self.args = args

lock_api = pymol.lock_api

def write_html_ref(file):
   lst = globals()
   f=open(file,'w')
   head = 'H2'
   f.write("<HTML><BODY><H1>Reference</H1>")
   kees = lst.keys()
   kees.sort()
   for a in kees:
      if hasattr(lst[a],'__doc__'):
         if a[0:1]!='_' and (a not in ['string','thread',
                                       'setup_global_locks',
                                       'real_system',
                                       ]):
            doc = lst[a].__doc__
            if type(doc) is types.StringType:
               if len(doc):
                  doc = string.strip(doc)
                  doc = string.replace(doc,"<","&lt;")
                  f.write("<HR SIZE=1><%s>"%head+a+"</%s>\n"%head)
                  f.write("<PRE>"+string.strip(doc)+"\n\n</PRE>")
   f.write("</BODY></HTML>")
   f.close()
            
def commands():
   '''
COMMANDS
 
   INPUT/OUTPUT  load      save      delete    quit
   VIEW          turn      move      clip      rock
                 show      hide      enable    disable
                 reset     refresh   rebuild
                 zoom      origin    orient
   MOVIES        mplay     mstop     mset      mdo
                 mpng      mmatrix   frame
                 rewind    middle    ending
                 forward   backward
   IMAGING       png       mpng
   RAY TRACING   ray      
   MAPS          isomesh   isodot
   DISPLAY       cls       viewport  splash    
   SELECTIONS    select    mask   
   SETTINGS      set       button
   ATOMS         alter     alter_state 
   EDITING       create    replace   remove    h_fill   remove_picked
                 edit      bond      unbond    h_add    fuse       
                 undo      redo      protect   cycle_valence  
   FITTING       fit       rms       rms_cur   pair_fit  
                 intra_fit intra_rms intra_rms_cur   
   COLORS        color     set_color
   HELP          help      commands
   DISTANCES     dist      
   STEREO        stereo
   SYMMETRY      symexp
   SCRIPTS       @         run

Try "help <command-name>".  Also see the following extra topics:
 
   "movies", "keyboard", "mouse", "selections",
   "examples", "launching", "editing", and "api".
'''
   help('commands')

def editing():
   '''
SUMMARY

PyMOL has a minimal but functional molecular structure editing
capability.  However, you will need to set up Tinker if you want to to
"clean-up" your structures after editing.  Furthermore, if you are
going to modify molecules other than proteins, then you will also need
a way of assigning atom types (Amber) on the fly.  Unfortunately, my
solution to that problem hasn't been published yet.

To edit a conformation or structure, you first need to enter editing
mode (see Mouse Menu).  Then you can pick an atom (CTRL-Middle click)
or a bond (CTRL-Right click).  Next, you can use the other
CTRL-key/click combinations listed on the right hand side of the
screen to adjust the attached fragments.  For example, CTRL-left click
will torsion fragments.

Editing structures is done through a series of CTRL key actions
applied to the currently selected atom or bonds. See "help edit_keys"
for the exact combinations.  To build structure, you usually just
replace hydrogens with methyl groups, etc., and then repeat.

NOTE
  
Only "lines" representations can be picked.
   
'''

def release():
   '''
RELEASE NOTES

PyMOL is a free, open, and expandable molecular graphics system
written by a computational scientist to enable molecular modeling and
visualization from directly within Python.  It will be of most benefit
to hybrid scientist/developers in the fields of structural biology,
molecular modeling, computational chemistry, and informatics who want
a completely unrestricted visualization tool capable of working
directly with their own programs via Python.  It will also be of benefit
to advanced non-developers familiar with similar programs such as
Midas, O, Grasp, X-PLOR and CNS.

Due to PyMOL's current "user-unfriendliness", this release is ONLY
appropriate for those who prefer to use text commands and scripts, and
for developers who want to integrate PyMOL's visualization and
molecular editing capabilities with their own work.

PyMOL currently includes a diverse command language, a powerful
application programmers interface (API), and a variety of mouse and
keyboard driven functionality for viewing, animation, rendering, and
molecular editing.  However, this release of PyMOL does NOT include an
adequate graphical user interface, menu bar, test suite, or a complete 
help system.  Such enhancements are in progress, but proceed
at a slow pace.  A manual is now available on the web site.

Two external GUI development options are supported for PyMOL:
"Tkinter" and "wxPython".  Developers can take their pick.  I am
committed to insuring that PyMOL will work with both of them, but it
is unlikely that I will have time to develop a complete external GUI
myself any time soon using either toolkit. 

Note that only Tkinter is supported under Windows with the default
PyMOL and Python distributions, so for maximum ease of installation
under Windows, stick with Tkinter (Tcl/Tk).

Warren L. DeLano (2/21/2001), warren@delanoscientific.com
'''

def edit_keys():
   '''
EDITING KEYS 

   These are defaults, which can be redefined.  Note that while
entering text on the command line, some of these control keys take on
text editing functions instead (CTRL - A, E, and K)

ATOM REPLACEMENT
 
   CTRL-C    Replace picked atom with carbon   (C)
   CTRL-N    Replace picked atom with nitrogen (N)
   CTRL-O    Replace picked atom with oxygen   (O)
   CTRL-S    Replace picked atom with sulpher  (S)
   CTRL-G    Replace picked atom with hydrogen (H)
   CTRL-F    Replace picked atom with fluorene (F)
   CTRL-L    Replace picked atom with chlorine (Cl)
   CTRL-B    Replace picked atom with bromine  (Br)
   CTRL-I    Replace picked atom with iodine   (I)
   
ATOM MODIFICATION
  
   CTRL-J    Set charge on picked atom to -1
   CTRL-K    Set charge on picked atom to +1
   CTRL-D    Remove atom or bond (DELETE works too).
   CTRL-Y    Add a hydrogen to the current atom
   CTRL-F    Adjust hydrogens to match valence.

BONDS

   CTRL-T    Connect atoms in the (lb) and (rb) selections.
   CTRL-W    Cycle the bond valence on the picked bond.

MISC

   CTRL-Z    undo the previous conformational change.
             (you can not currently undo modifications).
   CTRL-A    redo the previous conformational change.

'''
   help('edit_keys')


def at_sign():
   '''
DESCRIPTION
 
"@" sources a PyMOL command script as if all of the commands in the
file were typed into the PyMOL command line.

USAGE
  
 @ <script-file>

PYMOL API

Not directly available. Instead, use cmd.do("@...").
 
'''
   help(at_sign)

def run():
   '''
DESCRIPTION
 
"run" executes an external Python script in a local name space,
the global namespace, or in its own namespace (as a module).

USAGE
  
   run <python-script> [, (local | global | module) ]

PYMOL API

Not directly available.  Instead, use cmd.do("run ...").

NOTES

The default mode for run is "global".
 
Due to an idiosyncracy in Pickle, you can not pickle objects
directly created at the main level in a script run as "module",
(because the pickled object becomes dependent on that module).
Workaround: delegate construction to an imported module.

'''
   help(run)

def spawn():
   '''
DESCRIPTION
 
"spawn" launches a Python script in a new thread which will run
concurrently with the PyMOL interpreter. It can be run in its own
namespace (like a Python module, default), a local name space, or in
the global namespace.
 
USAGE
  
 run <python-script> [, (local | global | module )]

PYMOL API

 Not directly available.  Instead, use cmd.do("spawn ...").

NOTES

 The default mode for spawn is "module".
 
 Due to an idiosyncracy in Pickle, you can not pickle objects
 directly created at the main level in a script run as "module",
 (because the pickled object becomes dependent on that module).
 Workaround: delegate construction to an imported module.

'''
   help(spawn)

def api():
   '''
DESCRIPTION
 
The PyMOL Python Application Programming Interface (API) should be
accessed exclusively through the "cmd" module.  Nearly all
command-line functions have a corresponding API method.
 
USAGE
 
   from pymol import cmd
   <result> = cmd.<methods>( <args> ) 
    
API-ONLY METHODS
 
   KEY BINDING   set_key
   (more documentation to come)
 
NOTES
 
Although the PyMOL core is not multi-threaded, the API is threads-safe
and can be called asynchronously by external python programs.  PyMOL
handles the necessary locking to insure that internal states do not
get corrupted.
   '''
   help('api')

def keyboard():
   '''
 KEYBOARD COMMANDS and MODIFIERS
 
   TAB          Toggle onscreen text.
   INSERT       Toggle rocking.
 
   LEFT ARROW, RIGHT ARROW  Go backward or forward one frame.
   HOME,       END          Go to the beginning or end of a movie.
 
 ATOM SELECTIONS (Only work on "lines and nonbonded" representations!)
  
   CTRL/left mouse click    Pick atom and store as selection (lb).
   CTRL/middle mouse click  Pick atom and store as selection (mb).
   CTRL/right mouse click   Pick atom and store as selection (rb).
 
   CTRL-SHIFT/left mouse click    Pick atom and add to selection (lb).
   CTRL-SHIFT/middle mouse click  Pick atom and add to selection (mb).
   CTRL-SHIFT/right mouse click   Pick atom and add to selection (rb).

 EDITING 

   type "help edit_keys" for keyboard shortcuts used in editing.
   
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
 
   (normal mode - when not editing)

   CTRL-left mouse click    Pick atom and store as selection (lb).
   CTRL-middle mouse click  Pick atom and store as selection (mb).
   CTRL-right mouse click   Pick atom and store as selection (rb).
 
   CTRL-SHIFT-left mouse click    Pick atom and add to selection (lb).
   CTRL-SHIFT-middle mouse click  Pick atom and add to selection (mb).
   CTRL-SHIFT-right mouse click   Pick atom and add to selection (rb).
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
 
   pymol.com [-ciqstwx] <file.xxx> [-p <file.py> ] ...
 
   -c   Command line mode, no GUI.  For batch opeations.

   -s   Enable stereo mode (not currently autodetected).

   -r <file.py>[,global|local|module] Run a python program in on startup.
   -l <file.py>[,global|local|module] Spawn a python program in new thread.
   
   -i   Disable the internal OpenGL GUI (object list, menus, etc.)
   -x   Disable the external GUI module.
   -t   Use Tcl/Tk based external GUI module (pmg_tk).
   -w   Use wxPython based external GUI module (pmg_wx).
	-q   Quiet launch.  Suppress splash screen.

   <file> can have one of the following extensions, and all 
   files provided will be loaded or run after PyMOL starts.
    
      .pml    PyMOL command script to be run on startup
      .py     Python program to be run on startup
       
      .pdb    Protein Data Bank format file to be loaded on startup
      .mmod   Macromodel format to be loaded on startup
      .mol    MDL MOL file to be loaded on startup
      .xplor  X-PLOR Map file to be loaded on startup
      .pkl    Pickled ChemPy Model (class "chempy.model.Indexed")
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
   
   Python built-in glob module can be useful for loading movies.

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
      flag <number>               f;
      alt <code>                  -
      numeric_type <numeric type> nt; <numeric type>
      text_type <text type>       tt; <text type>
      b <operator> <value>        -
      formal_charge <op> <value>  fc; <operator> <value>
      partial_charge <op> <value> pc; <operator> <value>
      id <original-index>         -
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
      gap <distance>              -
      in <selection>              -
      like <selection>            l;
 
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
TO DOCUMENT
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

def spheroid(*arg):
   '''
TO DOCUMENT, EXPERIMENTAL
'''
   try:
      print "Warning: 'spheroid' is experimental, incomplete, and unstable."
      lock()
      if len(arg)==0:
         r = _cmd.spheroid("")
      else:
         r = _cmd.spheroid(arg[0])
   finally:
      unlock()
   return r

def cls():
   '''
TO DOCUMENT
'''
   r = None
   try:
      lock()
      r = _cmd.cls()
   finally:
      unlock()
   return r

def fragment(name):
   '''
TO DOCUMENT
'''
   try:
      load_model(fragments.get(name),name)
   except:
      print "Error: unable to load fragment %s" % name

def get_dihedral(s1,s2,s3,s4):
   '''
TO DOCUMENT
'''
   r = None
   try:
      lock()
      r = _cmd.get_dihe(str(s1),str(s2),str(s3),str(s4),0)
   finally:
      unlock()
   return r

def set_dihedral(s1,s2,s3,s4,deg):
   '''
NONFUNCTIONAL
'''
   try:
      lock()
      r = _cmd.set_dihe(str(s1),str(s2),str(s3),str(s4),float(deg),0)
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

def edit_mode(*arg):
   '''
TO DOCUMENT
'''
   try:
      lock()
      r = _cmd.get_setting("button_mode")
      r = int(r)
      if not len(arg):
         if r:
            _cmd.set("button_mode","0")
         else:
            _cmd.set("button_mode","1")            
      else:
         if arg[0]=='on':
            _cmd.set("button_mode","1")
         if arg[0]=='off':
            _cmd.set("button_mode","0")
      config_mouse()
   finally:
      unlock()
   pass

def config_mouse(quiet=0):
   '''
TO DOCUMENT
'''
   # NOTE: PyMOL automatically runs this routine upon start-up
   try:
      lock()
      r = _cmd.get_setting("button_mode")
      r = int(r)
      if not r:
         # visualization
         button('l','','rota')
         button('m','','move')
         button('r','','movz')
         button('l','shft','rotz')
         button('m','shft','move')
         button('r','shft','clip')                  
         button('l','ctrl','lb')
         button('m','ctrl','mb')
         button('r','ctrl','rb')                  
         button('l','ctsh','+lb')
         button('m','ctsh','orig')
         button('r','ctsh','+rb')
         if not quiet:
            print " Mouse: configured for visualization."
      else:
         # editing
         button('l','','rota')
         button('m','','move')
         button('r','','movz')
         button('l','shft','rotf')
         button('m','shft','movf')
         button('r','shft','clip')                  
         button('l','ctrl','torf')
         button('m','ctrl','pkat')
         button('r','ctrl','pkbd')                  
         button('l','ctsh','lb')
         button('m','ctsh','orig')
         button('r','ctsh','rb')
         if not quiet:
            print " Mouse: configured for editing."
   finally:
      unlock()


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
   returns the average distance between all atoms/frames
 
NOTES
 
   "dist" alone will show distances between selections (lb) and (rb)
   created by left and right button atom picks (hold down CTRL)
'''
   # the following needs to be seriously reworked...
   #
   la = len(arg)
   argst = None
   nam = None
   if la<3:
      if la>0:
         if arg[0][0]!='(':
            nam = arg[0]
            if la>1:
               argst = string.join(arg[1:],',')
         else:
            argst = string.join(arg,',')
      if not argst:
         argst = "(lb),(rb)"
      if not nam:
         try:
            lock()
            cnt = _cmd.get("dist_counter") + 1.0
            _cmd.set("dist_counter","%1.0f" % cnt)
            nam = "dist%02.0f" % cnt
         finally:
            unlock()
   else:
      nam = arg[0]
      argst = string.join(arg[1:],',')
   arg = _split(argst,',')
   la = len(arg)
   if la<2:
      print "Error: invalid arguments for dist command."
      raise QuietException
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

def bond(*arg):
   '''
DESCRIPTION
 
"bond" creates a new bond between two selections, each of
which should contain one atom.
 
USAGE
 
 
PYMOL API
 
NOTES
 
'''
   la = len(arg)
   if la==1:
      print "Error: invalid arguments for bond command."
      raise QuietException
   else:
      if la>1:
         sel1 = arg[0]
         sel2 = arg[1]
      else:
         sel1 = "(lb)"
         sel2 = "(rb)"
      order = 1
      if(la>2):
         order = int(arg[2])
      try:
         lock()
         r = _cmd.bond(sel1,sel2,order,1)
      finally:
         unlock()
   return r

def invert(*arg):
   '''
TO DOCUMENT
'''
   la = len(arg)
   if la:
      sel1=arg[0]
   else:
      sel1="(lb)"
   if la>1:
      sel2=arg[1]
   else:
      sel2="(rb)"
   try:
      lock()
      r = _cmd.invert(sel1,sel2,0)
   finally:
      unlock()
   return r

def unbond(*arg):
   '''
DESCRIPTION
 
"unbond" removes all bonds between two selections.

USAGE
  
PYMOL API
 
NOTES
 
'''
   la = len(arg)
   if la==1:
      print "Error: invalid arguments for unbond command."
      raise QuietException
   else:
      if la>1:
         sel1 = arg[0]
         sel2 = arg[1]
      else:
         sel1 = "(lb)"
         sel2 = "(rb)"
      try:
         lock()
         r = _cmd.bond(sel1,sel2,0,0)
      finally:
         unlock()
   return r

def show_help(cmd):
   print "PyMOL>help %s\n" % cmd
   help(cmd)
   print "(Hit TAB to hide)"

def undo():
   '''
TO DOCUMENT
'''
   try:
      lock()
      _cmd.undo(-1)
   finally:
      unlock()

def redo():
   '''
TO DOCUMENT
'''
   try:
      lock()
      _cmd.undo(1)
   finally:
      unlock()

def help(*arg):
   '''
USAGE
 
   help <command>
   '''
   set("text","1")
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
         print " \n",string.strip(doc),"\n"
      else:
         print "Error: sorry no help available on that command."
   elif help_only.has_key(cmd):
      doc = help_only[cmd][0].__doc__
      if doc:
         print " \n",string.strip(doc),"\n"
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
      raise QuietException
   elif len(arg)<3:
      nam= arg[0]
      argst = arg[1]
   else:
      nam= arg[0]
      arg = map(str,arg)
      argst = string.join(arg[1:],',')
   arg = _split(argst,',')
   la = len(arg)
   if la<3:
      print "Error: invalid arguments for symexp command."
      raise QuietException
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

def isomesh(*arg):
   '''
DESCRIPTION
 
"isomesh" creates a mesh isosurface object from a map object.
 
USAGE
 
   isomesh name = map-object, level [,(selection) [,buffer] ] 
   '''
   if len(arg)<2:
      print "Error: invalid arguments for isomesh command."
      raise QuietException
   elif len(arg)<3:
      nam= arg[0]
      argst = arg[1]
   else:
      nam= arg[0]
      arg = map(str,arg)
      argst = string.join(arg[1:],',')   
   arg = _split(argst,',')
   la = len(arg)
   if la<1:
      print "Error: invalid arguments for isomesh command."
      raise QuietException
   else:
      maap=arg[0]
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
         r = _cmd.isomesh(nam,0,maap,mopt,optarg1,optarg2,lvl,0)
      finally:
         unlock()
   return r

def isodot(*arg):
   '''
DESCRIPTION
 
"isomesh" creates a dot isosurface object from a map object.
 
USAGE
 
   isodot name = map-object, level [,(selection) [,buffer] ] 
   '''
   if len(arg)<2:
      print "Error: invalid arguments for isodot command."
      raise QuietException
   elif len(arg)<3:
      nam= arg[0]
      argst = arg[1]
   else:
      nam= arg[0]
      arg = map(str,arg)
      argst = string.join(arg[1:],',')   
   arg = _split(argst,',')
   la = len(arg)
   if la<1:
      print "Error: invalid arguments for isodot command."
      raise QuietException
   else:
      maap=arg[0]
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
         r = _cmd.isomesh(nam,0,maap,mopt,optarg1,optarg2,lvl,1)
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

def label(*arg):
   '''
DESCRIPTION
 
"label" labels one or more atoms properties over a selection using
the python evaluator with a separate name space for each atom.  The
symbols defined in the name space are:
 
   name, resn, resi, chain, q, b, segi, type (ATOM,HETATM) 
   formal_charge, partial_charge, numeric_type, text_type
   
All strings in the expression must be explicitly quoted.  This
operation typically takes several seconds per thousand atoms altered.
 
USAGE

   label (selection),expression
   label expression
   
EXAMPLES
  
   label (chain A),chain
   label (n;ca),"%s-%s" % (resn,resi)
   label (resi 200),"%1.3f" % partial_charge
   '''
   if len(arg)<2:
      sele = "(all)"
      expr = arg[0]
   else:
      sele = arg[0]
      expr = arg[1]
   try:
      lock()
      r = _cmd.label(sele,'label='+expr)
   finally:
      unlock()   
   return r

def alter(sele,expr):
   '''
DESCRIPTION
 
"alter" changes one or more atomic properties over a selection using
the python evaluator with a separate name space for each atom.  The
symbols defined in the name space are:
 
   name, resn, resi, chain, alt,
   q, b, segi, and type (ATOM,HETATM),
   partial_charge, formal_charge,
   text_type, numeric_type
   
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
      r = _cmd.alter(sele,expr,0)
   finally:
      unlock()   
   return r


def iterate(sele,expr):
   '''
DESCRIPTION
 
"iterate" iterates over an expression with a separate name space
for each atom.  However, unlike the "alter" command, atomic properties
can not be altered.  Thus, "iterate" is more efficient than "alter".

It can be used to perform operations and aggregations using
atomic selections, and store the results in any global object,
such as the predefined "stored" object.

The local namespace for "iterate" contains the following names

   name, resn, resi, chain, alt,
   q, b, segi, and type (ATOM,HETATM),
   partial_charge, formal_charge,
   text_type, numeric_type
 
All strings in the expression must be explicitly quoted.  This
operation typically takes a second per thousand atoms.
 
USAGE
 
   iterate (selection),expression
 
EXAMPLES

   stored.net_charge = 0
   iterate (all),stored.net_charge = stored.net_charge + partial_charge

   stored.names = []
   iterate (all),stored.names.append(name)
   
   '''
   try:
      lock()
      r = _cmd.alter(sele,expr,1)
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
      r = _cmd.alter_state(int(state)-1,sele,expr,0)
   finally:
      unlock()   
   return r



def iterate_state(state,sele,expr):
   '''
DESCRIPTION
 
"iterate_state" is to "alter_state" as "iterate" is to "alter"
 
USAGE
 
   iterate_state state,(selection),expression
 
EXAMPLES

   stored.sum_x = 0.0
   iterate 1,(all),stored.sum_x = stored.sum_x + x
   '''
   try:
      lock()
      r = _cmd.alter_state(int(state)-1,sele,expr,1)
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
            print "Error: stereo not available"
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
   state = [ 1,1 ]
   adjust = 0.0
   if len(arg)>3:
      adjust = float(arg[3])
   if len(arg)==3:
      state[0]=int(arg[2][0])
      if state[0]<1: state[0]=1;
      state[1]=int(arg[2][1])
      if state[1]<1: state[1]=1
   try:
      lock()
      r = _cmd.overlap(arg[0],arg[1],state[0]-1,state[1]-1,float(adjust))
   finally:
      unlock()
   return r

def distance(*arg):
   '''
OBSOLETE - TO BE REMOVED
   '''
   la = len(arg)
   if la==0:
      a="lb"
      b="rb"
   elif la==1:
      a=arg[0]
      b="lb"
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
   thred = thread.get_ident()
   if (thred == pymol.glutThread):
      _cmd.flush_now()
      lock_api.release()
   else:
      lock_api.release()
      while _cmd.wait_queue(): # wait till our instruction (if any)
         e = threading.Event() # has been executed before continuing
         e.wait(0.05)
         del e

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

def system(a):
   '''
TO DOCUMENT
   '''
   r = _cmd.system(a)
   return r

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
      r = _cmd.intrafit(str(arg[0]),b,0)
   finally:
      unlock()
   return r

def update(a,b):
   '''
DESCRIPTION
  
   "update" transfers coordinates from one selection to another.
USAGE
 
   update (selection),(selection)
 
EXAMPLES
 
   update target,(variant)

NOTES

   Currently, this applies across all pairs of states.  Fine
   control will be added later.
   
'''
   if a[0]!='(': a="(%"+a+")"
   if b[0]!='(': b="(%"+b+")"
   try:
      lock()   
      r = _cmd.update(a,b,-1,-1)
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

def pair_fit(*arg):
   '''
DESCRIPTION
  
   "pair_fit" fits a set of atom pairs between two models.  Each atom
   in each pair must be specified individually, which can be tedious
   to enter manually.  Script files are recommended when using this
   command.
   
USAGE
 
   pair_fit (selection), (selection), [ (selection), (selection) [ ...] ]
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


def remove(sele):
   '''
DESCRIPTION
  
   "remove" eleminates a selection of atoms from models.
      
USAGE
 
   remove (selection)
 
PYMOL API
  
   cmd.remove( string selection )
    
EXAMPLES
 
   remove ( resi 124 )
   
'''
   r = 1
   try:
      lock()   
      r = _cmd.remove(sele)
   finally:
      unlock()
   return r

def remove_picked():
   '''
DESCRIPTION
  
   "remove_picked" removes the atom or bond currently
   picked for editing.
      
USAGE
 
   remove_picked
 
PYMOL API
  
   cmd.remove_picked()

NOTES

   This function is usually connected to the
   DELETE key and "CTRL-D".
    
'''
   r = 1
   try:
      lock()   
      r = _cmd.remove_picked()
   finally:
      unlock()
   return r

def cycle_valence():
   '''
DESCRIPTION
  
   "cycle_valnce" cycles the valence on the currently selected bond.
      
USAGE
 
   cycle_valence
 
PYMOL API
  
   cmd.remove_picked()

NOTES

   This function is usually connected to the
   DELETE key and "CTRL-D".
    
'''
   r = 1
   try:
      lock()   
      r = _cmd.cycle_valence()
   finally:
      unlock()
   return r


def attach(name,geom,valence):
   '''
DESCRIPTION
  
   "attach" adds a single atom onto the picked atom.
      
USAGE
 
   attach name, geometry, valence
 
PYMOL API
  
   cmd.attach( name, geometry, valence )

NOTES

   see code for details

'''
   r = 1
   try:
      lock()   
      r = _cmd.attach(str(name),int(geom),int(valence))
   finally:
      unlock()
   return r

def fuse(*arg):
   '''
DESCRIPTION
  
   "fuse" joins two objects into one by forming a bond. 
   The first object is moved so as to form an approximately
   resonable bond with the first. 
      
USAGE
 
   fuse (selection), (selection)
 
PYMOL API
  
   cmd.fuse( string selection, string selection )

NOTES

   Each selection must include a single atom in each object.
   The atoms can both be hydrogens, in which case they are
   eliminated, or they can both be non-hydrogens, in which
   case a bond is formed between the two atoms.

'''
   la = len(arg)
   if la==1:
      print "Error: invalid arguments for fuse command."
      raise QuietException
   else:
      if la>1:
         sel1 = arg[0]
         sel2 = arg[1]
      else:
         sel1 = "(lb)"
         sel2 = "(rb)"
      try:
         lock()
         r = _cmd.fuse(sel1,sel2)
      finally:
         unlock()
   return r

def edit(*arg):
   '''
DESCRIPTION
  
   "edit" picks an atom or bond for editing.
      
USAGE
 
   edit (selection) [ ,(selection) ]
 
PYMOL API
  
   cmd.edit( string selection  [ ,string selection ] )

NOTES

   If only one selection is provided, an atom is picked.
   If two selections are provided, the bond between them
   is picked (if one exists).

'''
   sel0=''
   sel1=''
   sel2=''
   sel3=''
   la = len(arg)
   if la>=1:
      sel0=arg[0]
   if la>=2:
      sel1=arg[1]
   if la>=3:
      sel2=arg[2]
   if la>=4:
      sel3=arg[3]
   r = 1
   try:
      lock()   
      r = _cmd.edit(sel0,sel1,sel2,sel3)
   finally:
      unlock()
   return r

def torsion(deg):
   '''
DESCRIPTION
  
   "torsion" rotates the torsion on the bond currently
   picked for editing.  The rotated fragment will correspond
   to the first atom specified when picking the bond (or the
   nearest atom, if picked using the mouse).
      
USAGE
 
   torsion angle
 
PYMOL API
  
   cmd.torsion( float angle )

'''
   try:
      lock()   
      r = _cmd.torsion(float(deg))
   finally:
      unlock()
   return r
   

def h_fill():
   '''
DESCRIPTION
  
   "h_fill" removes and replaces hydrogens on the atom
   or bond picked for editing.  
      
USAGE
 
   h_fill
 
PYMOL API
  
   cmd.h_fill()

NOTES
   
   This is useful for fixing hydrogens after changing
   bond valences.
'''
   r = 1
   try:
      lock()   
      r = _cmd.h_fill()
   finally:
      unlock()
   return r

def h_add(*arg):
   '''
DESCRIPTION
  
   "h_add" uses a primitive algorithm to add hydrogens
   onto a molecule.
      
USAGE
 
   h_add (selection)
 
PYMOL API
  
   cmd.h_add( string selection )

'''
   r = 1
   if len(arg):
      sele = arg[0]
   else:
      sele = "(all)"
   try:
      lock()   
      r = _cmd.h_add(sele)
   finally:
      unlock()
   return r
   
def protect(*arg):
   '''
TO DOCUMENT
'''
   
   if len(arg):
      a=arg[0]
   else:
      a="(all)"
   try:
      lock()   
      r = _cmd.protect(a,1)
   finally:
      unlock()
   return r

def unprotect(*arg):
   '''
TO DOCUMENT
'''
   
   if len(arg):
      a=arg[0]
   else:
      a="(all)"
   try:
      lock()   
      r = _cmd.protect(a,0)
   finally:
      unlock()
   return r

def mask(*arg):
   '''
TO DOCUMENT
'''
   if len(arg):
      a=arg[0]
   else:
      a="(all)"
   try:
      lock()   
      r = _cmd.mask(a,1)
   finally:
      unlock()
   return r

def unmask(*arg):
   '''
TO DOCUMENT
'''
   if len(arg):
      a=arg[0]
   else:
      a="(all)"
   try:
      lock()   
      r = _cmd.mask(a,0)
   finally:
      unlock()
   return r

def replace(name,geom,valence):
   '''
TO DOCUMENT
'''
   r = 1
   try:
      lock()
      r = _cmd.replace(str(name),int(geom),int(valence))
   finally:
      unlock()
   return r

def zoom(*arg):
   '''
DESCRIPTION
  
   "zoom" scales and translates the window and the origin to cover the
   atom selection.
      
USAGE
 
   zoom object-or-selection [,buffer]
   zoom (selection) [,buffer]
 
PYMOL API

   cmd.zoom( string object-or-selection [,float buffer] )
   '''
   b=0.0
   if len(arg):
      a=arg[0]
      if len(arg)>1:
         b=float(arg[1])
   else:
      a="all"
   try:
      lock()   
      r = _cmd.zoom(a,b)
   finally:
      unlock()
   return r

def rename(*arg):
   '''
DESCRIPTION
  
   "rename" creates new atom names which are unique within residues.
      
USAGE

   CURRENT
      rename object-name [ ,force ]
      
      force = 0 or 1 (default: 0)
      
   PROPOSED
      rename object-or-selection,force   

PYMOL API

   CURRENT
      cmd.rename( string object-name, int force )

   PROPOSED
      cmd.rename( string object-or-selection, int force )

NOTES

   To regerate only some atom names in a molecule, first clear them
   with an "alter (sele),name=''" commmand, then use "rename"

'''
   force = 0
   if len(arg):
      a=arg[0]
   else:
      a="(all)"
   if len(arg)==2:
      force=int(arg[1])
   try:
      lock()   
      r = _cmd.rename(a,force)
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
   if thread.get_ident() == pymol.glutThread:
      return 1
   else:
      return 0

def refresh():
   '''
DESCRIPTION
  
   "refresh" causes the scene to be refresh as soon as possible.

USAGE

   refresh

PYMOL API
 
   cmd.refresh()
   '''
   if thread.get_ident() == pymol.glutThread:
      r = _cmd.refresh_now()
   else:
      try:
         lock()
         r = _cmd.do("cmd._refresh()")
      finally:
         unlock()
   return r

def _refresh(swap_buffers=1):
   '''
   INTERNAL - can only be safely called by GLUT thread 
   '''
   try:
      lock()
      if thread.get_ident() == pymol.glutThread:
         if swap_buffers:
            r = _cmd.refresh_now()
         else:
            r = _cmd.refresh()
      else:
         print "Error: Ignoring an unsafe call to cmd._refresh"
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
      r = _cmd.set(a,str(b))
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
   if thread.get_ident() == pymol.glutThread:
      try: 
         lock()
         r = _cmd.quit()
      finally:
         unlock()
   else:
      try:
         lock()
         r = _cmd.do("_cmd.quit()")
      finally:
         unlock()
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
   '''
   INTERNAL - can only be safely called by GLUT thread 
   '''
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

def _ctrl(k):
   if ctrl.has_key(k):
      if ctrl[k][0]:
         if ctrl[k][1]:
            apply(ctrl[k][0],ctrl[k][2])
         else:
            apply(ctrl[k][0],())
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
   '''
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
      r = _cmd.setframe(4,0)
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
  
   "save" writes selected atoms to a file.  The file format is
   autodetected if the extesion is ".pdb" or ".pkl"
 
USAGE
 
   save filename [,(selection) [,state] ]
 
PYMOL API
  
   cmd.save(filename, selection, state)
   '''
   r = 1
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
         format = 'mol'
      elif re.search("\.pkl$",fname):
         format = 'pkl'
   if format=='pdb':
      fname = os.path.expanduser(fname)
      fname = os.path.expandvars(fname)
      f=open(fname,"w")
      if f:
         try:
            lock()
            f.write(_cmd.get_pdb(sele,int(state)-1))
         finally:
            unlock()
            f.close()
         r = None
         print " Save: wrote \""+fname+"\"."
   elif format=='pkl':
      io.pkl.toFile(get_model(sele),fname)
      print " Save: wrote \""+fname+"\"."
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


def get_area(*arg):
   '''
   PRE-RELEASE functionality - API will change
   '''
   la=len(arg)
   sele = "(all)"
   state = 1
   load_b = 0
   if la>0:
      sele = arg[0]
   if la>1:
      state = int(arg[1])
   if la>2:
      load_b = int(arg[2])
   try:
      lock()
      r = _cmd.get_area(sele,int(state)-1,load_b)
   finally:
      unlock()
   return r

def get_names(*arg):
   '''
DESCRIPTION
  
   "get_names" returns a list of object and/or selection names.
 
PYMOL API
 
   cmd.get_names( [string: "objects"|"selections"|"all"] )
 
NOTES
 
   The default behavior is to return only object names.
   
   '''
   la = len(arg)
   mode = 1
   if la:
      if arg[0]=='objects':
         mode = 1
      elif arg[0]=='selections':
         mode = 2
      elif arg[0]=='all':
         mode = 0
   try:
      lock()
      r = _cmd.get_names(mode)
   finally:
      unlock()
   return r

def get_type(name):
   '''
DESCRIPTION
  
   "get_type" returns a string describing the named object or
    selection or the string "nonexistent" if the name in unknown.
 
PYMOL API
 
   cmd.get_type(string object-name)
 
NOTES

   Possible return values are
   
   "object:molecule"
   "object:map"
   "object:mesh"
   "object:distance"
   "selection"
   
   '''
   try:
      lock()
      r = _cmd.get_type(name)
   finally:
      unlock()
   return r
   
   
def id_atom(*arg):
   '''
TO DOCUMENT
   '''
   r = -1
   la = len(arg)
   l = apply(identify,arg)
   ll = len(l)
   if not ll:
      if la:
         print "Error: atom %s not found by id_atom." % arg[0]
      else:
         print "Error: atom not found by id_atom."
      raise QuietException
   elif ll>1:
      if la:
         print "Error: multiple atoms %s found by id_atom." % arg[0]
      else:
         print "Error: multiple atoms found by id_atom."
      raise QuietException
   else:
      r = l[0]
   return r

def identify(*arg):
   '''
DESCRIPTION
  
   "identify" returns a list of atom IDs corresponding to the ID code
   of atoms in the selection
 
PYMOL API
 
   list = cmd.identify( [selection] )
 
   '''
   r = []
   try:
      lock()
      sele = "(all)"
      if len(arg)==1:
         sele = arg[0]
      r = _cmd.identify(sele,0) # 0 = default mode
   finally:
      unlock()
   return r

def get_extent(*arg):
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
      r = _cmd.get_min_max(sele,int(state)-1)
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

def load_coords(*arg):
   '''
TO DOCUMENT
   '''
   r = 1
   try:
      lock()
      ok = 1
      ftype = loadable.model
      state = -1
      model = arg[0];
      if len(arg)<2:
         ok=0
      if len(arg)>=2:
         oname = string.strip(arg[1])
      if len(arg)>=3:
         state = int(arg[2])-1
      if ok:
         r = _cmd.load_coords(oname,model,
                              state,ftype)
      else:
         print "Error: invalid arguments."
   finally:
      unlock()
   return r

def finish_object(obj):
   '''
TO DOCUMENT
   '''
   r = 1
   try:
      lock()   
      r = _cmd.finish_object(obj)
   finally:
      unlock()
   return r

def load_object(*arg): # assume first argument is the object type
   '''
TO DOCUMENT
   '''
   r = 1
   try:
      lock()   
      ftype = arg[0]
      state = -1
      finish = 1
      discrete = 0
      object = arg[1];
      la = len(arg)
      if la>2:
         oname = string.strip(arg[2])
      if la>3:
         state = int(arg[3])-1
      if la>4:
         finish = int(arg[4])
      if la>5:
         discrete = int(arg[5])
      if la>1:
         r = _cmd.load_object(oname,object,state,ftype,finish,discrete)
      else:
         print "Error: invalid arguments."
   finally:
      unlock()
   return r
   
def load_brick(*arg):
   '''
TO DOCUMENT
'''
   
   lst = [loadable.brick]
   lst.extend(list(arg))
   return apply(load_object,lst)

def load_map(*arg):
   '''
TO DOCUMENT
'''
   
   lst = [loadable.map]
   lst.extend(list(arg))
   return apply(load_object,lst)

def load_callback(*arg):
   '''
TO DOCUMENT
'''
   
   lst = [loadable.callback]
   lst.extend(list(arg))
   return apply(load_object,lst)

def load_cgo(*arg):
   '''
TO DOCUMENT
'''
   
   lst = [loadable.cgo]
   lst.extend(list(arg))
   return apply(load_object,lst)

def load_model(*arg):
   '''
DESCRIPTION
  
   "load_model" reads a ChemPy model into an object

 
PYMOL API
  
   cmd.load_model(model, object [,state [,finish [,discrete ]]])
   '''
   lst = [loadable.model]
   lst.extend(list(arg))
   return apply(load_object,lst)

def _load(oname,finfo,state,ftype,finish,discrete):
   r = 1
   if ftype not in (loadable.model,loadable.brick):
      if ftype == loadable.r3d:
         import cgo
         obj = cgo.from_r3d(finfo)
         if obj:
            _cmd.load_object(oname,obj,state,loadable.cgo,finish,discrete)
         else:
            print " load: couldn't load raster3d file."
      else:
         r = _cmd.load(oname,finfo,state,ftype,finish,discrete)
   else:
      try:
         x = io.pkl.fromFile(finfo)
         if isinstance(x,types.ListType) or isinstance(x,types.TupleType):
            for a in x:
               r = _cmd.load_object(oname,a,state,ftype,0,discrete)
               if(state>0):
                  state = state + 1
            _cmd.finish_object(oname)
         else:
            r = _cmd.load_object(oname,x,state,ftype,finish,discrete)            
      except:
         print 'Error: can not load file "%s"' % finfo
   return r

def load(*arg):
   '''
DESCRIPTION
  
   "load" reads several file formats.  The file extension is used to
   determine the format.  PDB files must end in ".pdb", MOL files must
   end in ".mol", Macromodel files must end in ".mmod".  and XPLOR
   maps must end in ".xplor".

   Pickled ChemPy models with a ".pkl" can also be directly read.
 
   If an object is specified, then the file is load into that object.
   Otherwise, an object is created with the same name as the file
   prefix.
 
USAGE
 
   load filename [,object [,state [,type [,finish [,discrete ]]]]]
 
PYMOL API
  
   cmd.load( filename [,object [,state [,type [,finish [,discrete ]]]]] 
   '''
   r = 1
   try:
      lock()   
      ftype = 0
      state = -1
      finish = 1
      discrete = 0
      if len(arg)>4:
         finish=int(arg[4])
      if len(arg)>5:
         discrete=int(arg[5])
      fname = arg[0];
      fname = os.path.expanduser(fname)
      fname = os.path.expandvars(fname)
      if re.search("\.pdb$",arg[0]):
         ftype = loadable.pdb
      elif re.search("\.mol$",arg[0]):
         ftype = loadable.mol
      elif re.search("\.mmod$",arg[0]):
         ftype = loadable.mmod
      elif re.search("\.mmd$",arg[0]):
         ftype = loadable.mmod
      elif re.search("\.xplor$",arg[0]):
         ftype = loadable.xplor
      elif re.search("\.pkl$",arg[0]):
         ftype = loadable.model
      elif re.search("\.r3d$",arg[0]):
         ftype = loadable.r3d
      elif re.search("\.sdf$",arg[0]):
         oname = re.sub("[^/]*\/","",arg[0])
         oname = re.sub("\.sdf$","",oname)
         ftype = loadable.molstr
         sdf = SDF(arg[0])
         while 1:
            rec = sdf.read()
            if not rec: break
            r = _load(oname,string.join(rec.get('MOL'),''),state,ftype,0,1)
         del sdf
         _cmd.finish_object(oname)
         do("zoom (%s)"%oname)
         ftype = -1
      if ftype>=0:
         ok = 1
         if len(arg)<1:
            ok = 0
         if len(arg)==1:
            oname = re.sub("[^/]*\/","",arg[0])
            oname = re.sub("\.pdb$|\.mol$|\.mmod$|\.mmd$|\.xplor$|\.pkl$",
                           "",oname)
         else:
            oname = string.strip(arg[1])
         if len(arg)>2:
            state = int(arg[2])-1
         if len(arg)>3:
            if hasattr(loadable,str(arg[3])):
               ftype = getattr(loadable,str(arg[3]))
            else:
               ftype = int(arg[3])
         if ok:
            r = _load(oname,fname,state,ftype,finish,discrete)
         else:
            print "Error: invalid arguments for dist command."
   finally:
      unlock()
   return r

def read_molstr(*arg):
   '''
DESCRIPTION
  
   "read_molstr" reads an MDL MOL format file as a string
   
PYMOL API ONLY
 
   cmd.read_molstr( string MOL-content, string object name 
   [ ,int state [ ,int finish [ ,int discrete ] ] ] )

NOTES

   "state" is a 1-based state index for the object.

   "finish" is a flag (0 or 1) which can be set to zero to improve
   performance when loading large numbers of objects, but you must
   call "finish_object" when you are done.

   "discrete" is a flag (0 or 1) which tells PyMOL that there will be
   no overlapping atoms in the PDB files being loaded.  "discrete"
   objects save memory but can't be edited.
'''
   r = 1
   try:
      lock()
      ftype = 3
      if len(arg)>1:
         oname = string.strip(arg[1])
      if len(arg)==2:
         r = _cmd.load(oname,arg[0],-1,ftype,1,1)
      elif len(arg)==3:
         r = _cmd.load(oname,arg[0],int(arg[2])-1,ftype,1,1)
      elif len(arg)==4:
         r = _cmd.load(oname,arg[0],int(arg[2])-1,ftype,int(arg[3]),discrete)         
      elif len(arg)==5:
         r = _cmd.load(oname,arg[0],int(arg[2])-1,ftype,int(arg[3]),int(arg[4]))
      else:
         print "argument error."
   finally:
      unlock()
   return r

def read_mmodstr(*arg):
   '''
TO DOCUMENT
'''
   r = 1
   try:
      lock()   
      ftype = 6
      if len(arg)==2:
         oname = string.strip(arg[1])
         r = _cmd.load(oname,arg[0],-1,ftype,1,1)
      elif len(arg)==3:
         oname = string.strip(arg[1])
         r = _cmd.load(oname,arg[0],int(arg[2])-1,ftype,1,1)
      else:
         print "argument error."
   finally:
      unlock()
   return r

def read_pdbstr(*arg):
   '''
DESCRIPTION
  
   "read_pdbstr" reads a pdb file as a string
   
PYMOL API ONLY
 
   cmd.read_pdbstr( string pdb-content, string object name 
   [ ,int state [ ,int finish [ ,int discrete ] ] ] )

NOTES

   "state" is a 1-based state index for the object.

   "finish" is a flag (0 or 1) which can be set to zero to improve
   performance when loading large numbers of objects, but you must
   call "finish_object" when you are done.

   "discrete" is a flag (0 or 1) which tells PyMOL that there will be
   no overlapping atoms in the PDB files being loaded.  "discrete"
   objects save memory but can't be edited.
'''
   r = 1
   finish = 1
   discrete = 0
   if len(arg)>3:
      finish=int(arg[3])
   if len(arg)>4:
      discrete=int(arg[4])
   try:
      lock()   
      ftype = loadable.pdbstr
      if len(arg)==2:
         oname = string.strip(arg[1])
         r = _cmd.load(oname,arg[0],-1,ftype,finish,discrete)
      elif len(arg)>=3:
         oname = string.strip(arg[1])
         r = _cmd.load(oname,arg[0],int(arg[2])-1,ftype,finish,discrete)
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
 
   select near = (ll expand 8)
   select bb = (name ca,n,c,o )
   '''
   try:
      quiet=0
      lock()   
      if len(arg)==1:
         sel_cnt = _cmd.get("sel_counter") + 1.0
         _cmd.set("sel_counter","%1.0f" % sel_cnt)
         sel_name = "sel%02.0f" % sel_cnt
         sel = arg[0]
      else:
         sel_name = arg[0]
         sel = arg[1]
      if len(arg)==3:
         quiet=int(arg[2])
      r = _cmd.select(sel_name,sel,quiet)
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

def rebuild():
   '''
DESCRIPTION

   "rebuild" forces PyMOL to recreate all geometric objects in
   case any of them have gone out of sync.

USAGE
   
   rebuild

PYMOL API

   cmd.rebuild()

'''
   r = 1
   try:
      lock()
      r = _cmd.rebuild()
   finally:
      unlock()
      
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
 
   The available representations are:
    
      lines     spheres   mesh      ribbon
      sticks    dots      surface   labels
      nonbonded nb_spheres
   
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
 
   The available representations are:
    
      lines     spheres   mesh      ribbon
      sticks    dots      surface   labels
      nonbonded nb_spheres
   
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

def paste():
   r=1
   lst = []
   if hasattr(pymol,"machine_get_clipboard"):
      lst = pymol.machine_get_clipboard()
   if len(lst):
      while 1:
         if len(lst[-1]):
            if ord(lst[-1][-1])>32: # trim off final CR
               break;
         else:
            break;
         lst[-1]=lst[-1][:-1]
      _cmd.paste(lst)      
   return r

def button(but,mod,act):
   '''
DESCRIPTION
  
   "button" can be used to redefine what the mouse buttons do.
   
USAGE
 
   button <button-code>,<modifier>,<action>
 
PYMOL API
 
   cmd.button( string button-code, string modifier, string action )
 
NOTES

   button-code: L, M, R
   modifers:    None, Shft, Ctrl, CtSh
   actions:     Rota, Move, MovZ, Clip, RotZ, ClpN, ClpF
                lb,   mb,   rb,   +lb,  +mb,  +rb,
                PkAt, PkBd, RotF, TorF, MovF, Orig

   Switching from visualization to editing mode will redefine the
   buttons, so do not use the built-in switch if you want to preserve
   your custom configuration.

'''
   r=1
   try:
      lock()
      but = string.lower(but)
      but = but[0]
      mod = string.lower(mod)
      act = string.lower(act)
      if not button_code.has_key(but):
         print "Error: unrecognized button name '%s'." % but
      elif not but_mod_code.has_key(mod):
         print "Error: unrecognized button modifier '%s'." % mod
      elif not but_act_code.has_key(act):
         print "Error: unrecognized button action '%s'." % act
      else:
         but_code = button_code[but] + 3*but_mod_code[mod]
         act_code = but_act_code[act]
         r = _cmd.button(but_code,act_code)
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

def check(*arg):
   '''
UNSUPPORTED

This function relies on code that is not currently part of PyMOL/ChemPy
   '''
   # NOTE: the realtime module relies on code that is not yet part of PyMOL/ChemPy
   from chempy.tinker import realtime
   la = len(arg)
   if not la:
      arg = get_names("objects")
      arg = arg[0:1]
      la = len(arg)
   if la:
      if la>1:
         preserve=int(arg[1])
      else:
         preserve=0
      realtime.assign("("+arg[0]+")",preserve)
      realtime.setup("("+arg[0]+")")

def fast_minimize(*arg):
   '''
DESCRIPTION

   Runs Tinker's minimize routine in a separate thread without first 
   executing the setup procedure.  

USAGE

   fast_minimize <object-name>,<max-iterations>,<final-gradient>,<interval>

PYMOL API

   Direct API calls not recommended due to threading issues.

NOTE

   Assumes that the molecular state has already been set up and that no
   subsequent changes have been made to atom types.

   '''
   from chempy.tinker import realtime  
   grad  = 0.01
   iter = 500
   interval = 50
   la = len(arg)
   if not la:
      arg = get_names("objects")
      arg = arg[0:1]
      la = len(arg)
   if la:
      sele  = "("+arg[0]+")"
      if la>1:
         iter = int(arg[1])
      if la>2:
         grad = float(arg[2])
      if la>3:
         interval = int(arg[3])
      t = threading.Thread(target=realtime.mini,args=(iter,grad,interval,arg[0]))
      t.setDaemon(1)
      t.start()
   
def minimize(*arg):
   '''
DESCRIPTION

   Minimizes the indicated object in a separate thread after running an automated
   topology and parameter setup routine which generates a custom parameter file
   for Tinker based on the current text atom types ("text_type" field).

USAGE

   minimize <object-name>,<max-iterations>,<final-gradient>,<interval>

PYMOL API

   Direct API calls not recommended due to threading issues.  Instead, call the 
   ChemPy Tinker (chempy.tinker) module directly.

NOTES

   Only all-atom Amber 94/99 atom types are currently supported.

   '''
   # NOTE: the realtime module relies on code that is not yet part of PyMOL/ChemPy
   from chempy.tinker import realtime  
   grad  = 0.01
   iter = 500
   interval = 50
   la = len(arg)
   if not la:
      arg = get_names("objects")
      arg = arg[0:1]
      la = len(arg)
   if la:
      sele  = "("+arg[0]+")"
      if la>1:
         iter = int(arg[1])
      if la>2:
         grad = float(arg[2])
      if la>3:
         interval = int(arg[3])
      if realtime.setup(sele):
         t = threading.Thread(target=realtime.mini,args=(iter,grad,interval,arg[0]))
         t.setDaemon(1)
         t.start()
      else:
         print " minimize: missing parameters, can't continue"

def cd(dir):
   '''
DESCRIPTION

   Changing the current working directory.

USAGE
   
   cd <path>

   '''
   os.chdir(dir)

def pwd():
   '''
DESCRIPTION

   Print current working directory.

USAGE
   
   pwd

   '''
   print os.getcwd()
   

def ls(*arg):
   '''
DESCRIPTION

   List contents of the current working directory.

USAGE
   
   ls [<pattern>]
   dir [<pattern>]

   '''
   pat = "*"
   if len(arg):
      pat = arg[0]
   for a in glob(pat):
      print a   
      
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
   'bond'          : [bond         , 0 , 3 , ',' , 0 ],
   'button'        : [button       , 3 , 3 , ',' , 0 ],
   'cd'            : [cd           , 1 , 1 , ',' , 0 ],  
   'check'         : [check        , 0 , 2 , ',' , 0 ],
   'clip'          : [clip         , 2 , 2 , ',' , 0 ],
   'cls'           : [cls          , 0 , 0 , ',' , 0 ],
   'color'         : [color        , 1 , 2 , ',' , 0 ],
   'commands'      : [commands     , 0 , 0 , ',' , 0 ],
   'copy'          : [copy         , 2 , 2 , '=' , 0 ],
   'count_states'  : [count_states , 0 , 1 , ',' , 0 ],
   'cycle_valence' : [cycle_valence, 0 , 0 , ',' , 0 ],
   'create'        : [create       , 2 , 2 , '=' , 0 ],   
   'delete'        : [delete       , 1 , 1 , ',' , 0 ],
   'dir'           : [ls          , 0 , 1 , ',' , 0 ],  
   'disable'       : [disable      , 0 , 1 , ',' , 0 ],
   'dist'          : [dist         , 0 , 2 , '=' , 0 ],
   'distance'      : [distance     , 0 , 2 , '=' , 0 ],
   'dump'          : [dump         , 2 , 2 , ',' , 0 ],
   'edit'          : [edit         , 1 , 4 , ',' , 0 ],
   'edit_mode'     : [edit_mode    , 0 , 1 , ',' , 0 ],
   'enable'        : [enable       , 0 , 1 , ',' , 0 ],
   'ending'        : [ending       , 0 , 0 , ',' , 0 ],
   'export_dots'   : [export_dots  , 2 , 2 , ',' , 0 ],
   'fast_minimize' : [fast_minimize, 1,  4 , ',' , 0 ],
   'fit'           : [fit          , 2 , 2 , ',' , 0 ],
   'flag'          : [flag         , 2 , 2 , '=' , 0 ],
   'fork'          : [spawn        , 1 , 2 , ',' , 3 ],
   'forward'       : [forward      , 0 , 0 , ',' , 0 ],
   'fragment'      : [fragment     , 1 , 1 , ',' , 0 ],
   'fuse'          : [fuse         , 0 , 2 , ',' , 0 ],
   'frame'         : [frame        , 1 , 1 , ',' , 0 ],
   'h_add'         : [h_add        , 0 , 1 , ',' , 0 ],
   'help'          : [help         , 0 , 1 , ',' , 0 ],
   'hide'          : [hide         , 0 , 2 , ',' , 0 ],
   'intra_fit'     : [intra_fit    , 1 , 2 , ',' , 0 ],
   'intra_rms'     : [intra_rms    , 1 , 2 , ',' , 0 ],
   'intra_rms_cur' : [intra_rms_cur, 1 , 2 , ',' , 0 ],
   'invert'        : [invert       , 0 , 2 , ',' , 0 ],
   'isodot'        : [isodot       , 2 , 2 , '=' , 0 ],   
   'isomesh'       : [isomesh      , 2 , 2 , '=' , 0 ],
   'iterate'       : [iterate      , 2 , 2 , ',' , 0 ],
   'iterate_state' : [iterate_state, 3 , 3 , ',' , 0 ],
   'label'         : [label        , 1 , 2 , ',' , 0 ],
   'load'          : [load         , 1 , 6 , ',' , 0 ],
   'ls'            : [ls           , 0 , 1 , ',' , 0 ],  
   'mask'          : [mask         , 0 , 1 , ',' , 0 ],
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
   'minimize'      : [minimize     , 0 , 4 , ',' , 0 ],
   'mmatrix'       : [mmatrix      , 1 , 1 , ',' , 0 ],
   'origin'        : [origin       , 1 , 1 , ',' , 0 ],
   'orient'        : [orient       , 0 , 1 , ',' , 0 ],
   'overlap'       : [overlap      , 2 , 3 , ',' , 0 ],
   'pair_fit'      : [pair_fit     , 2 ,98 , ',' , 0 ],
   'protect'       : [protect      , 0 , 1 , ',' , 0 ],
   'pwd'           : [pwd          , 0 , 0 , ',' , 0 ],
   'ray'           : [ray          , 0 , 0 , ',' , 0 ],
   'rebuild'       : [rebuild      , 0 , 0 , ',' , 0 ],
   'redo'          : [redo         , 0 , 0 , ',' , 0 ],
   'refresh'       : [refresh      , 0 , 0 , ',' , 0 ],
   'remove'        : [remove       , 1 , 1 , ',' , 0 ],
   'remove_picked' : [remove_picked, 1 , 1 , ',' , 0 ],
   'rename'        : [rename       , 1 , 2 , ',' , 0 ],
   'replace'       : [replace      , 3 , 3 , ',' , 0 ],
   'reset'         : [reset        , 0 , 0 , ',' , 0 ],
   'rewind'        : [rewind       , 0 , 0 , ',' , 0 ],
   'rock'          : [rock         , 0 , 0 , ',' , 0 ],
   'run'           : [run           , 1 , 2 , ',' , 2 ],
   'rms'           : [rms          , 2 , 2 , ',' , 0 ],
   'rms_cur'       : [rms_cur      , 2 , 2 , ',' , 0 ],
   'save'          : [save         , 0 , 4 , ',' , 0 ],
   'select'        : [select       , 1 , 2 , '=' , 0 ],
   'set'           : [set          , 2 , 2 , '=' , 0 ],
   'set_color'     : [set_color    , 2 , 2 , '=' , 0 ],
   'set_key'       : [set_key      , 2 , 1 , ',' , 0 ], # API only
   'show'          : [show         , 0 , 2 , ',' , 0 ],
   'sort'          : [sort         , 0 , 1 , ',' , 0 ],
   'spawn'         : [spawn        , 1 , 2 , ',' , 3 ],
   'spheroid'      : [spheroid     , 0 , 1 , ',' , 0 ],
   'splash'        : [splash       , 0 , 0 , ',' , 0 ],
   '_special'      : [_special     , 3 , 3 , ',' , 0 ],
   'stereo'        : [stereo       , 1 , 1 , ',' , 0 ],
   'symexp'        : [symexp       , 2 , 2 , '=' , 0 ],
   'system'        : [system       , 1 , 1 , ',' , 0 ],
   'test'          : [test         , 0 , 0 , ',' , 0 ],
   'torsion'       : [torsion      , 1 , 1 , ',' , 0 ],
   'turn'          : [turn         , 2 , 2 , ',' , 0 ],
   'quit'          : [quit         , 0 , 0 , ',' , 0 ],
   '_quit'         : [_quit        , 0 , 0 , ',' , 0 ],
   'png'           : [png          , 1 , 1 , ',' , 0 ],
   'unbond'        : [unbond       , 0 , 3 , ',' , 0 ],
   'undo'          : [undo         , 0 , 0 , ',' , 0 ],
   'unmask'        : [unmask       , 0 , 1 , ',' , 0 ],
   'unprotect'     : [unprotect    , 0 , 1 , ',' , 0 ],
   'update'        : [update       , 2 , 2 , ',' , 0 ],
   'viewport'      : [viewport     , 2 , 2 , ',' , 0 ],
   'zoom'          : [zoom         , 0 , 2 , ',' , 0 ],
   }

help_only = {  # for API-only features
   'selections'    : [selections   , 0 , 0 , ',' , 0 ],
   'keyboard'      : [keyboard     , 0 , 0 , ',' , 0 ],
   'mouse'         : [mouse        , 0 , 0 , ',' , 0 ],
   'examples'      : [examples     , 0 , 0 , ',' , 0 ],
   'read_molstr'   : [read_molstr  , 0 , 0 , ',' , 0 ],
   'release'       : [release      , 0 , 0 , ',' , 0 ],   
   'launching'     : [launching    , 0 , 0 , ',' , 0 ],
   'load_model'    : [load_model   , 0 , 0 , ',' , 0 ],
   'movies'        : [movies       , 0 , 0 , ',' , 0 ],
   'editing'       : [editing      , 0 , 0 , ',' , 0 ],  
   'edit_keys'     : [edit_keys    , 0 , 0 , ',' , 0 ],
   'get_names'     : [get_names    , 0 , 0 , ',' , 0 ],
   'get_type'      : [get_type     , 0 , 0 , ',' , 0 ],  
   '@'             : [at_sign      , 0 , 0 , ',' , 0 ],  
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
   'labels'        : 8,
   'nonbonded'     : 9,
   'nb_spheres'    : 10,
}

button_code = {
   'l' : 0,
   'm' : 1,
   'r' : 2,
   }

but_mod_code = {
   ''      : 0,
   'none'  : 0,
   'shft'  : 1,
   'ctrl'  : 2,
   'ctsh'  : 3
   }

but_act_code = {
   'rota' :  0 ,
   'move' :  1 ,
   'movz' :  2 ,
   'clip' :  3 ,
   'rotz' :  4 ,
   'clpn' :  5 ,
   'clpf' :  6 ,
   'lb'   :  7 ,
   'mb'   :  8 ,
   'rb'   :  9 ,
   '+lb'  : 10 ,
   '+mb'  : 11 ,
   '+rb'  : 12 ,
   'pkat' : 13 ,
   'pkbd' : 14 ,
   'rotf' : 15 ,
   'torf' : 16 ,
   'movf' : 17 ,
   'orig' : 18 ,
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
   11       : [ 'F11'       , redo                   , 0 , None ],
   12       : [ 'F12'       , undo                   , 0 , None ],
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

ctrl = {
   'A' : [ redo                   , 0 , None ],
   'B' : [ replace                , 1 , ('Br',1,1) ],
   'C' : [ replace                , 1 , ('C',4,4) ],
   'D' : [ remove_picked          , 0 , None ],   
   'F' : [ replace                , 1 , ('F',1,1) ],   
   'G' : [ replace                , 1 , ('H',1,1) ],
   'I' : [ replace                , 1 , ('I',1,1) ],
   'J' : [ alter                  , 1 , ('ed1','formal_charge=-1.0') ],
   'K' : [ alter                  , 1 , ('ed1','formal_charge =1.0') ],
   'L' : [ replace                , 1 , ('Cl',1,1)],   
   'N' : [ replace                , 1 , ('N',4,3) ],
   'O' : [ replace                , 1 , ('O',4,2) ],   
   'P' : [ replace                , 1 , ('P',4,1) ],
   'Q' : [ h_add                  , 1 , ("ed1",) ],   
   'R' : [ h_fill                 , 0 , None ],   
   'S' : [ replace                , 1 , ('S',4,2) ],
   'T' : [ bond                   , 0 , None ],   
   'U' : [ alter                  , 1 , ('ed1','formal_charge =0.0') ],
   'W' : [ cycle_valence          , 0 , None ],   
   'X' : [ None                   , 0 , None ],
   'Y' : [ attach                 , 1 , ('H',1,1) ],
   'Z' : [ undo                   , 0 , None ],   
   }

class loadable:
   pdb = 0
   mol = 1
   molstr = 3
   mmod = 4
   xplor = 7
   model = 8
   pdbstr = 9    
   brick = 10    # chempy.brick object
   map = 11      # chempy.map object
   callback = 12 # pymol callback obejct
   cgo = 13      # compiled graphic object
   r3d = 14      # r3d, only used within cmd.py
   
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


