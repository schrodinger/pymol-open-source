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

# NOTE: this cmd module has grown absurdly huge, therefore expect this
# file to be broken down sometime after the 0.51 release into more
# manageable pieces.  Note that cmd.<whatever> will still work -- the
# symbols will be mapped into this namespace even after the code modules
# are separated.

# NEW CALL RETURN CONVENTIONS for _cmd.so C-layer
#
# (1) Calls into C (_cmd) should return results/status and print errors
#     and feedback (according to mask) BUT NEVER RAISE EXCEPTIONS.
#
# (2) In the absence of an expected return value, truth applies:
#        Success is 1, true 
#        Failure is 0, false. None, NULL
#
#     NOTE: if you need more informative failure, then code and
#           document an exception to this rule for your functions!
#
# (3) If _cmd produces a specific return result, be sure to include an
#     error result as one of the possibilities outside the range of the
#     expected return value.  For example, a negative distance or count
#
# (4) cmd.py API wrappers can then raise exceptions based on truth
#     and should return truth for success or None for failure
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
import imp
import parsing
import sys
import copy
import selector
import operator
import time

from shortcut import Shortcut

from glob import glob

from chempy import io
from chempy.sdf import SDF,SDFRec
from chempy import fragments

import setting

file_ext_re= re.compile(string.join([
   "\.pdb$|\.ent$|\.mol$|",
   r"\.PDB$|\.ENT$|\.MOL$|",
   r"\.mmod$|\.mmd$|\.dat$|\.out$|",
   r"\.MMOD$|\.MMD$|\.DAT$|\.OUT$|",
   r"\.xplor$|\.pkl$|\.sdf$|", 
   r"\.XPLOR$|\.PKL$|\.SDF$|",                        
   r"\.r3d$|\.xyz$|\.xyz_[0-9]*$|", 
   r"\.R3D$|\.XYZ$|\.XYZ_[0-9]*$|",
   r"\.cc1$|\.cc2$|", # ChemDraw 3D
   r"\.CC1$|\.CC2$|",
   r"\.pmo$|", # Experimental molecular object format
   r"\.PMO$|",
   r"\.ccp4$|\.CCP4$" # CCP4   
   ],''))

safe_oname_re = re.compile(r"\+|\(|\)|\||\&|\!|\,")  # quash reserved characters
   
QuietException = parsing.QuietException

# the following lock is used by both C and Python to insure that no more than
# one active thread enters PyMOL at a given time. 
# 
lock_api = pymol.lock_api
lock_api_c = pymol.lock_api_c

toggle_dict = {'on':1,'off':0,'1':1,'0':0}
toggle_sc = Shortcut(toggle_dict.keys())

def is_string(obj):
   return isinstance(obj,types.StringType)


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
                                       'real_system','sys','imp','glob'
                                       ]):
            doc = lst[a].__doc__
            if is_string(doc):
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
                 view      get_view  set_view
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
                 undo      redo      protect   cycle_valence  attach
   FITTING       fit       rms       rms_cur   pair_fit  
                 intra_fit intra_rms intra_rms_cur   
   COLORS        color     set_color
   HELP          help      commands
   DISTANCES     dist      
   STEREO        stereo
   SYMMETRY      symexp
   SCRIPTS       @         run
   LANGUAGE      alias     extend

Try "help <command-name>".  Also see the following extra topics:
 
   "movies", "keyboard", "mouse", "selections",
   "examples", "launching", "editing", and "api".
'''
   help('commands')

def editing():
   '''
SUMMARY

PyMOL has a rudimentary, but quite functional molecular structure
editing capability.  However, you will need to use an external mimizer
to "clean-up" your structures after editing.  Furthermore, if you are
going to modify molecules other than proteins, then you will also need
a way of assigning atom types on the fly.

To edit a conformation or structure, you first need to enter editing
mode (see Mouse Menu).  Then you can pick an atom (CTRL-Middle click)
or a bond (CTRL-Right click).  Next, you can use the other
CTRL-key/click combinations listed on the right hand side of the
screen to adjust the attached fragments.  For example, CTRL-left click
will move fragments about the selected torsion.

Editing structures is done through a series of CTRL key actions
applied to the currently selected atom or bonds. See "help edit_keys"
for the exact combinations.  To build structures, you usually just
replace hydrogens with methyl groups, etc., and then repeat.  They are
no short-cuts currently available for building common groups, but that
is planned for later versions.

NOTE
  
Only "lines" and "sticks" representations can be picked using the
mouse, however other representations will not interfere with picking
so long as one of these representation is present underneath.
   
'''

def release():
   '''
RELEASE NOTES

PyMOL is a free, open, and expandable molecular graphics system
written by a computational scientist to enable molecular modeling from
directly within Python.  It will be of most benefit to hybrid
scientist/developers in the fields of structural biology,
computational chemistry, and informatics who seek an open and
unrestricted visualization tool for interfacing with their own
programs.  PyMOL will also be of benefit to advanced non-developers
familiar with similar programs such as Midas, O, Grasp, X-PLOR and
CNS.

Due to PyMOL's current "user-unfriendliness", this release is most
appropriate for those who prefer to use text commands and scripts, and
for developers who want to integrate PyMOL's visualization and
molecular editing capabilities with their own work.

PyMOL currently includes a diverse command language, a powerful
application programmers interface (API), and a variety of mouse and
keyboard driven functionality for viewing, animation, rendering, and
molecular editing.  A partial manual is now available on the web.

Two external GUI development options are supported for PyMOL:
"Tkinter" and "wxPython".  Developers can take their pick.  I am
committed to insuring that PyMOL will work with both of them, but it
is unlikely that I will have time to develop a complete external GUI
myself any time soon using either toolkit.

Note that only Tkinter is supported under Windows with the default
PyMOL and Python distributions, so for maximum ease of installation
under Windows, stick with Tkinter (Tcl/Tk).  For this reason, the
Tkinter-based GUI is going to be the default GUI for standard PyMOL
despite its drawbacks.

Warren L. DeLano (5/1/2001), warren@delanoscientific.com
'''

def edit_keys():
   '''
EDITING KEYS 

   These are defaults, which can be redefined.  Note that while
entering text on the command line, some of these control keys take on
text editing functions instead (CTRL - A, E, and K, and DELETE), so
you should clear the command line before trying to edit atoms.

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
   CTRL-R    Adjust hydrogens on atom/bond to match valence.
   CTRL-E    Inverts the picked stereo center, but you must first
             indicate the constant portions with the (lb) and (rb)
             selections.

   CTRL-T    Connect atoms in the (lb) and (rb) selections.
   CTRL-W    Cycle the bond valence on the picked bond.

UNDO and REDO of conformational changes (not atom changes!)

   CTRL-Z    undo the previous conformational change.
             (you can not currently undo atom modifications).
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
 
   "run" executes an external Python script in a local name space, the
   global namespace, or in its own namespace (as a module).

USAGE
  
   run python-script [, (local | global | module) ]

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
   namespace (like a Python module, default), a local name space, or
   in the global namespace.
 
USAGE
  
   run python-script [, (local | global | module )]

PYMOL API

   Not directly available.  Instead, use cmd.do("spawn ...").

NOTES

   The default mode for spawn is "module".
 
   Due to an idiosyncracy in Pickle, you can not pickle objects
   directly created at the main level in a script run as "module",
   (because the pickled object becomes dependent on that module).
   Workaround: delegate construction to an imported module.

   The best way to spawn processes at startup is to use the -l option
   (see "help launching").
'''
   help(spawn)

def api():
   '''
DESCRIPTION
 
   The PyMOL Python Application Programming Interface (API) should be
   accessed exclusively through the "cmd" module (never "_cmd"!).  Nearly
   all command-line functions have a corresponding API method.
 
USAGE
 
   from pymol import cmd
   result = cmd.<command-name>( argument , ... ) 
    
NOTES
 
   Although the PyMOL core is not multi-threaded, the API is
   thread-safe and can be called asynchronously by external python
   programs.  PyMOL handles the necessary locking to insure that
   internal states do not get corrupted.  This makes it very easy to
   build complicated systems which involve direct realtime visualization.

   '''
   
   help('api')

def keyboard():
   '''
KEYBOARD COMMANDS and MODIFIERS
 
   ESC          Toggle onscreen text.
   INSERT       Toggle rocking.

   LEFT ARROW, RIGHT ARROW    Go backward or forward one frame, or when
                              editing, go forward or back one character.
   HOME, END    Go to the beginning or end of a movie.

   Command Entry Field in the Interal GUI (black window)
   
   TAB          Complete commmand or filename (like in tcsh or bash).
   CTRL-A       Go to the beginning of the line.
   CTRL-E       Go to the end of the line.
   CTRL-K       Delete through to the end of the line.

   Command Entry Field on the External GUI (gray window).
   
   CTRL-C       These operating system-provided cut and paste functions
   CTRL-V       will only work in the external GUI command line.
   
EDITING 

   type "help edit_keys" for keyboard shortcuts used in editing.
   
   '''
   help('keyboard')


def mouse():
   '''
MOUSE CONTROLS
 
  The configuration can be changed using the "Mouse" menu.  The
  current configuration is described on screen with a small matrix on
  the lower right hand corner, using the following abbreviations:

   Buttons (Horizontal Axis)
   
      L        = left mouse click
      M        = middle mouse click
      R        = right mouse click

   Modifiers (Veritical axis on the matrix) 

      None     = no keys held down while clicking
      Shft     = hold SHIFT down while clicking
      Ctrl     = hold CTRL down while clicking
      CtSh     = hold both SHIFT and CTRL down while clicking
      
   Visualization Functions
   
      Rota     = Rotates camera about X, Y, and Z axes
      RotZ     = Rotates camera about the Z axis
      Move     = Translates along the X and Y axes
      MovZ     = Translates along Z axis
      Clip     = Y motion moves the near clipping plane while
      PkAt     = Pick an atom
      PkBd     = Pick a bond
      Orig     = Move origin to selected atom
      +lb      = Add an atom into the (lb) selection
      lb       = Define the (lb) selection with the indicated atom.
      rb       = Define the (rb) selection with the indicated atom.

   Editing Functions
   
      RotF     = Rotate fragment
      MovF     = Move fragment
      TorF     = Torsion fragment
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
   -i   Disable the internal OpenGL GUI (object list, menus, etc.)
   -x   Disable the external GUI module.
   -t   Use Tcl/Tk based external GUI module (pmg_tk).
   -w   Use wxPython based external GUI module (pmg_wx).
   -q   Quiet launch. Suppress splash screen.
   -p   Listen for commands on standard input.
   -e   Start in full-screen mode
   
   -f <# line> Controls display of commands and feedback in OpenGL (0=off).
   -r <file.py>[,global|local|module] Run a python program on startup.
   -l <file.py>[,global|local|module] Spawn a python program in new thread.
   -d <string> Run pymol command string upon startup.
   -u <script> Load and append to this PyMOL script or program file.
   -s <script> Save commands to this PyMOL script or program file.
   
   <file> can have one of the following extensions, and all 
   files provided will be loaded or run after PyMOL starts.
    
    .pml            PyMOL command script to be run on startup
    .py, .pym, .pyc Python program to be run on startup
    .pdb            Protein Data Bank format file to be loaded on startup
    .mmod           Macromodel format to be loaded on startup
    .mol            MDL MOL file to be loaded on startup
    .xplor          X-PLOR Map file (ASCII) to be loaded on startup
    .ccp4           CCP4 map file (BINARY) to be loaded on startup
    .pkl            Pickled ChemPy Model (class "chempy.model.Indexed")
    .r3d            Raster3D Object
    .cc1, .cc2      ChemDraw 3D cartesian coordinate file
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
 
   Selections are enclosed in parentheses and contain predicates,
   logical operations, object names, selection names and nested
   parenthesis: ( [... [(...) ... ]] )
 
      name <atom names>            n;<atom names>          
      resn <residue names>         r;<residue names>
      resi <residue identifiers>   i;<residue identifiers>
      chain <chain ID>             c;<chain identifiers>
      segi <segment identifiers>   s;<segment identifiers>
      elem <element symbol>        e;<element symbols>
      flag <number>                f;<number>
      alt <code>                   
      numeric_type <numeric type>  nt;<numeric type>
      text_type <text type>        tt;<text type>
      b <operator> <value>         
      q <operator> <value>         
      formal_charge <op> <value>   fc;<operator> <value>
      partial_charge <op> <value>  pc;<operator> <value>
      id <original-index>          
      hydrogen                     h;
      all                          *
      visible                      v;
      hetatm                       
      <selection> and <selection>  <selection>&<selection>
      <selection> or <selection>   <selection>|<selection>
      not <selection>              !<selection>
      byres <selection>            br;<selection>
      byobj <selection>            bo;<selection>
      around <distance>            a;<distance>
      expand <distance>            e;<distance>
      gap <distance>               
      in <selection>               
      like <selection>             l;<selection>
   '''
   help('selections')

def sort(object=""):
   '''
DESCRIPTION
 
   "sort" reorders atoms in the structure.  It usually only necessary
   to run this routine after an "alter" command which has modified the
   names of atom properties.  Without an argument, sort will resort
   all atoms in all objects.

USAGE
 
   sort [object]

PYMOL API

   cmd.sort(string object)

SEE ALSO

   alter
'''
   try:
      lock()
      r = _cmd.sort(str(object))
   finally:
      unlock()
   return r

def get_setting_updates(): # INTERNAL
   r = []
   if lock_attempt():
      try:
         r = _cmd.get_setting_updates()
      finally:
         unlock()
#   try:
#      lock()
#      r = _cmd.get_setting_updates()
#   finally:
#      unlock()
   return r

def log_open(fname='log.pml',mode='w'):
   try:
      try:
         if pymol._log_file!=None:
            pymol._log_file.close()
      except:
         pass
      pymol._log_file = open(fname,mode)
      if _feedback(fb_module.cmd,fb_mask.details): # redundant
         if mode!='a':
            print " Cmd: logging to '%s'."%fname
         else:
            print " Cmd: appending to '%s'."%fname            
      if(re.search(r"\.py$|\.PY$|\.pym$|.PYM$",fname)):
         set("logging",2,quiet=1)
      else:
         set("logging",1,quiet=1)
   except:
      print"Error: unable to open log file '%s'"%fname
      pymol._log_file = None
      set("logging",0,quiet=1)

def log(text,alt_text=None):
   if pymol._log_file!=None:
      mode = get_setting_legacy("logging")
      if mode:
         if mode==1:
            pymol._log_file.write(text)
         elif mode==2:
            if alt_text!=None:
               pymol._log_file.write(alt_text)
            else:
               pymol._log_file.write("cmd.do('''%s''')\n"%string.strip(text))
         pymol._log_file.flush()
         

_log = log # alias for set command which has local argument "log"

def log_close():
   if pymol._log_file!=None:
      pymol._log_file.close()
      set("logging",0,quiet=1)
      if _feedback(fb_module.cmd,fb_mask.details): # redundant
         print " Cmd: log closed."
      
def align(source,target): # EXPERIMENTAL, BUGGY
   r = None
   try:
      lock()
      r = _cmd.align(str(source),str(target))
   finally:
      unlock()
   return r

def transform_object(name,matrix,state=0,log=0,sele=''):
   r = None
   try:
      lock()
      r = _cmd.transform_object(str(name),int(state)-1,list(matrix),int(log),str(sele))
   finally:
      unlock()
   return r

def translate_atom(sele1,v0,v1,v2,state=0,mode=0,log=0):
   r = None
   sele1 = selector.process(sele1)
   try:
      lock()
      r = _cmd.translate_atom(str(sele1),float(v0),float(v1),float(v2),int(state)-1,int(mode),int(log))
   finally:
      unlock()
   return r

def get_setting_tuple(name,object='',state=0): # INTERNAL
   r = None
   if is_string(name):
      i = setting._get_index(name)
   else:
      i = int(name)
   if i<0:
      print "Error: unknown setting"
      raise QuietException
   try:
      lock()
      r = _cmd.get_setting_tuple(i,str(object),int(state)-1)
   finally:
      unlock()
   return r

def get_setting_text(name,object='',state=0):  # INTERNAL
   r = None
   if is_string(name):
      i = setting._get_index(name)
   else:
      i = int(name)
   if i<0:
      print "Error: unknown setting"
      raise QuietException
   try:
      lock()
      r = _cmd.get_setting_text(i,str(object),int(state)-1)
   finally:
      unlock()
   return r

def get_renderer():  # 
   try:
      lock()
      r = _cmd.get_renderer()
   finally:
      unlock()
   return r

def focus():  # BROKEN
   try:
      lock()
      r = _cmd.focus()
   finally:
      unlock()
   return r

def spheroid(object=""):  # EXPERIMENTAL
   try:
      print "Warning: 'spheroid' is experimental, incomplete, and unstable."
      lock()
      r = _cmd.spheroid(str(object))
   finally:
      unlock()
   return r

def cls(): 
   '''
DESCRIPTION
 
   "cls" clears the output buffer.

USAGE
 
   cls
'''
   r = None
   try:
      lock()
      r = _cmd.cls()
   finally:
      unlock()
   return r

def _adjust_coord(a,i,x):
   a.coord[i]=a.coord[i]+x
   return None

def fragment(name,object=None,origin=1,zoom=0):
   '''
DESCRIPTION
 
   "fragment" retrieves a 3D structure from the fragment library, which is currently
   pretty meager (just amino acids).
   
USAGE
 
   fragment name

'''
   try:
      save=get_setting_legacy('auto_zoom')
      set('auto_zoom',zoom,quiet=1)
      try:
         if object==None:
            object=name
         model = fragments.get(str(name))
         la = len(model.atom)
         if la:
            mean = map(lambda x,la=la:x/la,[
               reduce(operator.__add__,map(lambda a:a.coord[0],model.atom)),
               reduce(operator.__add__,map(lambda a:a.coord[1],model.atom)),
               reduce(operator.__add__,map(lambda a:a.coord[2],model.atom))])
            position = get_position()
            for c in range(0,3):
               mean[c]=position[c]-mean[c]
               map(lambda a,x=mean[c],c=c:_adjust_coord(a,c,x),model.atom)
            mean = map(lambda x,la=la:x/la,[
               reduce(operator.__add__,map(lambda a:a.coord[0],model.atom)),
               reduce(operator.__add__,map(lambda a:a.coord[1],model.atom)),
               reduce(operator.__add__,map(lambda a:a.coord[2],model.atom))])
         load_model(model,str(object))
      finally:
         set('auto_zoom',save,quiet=1)
   except:
      print "Error: unable to load fragment %s" % name

def wizard(name):
   '''
DESCRIPTION
 
   "wizard" launches on of the built-in wizards.  There are special
   Python scripts which work with PyMOL in order to obtain direct user
   interaction and easily peform complicated tasks.
   
USAGE
 
   wizard name

PYMOL API

   cmd.wizard(string name)

EXAMPLE

   wizard distance  # launches the distance measurement wizard
'''
   import wizard
   try:
      if not sys.modules.has_key(name):
         mod_tup = imp.find_module(name,wizard.__path__)
         mod_obj = imp.load_module(name,mod_tup[0],mod_tup[1],mod_tup[2])
      else:
         mod_obj = sys.modules[name]
      if mod_obj:
         oname = string.capitalize(name)
         if hasattr(mod_obj,oname):
            wiz = apply(getattr(mod_obj,oname))
            if wiz:
               set_wizard(wiz)
               do("refresh")
   except ImportError:
      print "Error: Sorry, couldn't find the '"+name+"' wizard."

def get_phipsi(sele1="(name ca)",state=-1):
   # preprocess selections
   sele1 = selector.process(sele1)
   #   
   r = None
   try:
      lock()
      r = _cmd.get_phipsi(str(sele1),int(state)-1)
   finally:
      unlock()
   return r

def get_position():
   r = None
   try:
      lock()
      r = _cmd.get_position()
   finally:
      unlock()
   return r
   
def get_dihedral(atom1,atom2,atom3,atom4,state=1):
   # preprocess selections
   atom1 = selector.process(atom1)
   atom2 = selector.process(atom2)
   atom3 = selector.process(atom3)
   atom4 = selector.process(atom4)
   #   
   r = None
   try:
      lock()
      r = _cmd.get_dihe(str(atom1),str(atom2),str(atom3),str(atom4),int(state)-1)
   finally:
      unlock()
   return r

def set_dihedral(atom1,atom2,atom3,atom4,angle,state=1):
   # preprocess selections
   atom1 = selector.process(atom1)
   atom2 = selector.process(atom2)
   atom3 = selector.process(atom3)
   atom4 = selector.process(atom4)
   #   
   try:
      lock()
      r = _cmd.set_dihe(str(atom1),str(atom2),str(atom3),str(atom4),float(angle),int(state)-1)
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

def edit_mode(mode=None):
   try:
      lock()
      r = _cmd.get_setting("button_mode")
      r = int(r)
      if mode==None:
         if r:
            _cmd.legacy_set("button_mode","0")
         else:
            _cmd.legacy_set("button_mode","1")            
      else:
         if mode=='on':
            _cmd.legacy_set("button_mode","1")
         if mode=='off':
            _cmd.legacy_set("button_mode","0")
      config_mouse()
   finally:
      unlock()
   pass

def get_setting_legacy(name): # INTERNAL, DEPRECATED
   r = None
   try:
      lock()
      r = _cmd.get_setting(name)
   finally:
      unlock()
   return r

def resume(fname):
   if os.path.exists(fname):
      if(re.search(r"\.py$|\.PY$|\.pym$|.PYM$",fname)):
         do("run %s"%fname)
      else:
         do("@%s"%fname)
   do("log_open %s,a"%fname)
   
def config_mouse(quiet=0): # INTERNAL
   # NOTE: PyMOL automatically runs this routine upon start-up
   try:
      lock()
      r = _cmd.get_setting("button_mode")
      r = int(r)
      if not r:
         # visualization
         button('l','none','rota')
         button('m','none','move')
         button('r','none','movz')
         button('l','shft','+lbx')
         button('m','shft','-lbx')
         button('r','shft','clip')                  
         button('l','ctrl','+lb')
         button('m','ctrl','pkat')
         button('r','ctrl','pkbd')                  
         button('l','ctsh','lb')
         button('m','ctsh','orig')
         button('r','ctsh','rb')
         if not quiet:
            print " Mouse: configured for visualization."
      else:
         # editing
         button('l','none','rota')
         button('m','none','move')
         button('r','none','movz')
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

def distance(name=None,selection1="(lb)",selection2="(rb)",cutoff=None,mode=None):
   '''
DESCRIPTION
 
   "distance" creates a new distance object between two
   selections.  It will display all distances within a cutoff.
 
USAGE
 
   distance 
   distance (selection1), (selection2)
   distance name = (selection1), (selection1) [,cutoff [,mode] ]
 
   name = name of distance object 
   selection1,selection2 = atom selections
   cutoff = maximum distance to display
   mode = 0 (default)

PYMOL API
 
   cmd.distance( string name, string selection1, string selection2,
          string cutoff, string mode )
   returns the average distance between all atoms/frames
 
NOTES

   The distance wizard makes measuring distances easier than using
   the "dist" command for real-time operations.

   "dist" alone will show distances between selections (lb) and (rb)
   created by left and right button atom picks.  CTRL-SHIFT/left-click
   on the first atom,  CTRL-SHIFT/right-click on the second, then run
   "dist".

'''
   # handle unnamed distance 

   if name!=None:
      if len(name):
         if name[0]=='(': # we're one argument off...
            if cutoff!=None:
               mode = cutoff
            if selection2!="(rb)":
               cutoff = selection2
            if selection1!="(lb)":
               selection2 = selection1
            selection1=name
            name = None

   # in unlabeled, then get next name in series
   
   if name!=None:
      nam=name
   else:
      try:
         lock()
         cnt = _cmd.get("dist_counter") + 1.0
         _cmd.legacy_set("dist_counter","%1.0f" % cnt)
         nam = "dist%02.0f" % cnt
      finally:
         unlock()

   # defaults
   if mode == None:
      mode = 0
   if cutoff == None:
      cutoff = -1.0
   # preprocess selections
   selection1 = selector.process(selection1)
   selection2 = selector.process(selection2)
   # now do the deed
   try:
      lock()
      r = _cmd.dist(str(nam),str(selection1),str(selection2),int(mode),float(cutoff))
   finally:
      unlock()
   return r

# LEGACY support for cmd.dist
def dist(*arg,**kw):
   return apply(distance,arg,kw)
#

def get_view(output=1):
   '''
DESCRIPTION
 
   "get_view" returns and optionally prints out the current view
   information in a format which can be embedded into a command
   script and used in subsequent calls to "set_view"
 
USAGE

   get_view
   
PYMOL API

   cmd.get_view(output=1)  
 
API USAGE

   cmd.get_view(0) # zero option suppresses output
'''
   
   r = None
   try:
      lock()
      r = _cmd.get_view()
   finally:
      unlock()
   if len(r):
      if get_setting_legacy("logging")!=0.0:
         print " get_view: matrix written to log file."
         log("_ set_view (\\\n","cmd.set_view((\\\n")
         log("_  %14.9f, %14.9f, %14.9f,\\\n"%r[0:3]  , "  %14.9f, %14.9f, %14.9f,\\\n"%r[0:3])
         log("_  %14.9f, %14.9f, %14.9f,\\\n"%r[4:7]  , "  %14.9f, %14.9f, %14.9f,\\\n"%r[4:7])
         log("_  %14.9f, %14.9f, %14.9f,\\\n"%r[8:11] , "  %14.9f, %14.9f, %14.9f,\\\n"%r[8:11])
         log("_  %14.9f, %14.9f, %14.9f,\\\n"%r[16:19], "  %14.9f, %14.9f, %14.9f,\\\n"%r[16:19])
         log("_  %14.9f, %14.9f, %14.9f,\\\n"%r[19:22], "  %14.9f, %14.9f, %14.9f,\\\n"%r[19:22]) 
         log("_  %14.9f, %14.9f, %14.9f )\n"%r[22:25] , "  %14.9f, %14.9f, %14.9f ))\n"%r[22:25])
         if output<2: # suppress if we have a log file open
            output=0
      if output:
         print "### cut below here and paste into script ###"
         print "set_view (\\"
         print "  %14.9f, %14.9f, %14.9f,\\"%r[0:3]
         print "  %14.9f, %14.9f, %14.9f,\\"%r[4:7]
         print "  %14.9f, %14.9f, %14.9f,\\"%r[8:11]
         print "  %14.9f, %14.9f, %14.9f,\\"%r[16:19]
         print "  %14.9f, %14.9f, %14.9f,\\"%r[19:22]
         print "  %14.9f, %14.9f, %14.9f )"%r[22:25]
         print "### cut above here and paste into script ###"
      
   r = r[0:3]+r[4:7]+r[8:11]+r[16:25]
   return r

def set_view(view):
   '''
DESCRIPTION
 
   "set_view" sets viewing information for the current scene,
   including the rotation matrix, position, origin of rotation,
   clipping planes, and the orthoscopic flag.
 
USAGE

   set_view (...)  where ... is 18 floating point numbers
   
PYMOL API

   cmd.set_view(string-or-sequence view)  
 
'''
   
   r = None
   if is_string(view):
      try:
         view = eval(view)
      except:
         print "Error: bad view argument; should be a sequence of 18 floats."
         raise QuietException
   if len(view)!=18:
      print "Error: bad view argument; should be a sequence of 18 floats."
      raise QuietException
   else:
      try:
         lock()
         r = _cmd.set_view((
            float(view[ 0]),float(view[ 1]),float(view[ 2]),0.0,
            float(view[ 3]),float(view[ 4]),float(view[ 5]),0.0,
            float(view[ 6]),float(view[ 7]),float(view[ 8]),0.0,
            0.0,0.0,0.0,1.0,
            float(view[ 9]),float(view[10]),float(view[11]),
            float(view[12]),float(view[13]),float(view[14]),
            float(view[15]),float(view[16]),float(view[17])))
      finally:
         unlock()
   return r

def bond(atom1="(lb)",atom2="(rb)",order=1):
   '''
DESCRIPTION
 
   "bond" creates a new bond between two selections, each of
   which should contain one atom.
 
USAGE

   bond [atom1,atom2 [,order]]
   
PYMOL API

   cmd.bond(string atom1, string atom2)
   
NOTES

   The atoms must both be within the same object.
   
   The default behavior is to create a bond between the (lb) and (rb)
   selections.

SEE ALSO

   unbond, fuse, attach, replace, remove_picked
'''
   # preprocess selections
   atom1 = selector.process(atom1)
   atom2 = selector.process(atom2)
   try:
      lock()
      r = _cmd.bond(atom1,atom2,int(order),1)
   finally:
      unlock()
   return r

def invert(selection1="(lb)",selection2="(rb)"):
   '''
DESCRIPTION

   "invert" inverts the stereo-chemistry of the atom currently picked
   for editing (pk1).  Two additional atom selections must be provided
   in order to indicate which atoms remain stationary during the
   inversion process.

USAGE

   invert (selection1),(selection2)

PYMOL API

   cmd.api( string selection1="(lb)", string selection2="(lb)" )
   
NOTE

   The invert function is usually bound to CTRL-E in editing mode.

   The default selections are (lb) and (rb), meaning that you can pick
   the atom to invert with CTRL-middle click and then pick the
   stationary atoms with CTRL-SHIFT/left-click and CTRL-SHIFT/right-
   click, then hit CTRL-E to invert the atom.

'''
   # preprocess selections
   selection1 = selector.process(selection1)
   selection2 = selector.process(selection2)
   #
   try:
      lock()
      r = _cmd.invert(str(selection1),str(selection2),0)
   finally:
      unlock()
   return r

def unbond(atom1="(lb)",atom2="(rb)"):
   '''
DESCRIPTION
 
   "unbond" removes all bonds between two selections.

USAGE

   unbond atom1,atom2
   
PYMOL API

   cmd.unbond(selection atom1="(lb)",selection atom2="(rb)")
   
SEE ALSO

   bond, fuse, remove_picked, attach, detach, replace
 
'''
   # preprocess selections
   atom1 = selector.process(atom1)
   atom2 = selector.process(atom2)   
   try:
      lock()
      r = _cmd.bond(str(atom1),str(atom2),0,0)
   finally:
      unlock()
   return r

def show_help(cmd): # INTERNAL
   print "PyMOL>help %s" % cmd
   help(cmd)
   if get_setting_legacy("internal_feedback")>0.1:
      print "(Hit ESC to hide)"

def set_wizard(*arg): # INTERNAL
   r = None
   wiz = None
   if len(arg):
      wiz=arg[0]
   try:
      lock()
      r = _cmd.set_wizard(wiz)
   finally:
      unlock()
   return r

def refresh_wizard(*arg): # INTERNAL
   r = None
   wiz = None
   if len(arg):
      wiz=arg[0]
   try:
      lock()
      r = _cmd.refresh_wizard()
   finally:
      unlock()
   return r

def get_wizard(*arg): # INTERNAL
   r = None
   try:
      lock()
      r = _cmd.get_wizard()
   finally:
      unlock()
   return r

def undo():
   '''
DESCRIPTION

   "undo" restores the previous conformation of the object currently
   being edited.
   
USAGE
 
   undo

SEE ALSO

   redo, push_undo
'''
   r = None
   try:
      lock()
      r = _cmd.undo(-1)
   finally:
      unlock()
   return r

def push_undo(selection,state=0):
   '''
DESCRIPTION

   "push_undo" stores the currently conformations of objects in the
   selection onto their individual kill rings.
   
USAGE
 
   push_undo (all)

SEE ALSO

   undo, redo
'''
   # preprocess selections
   selection = selector.process(selection)
   #
   r = None
   try:
      lock()
      r = _cmd.push_undo(str(selection),int(state)-1)
   finally:
      unlock()
   return r

def redo():
   '''
DESCRIPTION

   "redo" reapplies the conformational change of the object currently
   being edited.
   
USAGE
 
   redo

SEE ALSO

   undo, push_undo
'''
   try:
      lock()
      _cmd.undo(1)
   finally:
      unlock()

def help(command = "commands"):
   '''
DESCRIPTION

   "help" prints out the online help for a given command.
   
USAGE
 
   help command
   '''
   if get_setting_legacy("internal_feedback")>0.1:
      set("text","1",quiet=1)
   cmd = help_sc.auto_err(command,'topics')   
   if keyword.has_key(cmd):
      doc = keyword[cmd][0].__doc__
      if doc:
         print "\n",string.strip(doc),"\n"
      else:
         print "Error: sorry no help available on that command."
   elif help_only.has_key(cmd):
      doc = help_only[cmd][0].__doc__
      if doc:
         print "\n",string.strip(doc),"\n"
      else:
         print "Error: sorry no help available on that command."      
   else:
      print "Error: unrecognized command"

def symexp(prefix,object,selection,cutoff):
   '''
DESCRIPTION
 
   "symexp" creates all symmetry related objects for the specified object
   that occurs within a cutoff about an atom selection.  The new objects
   are labeled using the prefix provided along with their crystallographic
   symmetry operation and translation.
 
USAGE
 
   symexp prefix = object, (selection), cutoff
 
PYMOL API
 
   cmd.symexp( string prefix, string object, string selection, float cutoff) 

SEE ALSO

   load
   '''
   # preprocess selection
   selection=selector.process(selection)
   #
   try:
      lock()
      r = _cmd.symexp(str(prefix),str(object),str(selection),float(cutoff))
   finally:
      unlock()
   return r

def map_set_border(name,level=0.0):
   '''
DESCRIPTION
 
   "map_set_border" is a function (reqd by PDA) which allows you to set the
   level on the edge points of a map

USAGE

   map_set_border <name>,<level>

NOTES

   unsupported.
SEE ALSO

   load
   '''
   try:
      lock()
      r = _cmd.map_set_border(str(name),float(level))
   finally:
      unlock()
   return r

def isomesh(name,map,level=1.0,selection='',buffer=0.0,state=0,carve=None):
   '''
DESCRIPTION
 
   "isomesh" creates a mesh isosurface object from a map object.
 
USAGE
 
   isomesh name, map, level [,(selection) [,buffer [,state [,carve ]]]]

   "name" is the name for the new mesh isosurface object.
   
   "map" is the name of the map object to use for computing the mesh.
   
   "level" is the contour level.

   "selection" is an atom selection about which to display the mesh with
      an additional "buffer" (if provided).

   "state" is the state into which the object should be loaded.

   "carve" is a radius about each atom in the selection for which to
      include density. If "carve" is not provided, then the whole
      brick is displayed.

NOTES

   If the mesh object already exists, then the new mesh will be
   appended onto the object as a new state (unless you indicate a state).

SEE ALSO

   isodot, load
   '''
   if selection!='':
      mopt = 1 # about a selection
   else:
      mopt = 0 # render the whole map
   # preprocess selection
   selection = selector.process(selection)
   #
   if carve==None:
      carve=-1.0
   try:
      lock()
      r = _cmd.isomesh(str(name),0,str(map),int(mopt),
                       str(selection),float(buffer),
                       float(level),0,int(state)-1,float(carve))
   finally:
      unlock()
   return r

def isodot(name,map,level=1.0,selection='',buffer=0.0,state=0,carve=0):
   '''
DESCRIPTION
 
"isodot" creates a dot isosurface object from a map object.
 
USAGE
 
   isodot name = map, level [,(selection) [,buffer [, state ] ] ] 

   "map" is the name of the map object to use.
   
   "level" is the contour level.
   
   "selection" is an atom selection about which to display the mesh with
      an additional "buffer" (if provided).

NOTES

   If the dot isosurface object already exists, then the new dots will
   be appended onto the object as a new state.

SEE ALSO

   load, isomesh
   '''
   if selection!='':
      mopt = 1 # about a selection
   else:
      mopt = 0 # render the whole map
   # preprocess selections
   selection = selector.process(selection)
   #
   try:
      lock()
      r = _cmd.isomesh(str(name),0,str(map),int(mopt),
                       str(selection),float(buffer),
                       float(level),1,int(state)-1,float(carve))
   finally:
      unlock()
   return r

def ready(): # INTERNAL
   return _cmd.ready()

def splash():
   '''
DESCRIPTION
 
   "splash" shows the splash screen information.

USAGE

   splash
   '''
   if get_setting_legacy("internal_feedback")>0.1:
      set("text","1",quiet=1)
   print
   try:
      lock()
      r = _cmd.splash()
   finally:
      unlock()
   return r

def bg_color(color="black"):
   '''
DESCRIPTION
 
   "bg_color" sets the background color

USAGE

   bg_color [color]

PYMOL API

   cmd.color(string color="black")
   
   '''
   try:
      lock()
      r = _cmd.bg_color(str(color))
   finally:
      unlock()
   return r

def copy(target,source):
   '''
DESCRIPTION
 
   "copy" creates a new object that is an identical copy of an
   existing object
 
USAGE

   copy target, source

   copy target = source         # (DEPRECATED)
 
PYMOL API
 
   cmd.copy(string target,string source)

SEE ALSO

   create
   '''
   try:
      lock()
      r = _cmd.copy(str(source),str(target))
   finally:
      unlock()
   return r

def label(selection="(all)",expression=""):
   '''
DESCRIPTION
 
   "label" labels one or more atoms properties over a selection using
   the python evaluator with a separate name space for each atom.  The
   symbols defined in the name space are:
 
      name, resn, resi, chain, q, b, segi, type (ATOM,HETATM) 
      formal_charge, partial_charge, numeric_type, text_type
   
   All strings in the expression must be explicitly quoted.  This
   operation typically takes several seconds per thousand atoms
   altered.

   To clear labels, simply omit the expression or set it to ''.
   
USAGE

   label (selection),expression
   
EXAMPLES
  
   label (chain A),chain
   label (n;ca),"%s-%s" % (resn,resi)
   label (resi 200),"%1.3f" % partial_charge
   '''
   # preprocess selection
   selection = selector.process(selection)
   #
   r = 1
   try:
      lock()
      if len(str(expression))==0:
         r= _cmd.label(str(selection),'')
      else:
         r = _cmd.label(str(selection),'label='+str(expression))
   finally:
      unlock()   
   return r

def alter(selection,expression):
   '''
DESCRIPTION
 
   "alter" changes one or more atomic properties over a selection
   using the python evaluator with a separate name space for each
   atom.  The symbols defined in the name space are:
 
      name, resn, resi, chain, alt, elem, q, b, segi,
      type (ATOM,HETATM), partial_charge, formal_charge,
      text_type, numeric_type, ID
   
   All strings must be explicitly quoted.  This operation typically
   takes several seconds per thousand atoms altered.

   WARNING: You should always issue a "sort" command on an object
   after modifying any property which might affect canonical atom
   ordering (names, chains, etc.).  Failure to do so will confound
   subsequent "create" and "byres" operations.
   
USAGE
 
   alter (selection),expression
 
EXAMPLES
 
   alter (chain A),chain='B'
   alter (all),resi=str(int(resi)+100)
   sort
   
SEE ALSO

   alter_state, iterate, iterate_state, sort
   '''
   # preprocess selections
   selection = selector.process(selection)
   #
   try:
      lock()
      r = _cmd.alter("("+str(selection)+")",str(expression),0)
   finally:
      unlock()   
   return r


def iterate(selection,expression):
   '''
DESCRIPTION
 
   "iterate" iterates over an expression with a separate name space
   for each atom.  However, unlike the "alter" command, atomic
   properties can not be altered.  Thus, "iterate" is more efficient
   than "alter".

   It can be used to perform operations and aggregations using atomic
   selections, and store the results in any global object, such as the
   predefined "stored" object.

   The local namespace for "iterate" contains the following names

      name, resn, resi, chain, alt, elem,
      q, b, segi, and type (ATOM,HETATM),
      partial_charge, formal_charge,
      text_type, numeric_type, ID
 
   All strings in the expression must be explicitly quoted.  This
   operation typically takes a second per thousand atoms.
 
USAGE
 
   iterate (selection),expression
 
EXAMPLES

   stored.net_charge = 0
   iterate (all),stored.net_charge = stored.net_charge + partial_charge

   stored.names = []
   iterate (all),stored.names.append(name)

SEE ALSO

   iterate_state, atler, alter_state
   '''
   # preprocess selection
   selection = selector.process(selection)
   #
   try:
      lock()
      r = _cmd.alter("("+str(selection)+")",str(expression),1)
   finally:
      unlock()   
   return r

def alter_state(state,selection,expression):
   '''
DESCRIPTION
 
   "alter_state" changes the atomic coordinates of a particular state
   using the python evaluator with a separate name space for each
   atom.  The symbols defined in the name space are:
 
      x,y,z
 
USAGE
 
   alter_state state,(selection),expression
 
EXAMPLES
 
   alter 1,(all),x=x+5

SEE ALSO

   iterate_state, alter, iterate
   '''
   # preprocess selection
   selection = selector.process(selection)
   #
   try:
      lock()
      r = _cmd.alter_state(int(state)-1,"("+str(selection)+")",str(expression),0)
   finally:
      unlock()   
   return r


view_dict = {}
view_sc = Shortcut(['store','recall'])
view_dict_sc = Shortcut([])

def view(key,action='recall'):
   '''
DESCRIPTION
 
   "view" makes it possible to save and restore viewpoints on a given
   scene within a single session.
 
USAGE
 
   view key[,action]
   view ?

   key can be any string
   action should be 'store' or 'recall' (default: 'recall')
   
PYMOL API

   cmd.view(string key,string action)
   
EXAMPLES

   view 0,store
   view 0

SEE ALSO

   get_view
   '''
   if key=='?':
      print " view: stored views:"
      parsing.dump_str_list(view_dict.keys())
   else:
      action = view_sc.auto_err(action,'action')
      if action=='recall':
         key = view_dict_sc.auto_err(key,'view')
         set_view(view_dict[key])
         if _feedback(fb_module.scene,fb_mask.blather): # redundant
            print " view: recalled."
      elif action=='store':
         view_dict_sc.append(key)
         view_dict[key]=get_view(0)
         if _feedback(fb_module.scene,fb_mask.actions):
            print " view: view stored as '%s'."%key

def iterate_state(state,selection,expression):
   '''
DESCRIPTION
 
   "iterate_state" is to "alter_state" as "iterate" is to "alter"
 
USAGE
 
   iterate_state state,(selection),expression
 
EXAMPLES

   stored.sum_x = 0.0
   iterate 1,(all),stored.sum_x = stored.sum_x + x

SEE ALSO

   iterate, alter, alter_state
   '''
   # preprocess selection
   selection = selector.process(selection)
   #
   try:
      lock()
      r = _cmd.alter_state(int(state)-1,str(selection),str(expression),1)
   finally:
      unlock()   
   return r

def _stereo(flag): # SGI-SPECIFIC - bad bad bad
   
   if flag:
      os.system("/usr/gfx/setmon -n 1024x768_96s")
   else:
      os.system("/usr/gfx/setmon -n 72hz")

def stereo(state='on'):
   '''
DESCRIPTION
 
   "stereo" activates or deactives stereo mode.  Currently only
   high-end stereo graphics are supported on the SGI (stereo in a
   window).
 
USAGE
 
   stereo on
   stereo off

PYMOL API

   cmd.stereo(string state="on")
   '''
   state = toggle_dict[toggle_sc.auto_err(str(state),'toggle')]
   r = None
   if state:
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
         else:
            print "Error: stereo not available"
      finally:
         unlock();
   return r
   
def overlap(selection1,selection2,state1=1,state2=1,adjust=0.0):
#
#   UNSUPPORTED FEATURE - LIKELY TO CHANGE
#   (for maximum efficiency, use smaller molecule as selection 1)
#
   # preprocess selections
   selection1 = selector.process(selection1)
   selection2 = selector.process(selection2)
   #
   r = 1
   try:
      lock()
      r = _cmd.overlap(str(selection1),str(selection2),
                       int(state1)-1,int(state2),
                       float(adjust))
   finally:
      unlock()
   return r

def setup_global_locks(): # INTERNAL, OBSOLETE?
   pass

def lock_c(): # INTERNAL
   lock_api_c.acquire(1)

def unlock_c(): # INTERNAL
   lock_api_c.release()

def lock(): # INTERNAL
   lock_api.acquire(1) 

def lock_attempt(): # INTERNAL
   return lock_api.acquire(blocking=0)

def unlock(): # INTERNAL
   if (thread.get_ident() == pymol.glutThread):
      lock_api.release()
      _cmd.flush_now()
   else:
      lock_api.release()
      if _cmd.wait_queue(): # commands waiting to be executed?
         w = 0.0025 # NOTE: affects API perf. for "do" and delayed-exec
         while 1:
            e = threading.Event() # using this for portable delay
            e.wait(w)
            del e
            if not _cmd.wait_queue():
               break
            if w > 0.1: # wait up 0.2 sec max for PyMOL to flush queue
               if _feedback(fb_module.cmd,fb_mask.debugging):
                  fb_debug.write("Debug: avoiding possible dead-lock?\n")
               break
            w = w * 2 # wait twice as long each time until flushed

def export_dots(object,state):  
# UNSUPPORTED
   try:
      lock()
      r = _cmd.export_dots(object,int(state)-1)
   finally:
      unlock()
   return r

def sync(timeout=1.0,poll=0.05):
   '''
DESCRIPTION

   "sync" is an API-only function which waits until all current
   commmands have been executed before returning.  A timeout
   can be used to insure that this command eventually returns.

PYMOL API

   cmd.sync(float timeout=1.0,float poll=0.05)

SEE ALSO

   frame
'''
   now = time.time()
   timeout = float(timeout)
   poll = float(poll)
   if _cmd.wait_queue(): # commands waiting to be executed?
      while 1:
         e = threading.Event() # using this for portable delay
         e.wait(poll)
         del e
         if not _cmd.wait_queue():
            break
         if (time.time()-now)>timeout:
            break
            
def count_states(selection="(all)"):
   '''
DESCRIPTION

   "count_states" is an API-only function which returns the number of
   states in the selection.

PYMOL API

   cmd.count_states(string selection="(all)")

SEE ALSO

   frame
'''
   # preprocess selection
   selection = selector.process(selection)
   #
   try:
      lock()
      r = _cmd.count_states(selection)
   finally:
      unlock()
   return r

def do(commands):
   # WARNING: don't call this routine if you already have the API lock
   '''
DESCRIPTION
 
   "do" makes it possible for python programs to issue simple PyMOL
   commands as if they were entered on the command line.
    
PYMOL API
 
   cmd.do( commands )
 
USAGE (PYTHON)
 
   from pymol import cmd
   cmd.do("load file.pdb")
   '''
   lst = string.split(commands,"\n")   
   try:
      lock()
      for a in lst:
         if(len(a)):
            r = _cmd.do(a)
   finally:
      unlock()
   return r

def turn(axis,angle):
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

SEE ALSO

   move
   '''
   try:
      lock()
      r = _cmd.turn(str(axis),float(angle))
   finally:
      unlock()
   return r

def ray(width=0,height=0):
   '''
DESCRIPTION
  
   "ray" creates a ray traced image of the current frame. This
   can take some time (up to several minutes, depending on image
   complexity).
      
USAGE
 
   ray [width,height]
 
PYMOL API
  
   cmd.ray(int width,int height)

SEE ALSO

   "help faster" for optimization tips
   '''
   try:
      lock()   
      r = _cmd.render(int(width),int(height))
   finally:
      unlock()
   return r

def faster():
   '''
RAY TRACING OPTIMIZATION

   1. Reduce object complexity to a minimum acceptable level.
         For example, try lowering:
            "cartoon_sampling" 
            "ribbon_sampling", and
            "surface_quality", as appropriate.

   2. Increase "hash_max" so as to obtain a voxel dimensions of
      0.3-0.6.  Proper tuning of "hash_max" can speed up
      rendering by a factor of 2-3 for non-trivial scenes.
      
      WARNING: memory usage depends on hash_max^3, so avoid
      pushing into virtual memory.  Roughly speaking:
      
         hash_max = 80  -->   ~9 MB hash + data
         hash_max = 160 -->  ~72 MB hash + data
         hash_max = 240 --> ~243 MB hash + data

      Avoid using virtual memory.
      
   3. Recompiling with optimizations on usually gives a 25-33%
      performance boost for ray tracing.
   
'''
   help('faster')

def system(command):
   '''
DESCRIPTION

   "system" executes a command in a subshell under Unix or Windows.

USAGE

   system command

PYMOL API

   cmd.system(string command)

SEE ALSO

   ls, cd, pwd
   '''
   r = _cmd.system(str(command))
   return r

def intra_fit(selection,state=0):
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

SEE ALSO

   fit, rms, rms_cur, intra_rms, intra_rms_cur, pair_fit
   '''
   # preprocess selection
   selection = selector.process(selection)
   #   
   r = -1.0
   try:
      lock()
      r = _cmd.intrafit(str(selection),int(state)-1,2)
   finally:
      unlock()
   return r

def intra_rms(selection,state=0):
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

SEE ALSO

   fit, rms, rms_cur, intra_fit, intra_rms_cur, pair_fit
   '''
   # preprocess selection
   selection = selector.process(selection)
   #   
   r = -1.0
   try:
      lock()
      r = _cmd.intrafit(str(selection),int(state)-1,1)
   finally:
      unlock()
   return r

def intra_rms_cur(selection,state=0):
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

SEE ALSO

   fit, rms, rms_cur, intra_fit, intra_rms, pair_fit
   '''
   # preprocess selection
   selection = selector.process(selection)
   #   
   r = -1.0
   try:
      lock()
      r = _cmd.intrafit(str(selection),int(state)-1,0)
   finally:
      unlock()
   return r

def update(target,source):
   '''
DESCRIPTION
  
   "update" transfers coordinates from one selection to another.
USAGE
 
   update (target-selection),(source-selection)
 
EXAMPLES
 
   update target,(variant)

NOTES

   Currently, this applies across all pairs of states.  Fine
   control will be added later.

SEE ALSO

   load
'''
   
   a=target
   b=source
   # preprocess selections
   a = selector.process(a)
   b = selector.process(b)
   #
   if a[0]!='(': a="("+str(a)+")"
   if b[0]!='(': b="("+str(b)+")"   
   try:
      lock()   
      r = _cmd.update(str(a),str(b),-1,-1)
   finally:
      unlock()
   return r
   
def fit(selection,target):
   '''
DESCRIPTION
  
   "fit" superimposes the model in the first selection on to the model
   in the second selection.  Only matching atoms in both selections
   will be used for the fit.
   
USAGE
 
   fit (selection), (target-selection)
 
EXAMPLES
 
   fit ( mutant and name ca ), ( wildtype and name ca )

SEE ALSO

   rms, rms_cur, intra_fit, intra_rms, intra_rms_cur
   '''
   a=str(selection)
   b=str(target)
   # preprocess selections
   a = selector.process(a)
   b = selector.process(b)
   #
   if a[0]!='(': a="("+a+")"
   if b[0]!='(': b="("+b+")"
   try:
      lock()   
      r = _cmd.fit("(%s in %s)" % (str(a),str(b)),
                  "(%s in %s)" % (str(b),str(a)),2)
   finally:
      unlock()
   return r

def rms(selection,target):
   '''
DESCRIPTION
  
   "rms" computes a RMS fit between two atom selections, but does not
   tranform the models after performing the fit.
   
USAGE
 
   rms (selection), (target-selection)
 
EXAMPLES
 
   fit ( mutant and name ca ), ( wildtype and name ca )

SEE ALSO

   fit, rms_cur, intra_fit, intra_rms, intra_rms_cur, pair_fit   
   '''
   a=str(selection)
   b=str(target)
   # preprocess selections
   a = selector.process(a)
   b = selector.process(b)
   #
   if a[0]!='(': a="(%"+a+")"
   if b[0]!='(': b="(%"+b+")"
   try:
      lock()   
      r = _cmd.fit("(%s in %s)" % (str(a),str(b)),
                  "(%s in %s)" % (str(b),str(a)),1)
   finally:
      unlock()
   return r

def rms_cur(selection,target):
   '''
DESCRIPTION
  
   "rms_cur" computes the RMS difference between two atom
   selections without performing any fitting.
   
USAGE
 
   rms_cur (selection), (selection)

SEE ALSO

   fit, rms, intra_fit, intra_rms, intra_rms_cur, pair_fit   
   '''
   a=str(selection)
   b=str(target)
   # preprocess selections
   a = selector.process(a)
   b = selector.process(b)
   #   
   if a[0]!='(': a="(%"+a+")"
   if b[0]!='(': b="(%"+b+")"
   try:
      lock()   
      r = _cmd.fit("(%s in %s)" % (str(a),str(b)),
                  "(%s in %s)" % (str(b),str(a)),0)
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

SEE ALSO

   fit, rms, rms_cur, intra_fit, intra_rms, intra_rms_cur
   '''
   new_arg = []
   for a in arg:
      new_arg.append(selector.process(a))
   try:
      lock()   
      r = _cmd.fit_pairs(new_arg)
   finally:
      unlock()
   return r

def expfit(a,b): # Huh?
   try:
      lock()   
      r = _cmd.fit(a,b,2)
   finally:
      unlock()
   return r


cartoon_dict = {
   'skip'        : -1,
   'automatic'   : 0,
   'loop'        : 1,
   'rectangle'   : 2,
   'oval'        : 3,
   'tube'        : 4,
   'arrow'       : 5,
   'dumbbell'    : 6,   
}

cartoon_sc = Shortcut(cartoon_dict.keys())

def cartoon(type,selection="(all)"):
   '''
DESCRIPTION
  
   "cartoon" changes the default cartoon for a set of atoms.
      
USAGE
 
   cartoon type, (selection)

   type = skip | automatic | loop | rectangle | oval | tube | arrow | dumbbell

PYMOL API
  
   cmd.cartoon(string type, string selection )

EXAMPLES

   cartoon rectangle,(chain A)
   cartoon skip,(resi 145:156)
   
NOTES

   the "automatic" mode utilizes ribbons according to the
   information in the PDB HELIX and SHEET records.

'''
   # preprocess selection
   selection = selector.process(selection)
   #
   type = cartoon_dict[cartoon_sc.auto_err(str(type),'type')];
   r = 1
   try:
      lock()   
      r = _cmd.cartoon(str(selection),int(type))
   finally:
      unlock()
   return r

def remove(selection):
   '''
DESCRIPTION
  
   "remove" eleminates a selection of atoms from models.
      
USAGE
 
   remove (selection)
 
PYMOL API
  
   cmd.remove( string selection )
    
EXAMPLES
 
   remove ( resi 124 )

SEE ALSO

   delete
'''
   # preprocess selection
   selection = selector.process(selection)
   #   
   r = 1
   try:
      lock()   
      r = _cmd.remove(selection)
   finally:
      unlock()
   return r


def remove_picked(hydrogens=1):
   '''
DESCRIPTION
  
   "remove_picked" removes the atom or bond currently
   picked for editing. 
      
USAGE
 
   remove_picked [hydrogens]
 
PYMOL API
  
   cmd.remove_picked(integer hydrogens=1)

NOTES

   This function is usually connected to the
   DELETE key and "CTRL-D".
   
   By default, attached hydrogens will also be deleted unless
   hydrogen-flag is zero.

SEE ALSO

   attach, replace
'''
   r = 1
   try:
      lock()   
      r = _cmd.remove_picked(int(hydrogens))
   finally:
      unlock()
   return r

def cycle_valence(h_fill=1):
   '''
DESCRIPTION
  
   "cycle_valence" cycles the valence on the currently selected bond.
      
USAGE
 
   cycle_valence [ h_fill ]
 
PYMOL API
  
   cmd.cycle_valence(int h_fill)

EXAMPLES

   cycle_valence
   cycle_valence 0
   
NOTES

   If the h_fill flag is true, hydrogens will be added or removed to
   satisfy valence requirements.
   
   This function is usually connected to the DELETE key and "CTRL-W".

SEE ALSO

   remove_picked, attach, replace, fuse, h_fill
'''
   r = 1
   try:
      lock()   
      r = _cmd.cycle_valence()
   finally:
      unlock()
   if h_fill:
      globals()['h_fill']()
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

   Immature functionality.  See code for details.

'''
   r = 1
   try:
      lock()   
      r = _cmd.attach(str(name),int(geom),int(valence))
   finally:
      unlock()
   return r

def fuse(selection1="(lb)",selection2="(rb)"):
   '''
DESCRIPTION
  
   "fuse" joins two objectss into one by forming a bond.  A copy of
   the object containing the second atom is moved so as to form an
   approximately resonable bond with the first, and is then merged
   with the first object.
      
USAGE
 
   fuse (selection1), (selection2)
 
PYMOL API
  
   cmd.fuse( string selection1="(lb)", string selection2="(lb)" )

NOTES

   Each selection must include a single atom in each object.
   The atoms can both be hydrogens, in which case they are
   eliminated, or they can both be non-hydrogens, in which
   case a bond is formed between the two atoms.

SEE ALSO

   bond, unbond, attach, replace, fuse, remove_picked
'''
   # preprocess selections
   selection1 = selector.process(selection1)
   selection2 = selector.process(selection2)
   #   
   try:
      lock()
      r = _cmd.fuse(str(selection1),str(selection2))
   finally:
      unlock()
   return r

def unpick(*arg):
   '''
DESCRIPTION

   "unpick" deletes the special "pk" atom selections (pk1, pk2, etc.)
   used in atom picking and molecular editing.

USAGE

   unpick

PYMOL API

   cmd.unpick()

SEE ALSO

   edit
   '''
   try:
      lock()   
      r = _cmd.unpick()
   finally:
      unlock()
   return r
   
def edit(selection1='',selection2='',selection3='',selection4='',pkresi=0):
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

SEE ALSO

   unpick, remove_picked, cycle_valence, torsion
'''
   # preprocess selections
   selection1 = selector.process(selection1)
   selection2 = selector.process(selection2)
   selection3 = selector.process(selection3)
   selection4 = selector.process(selection4)
   #
   r = 1
   try:
      lock()   
      r = _cmd.edit(str(selection1),str(selection2),
                    str(selection3),str(selection4),int(pkresi))
   finally:
      unlock()
   return r

def torsion(angle):
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

SEE ALSO

   edit, unpick, remove_picked, cycle_valence
'''
   try:
      lock()   
      r = _cmd.torsion(float(angle))
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

SEE ALSO

   edit, cycle_valences, h_add
'''
   r = 1
   try:
      lock()   
      r = _cmd.h_fill()
   finally:
      unlock()
   return r

def h_add(selection="(all)"):
   '''
DESCRIPTION
  
   "h_add" uses a primitive algorithm to add hydrogens
   onto a molecule.
      
USAGE
 
   h_add (selection)
 
PYMOL API
  
   cmd.h_add( string selection="(all)" )

SEE ALSO

   h_fill
'''
   # preprocess selection
   selection = selector.process(selection)
   #   
   r = 1
   try:
      lock()   
      r = _cmd.h_add(selection)
   finally:
      unlock()
   return r
   
def protect(selection="(all)"):
   '''
DESCRIPTION

   "protect" protects a set of atoms from tranformations performed
   using the editing features.  This is most useful when you are
   modifying an internal portion of a chain or cycle and do not wish
   to affect the rest of the molecule.
   
USAGE

   protect (selection)

PYMOL API

   cmd.protect(string selection)
   
SEE ALSO

   deprotect, mask, unmask, mouse, editing
'''
   # preprocess selection
   selection = selector.process(selection)
   #   
   try:
      lock()   
      r = _cmd.protect(str(selection),1)
   finally:
      unlock()
   return r

def deprotect(selection="(all)"):
   '''
DESCRIPTION

   "deprotect" reveres the effect of the "protect" command.

USAGE

   deprotect (selection)
   
PYMOL API

   cmd.deprotect(string selection="(all)")

SEE ALSO

   protect, mask, unmask, mouse, editing
'''
   # preprocess selection
   selection = selector.process(selection)
   #   
   try:
      lock()   
      r = _cmd.protect(str(selection),0)
   finally:
      unlock()
   return r

def mask(selection="(all)"):
   '''
DESCRIPTION
  
   "mask" makes it impossible to select the indicated atoms using the
   mouse.  This is useful when you are working with one molecule in
   front of another and wish to avoid accidentally selecting atoms in
   the background.

USAGE

   mask (selection)

PYMOL API

   cmd.mask( string selection="(all)" )
   
SEE ALSO

   unmask, protect, deprotect, mouse
'''
   # preprocess selection
   selection = selector.process(selection)
   #
   try:
      lock()   
      r = _cmd.mask(str(selection),1)
   finally:
      unlock()
   return r

def unmask(selection="(all)"):
   '''
DESCRIPTION

   "unmask" reverses the effect of "mask" on the indicated atoms.

PYMOL API

   cmd.unmask( string selection="(all)" )
   
USAGE

   unmask (selection)

SEE ALSO

   mask, protect, deprotect, mouse
'''
   # preprocess selection
   selection = selector.process(selection)
   #   
   try:
      lock()   
      r = _cmd.mask(str(selection),0)
   finally:
      unlock()
   return r

def replace(name,geometry,valence):
   '''
DESCRIPTION
  
   "replace" replaces the picked atom with a new atom.
      
USAGE
 
   replace name, geometry, valence
 
PYMOL API
  
   cmd.replace(string name, int geometry,int valence )

NOTES

   Immature functionality. See code for details.

SEE ALSO

   remove, attach, fuse, bond, unbond
'''
   r = 1
   try:
      lock()
      r = _cmd.replace(str(name),int(geometry),int(valence))
   finally:
      unlock()
   return r

def zoom(selection="all",buffer=0.0):
   '''
DESCRIPTION
  
   "zoom" scales and translates the window and the origin to cover the
   atom selection.
      
USAGE
 
   zoom object-or-selection [,buffer]
   zoom (selection) [,buffer]
 
PYMOL API

   cmd.zoom( string object-or-selection [,float buffer] )

SEE ALSO

   origin, orient
   '''
   # preprocess selection
   selection = selector.process(selection)
   #   
   try:
      lock()   
      r = _cmd.zoom(str(selection),float(buffer))
   finally:
      unlock()
   return r

def rename(object,force=0):
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

SEE ALSO

   alter
'''   
   try:
      lock()   
      r = _cmd.rename(str(object),int(force))
   finally:
      unlock()
   return r
   
def frame(frame):
   '''
DESCRIPTION
  
   "frame" sets the viewer to the indicated movie frame.
   
USAGE
 
   frame frame-number
 
PYMOL API
 
   cmd.frame( int frame_number )
 
NOTES
 
   Frame numbers are 1-based

SEE ALSO

   count_states
   '''
   try:
      lock()   
      r = _cmd.frame(int(frame))
   finally:
      unlock()
   return r

def move(axis,angle):
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

SEE ALSO

   turn
   '''
   try:
      lock()   
      r = _cmd.move(str(axis),float(angle))
   finally:
      unlock()
   return r

clip_action_sc = Shortcut([ 'near','far','move','slab' ])

def clip(mode,offset):
   '''
DESCRIPTION
  
   "clip" alterss the near and far clipping planes according to 
      
USAGE
  
   clip {near|far|move|slab}, distance
 
EXAMPLES
 
   clip near, -5  # moves near plane away from you by 5 A
   clip far, 10   # moves far plane towards you by 10 A
   clip slab, 20  # sets slab thickness to 20 A
   clip move, -5  # moves the slab away from you by 5 A
   
PYMOL API

   cmd.clip( string mode, float distance )

SEE ALSO

   zoom, reset
   '''
   mode = clip_action_sc.auto_err(str(mode),'mode')
   try:
      lock()   
      r = _cmd.clip(str(mode),float(offset))
   finally:
      unlock()
   return r

def origin(selection="(all)"):
   '''
DESCRIPTION
  
   "origin" sets the center of rotation about a selection
      
USAGE
 
   origin object-or-selection
   origin (selection)
 
PYMOL API
 
   cmd.origin( string object-or-selection )

SEE ALSO

   zoom, orient, reset
   '''
   # preprocess selection
   selection = selector.process(selection)
   #   
   try:
      lock()   
      r = _cmd.origin(selection)
   finally:
      unlock()
   return r

def orient(selection="(all)"):
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

SEE ALSO

   zoom, origin, reset
   '''
   # preprocess selection
   selection = selector.process(selection)
   #   
   try:
      lock()
      r = _cmd.orient(selection)
   finally:
      unlock()
   return r

def is_glut_thread(): # internal
   if thread.get_ident() == pymol.glutThread:
      return 1
   else:
      return 0

def refresh():
   '''
DESCRIPTION
  
   "refresh" causes the scene to be refresh as soon as it is safe to
   do so.

USAGE

   refresh

PYMOL API
 
   cmd.refresh()

SEE ALSO

   rebuild
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

def _refresh(swap_buffers=1):  # Only call with GLUT thread!
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

def dirty(): # OBSOLETE?
   try:
      lock()
      r = _cmd.dirty()
   finally:
      unlock()
   return r

def set(name,value,selection='',state=0,quiet=0,updates=1,log=0):
   '''
DESCRIPTION
  
   "set" changes one of the PyMOL state variables,
      
USAGE
 
   set name, value [,object-or-selection [,state ]]

   set name = value      # (DEPRECATED)

   WARNING: object and state specific settings are not yet fully
     implemented -- look for them in version 0.51.
 
PYMOL API
 
   cmd.set ( string name, string value,
             string selection='', int state=0,
             int quiet=0, int updates=1 )

NOTES

   The default behavior (with a blank selection) changes the global
   settings database.  If the selection is 'all', then the settings
   database in all individual objects will be changed.  Likewise, for
   a given object, if state is zero, then the object database will be
   modified.  Otherwise, the settings database for the indicated state
   within the object will be modified.

   If a selection is provided, then all objects in the selection will
   be affected. 
   
   '''
   r = None
   if log:
      _log("set %s=%s\n"%(name,value))
   index = setting._get_index(str(name))
   if(index<0):
      print "Error: unknown setting '%s'."%name
      raise QuietException
   else:
      success = 0
      try:
         lock()
         type = _cmd.get_setting_tuple(int(index),str(""),int(-1))[0]
         if type==None:
            print "Error: unable to get setting type."
            raise QuietException
         try:
            if type==1: # boolean
               v = (setting.boolean_dict[
                      setting.boolean_sc.auto_err(
                         str(value),"boolean")],)
            elif type==2: # int
               v = (int(value),)
            elif type==3: # float
               v = (float(value),)
            elif type==4: # float3 - some legacy handling req.
               if is_string(value):
                  if not ',' in value:
                     v = string.split(value)
                  else:
                     v = eval(value)
               else:
                  v = value
               v = (float(v[0]),float(v[1]),float(v[2]))
            v = (type,v)
            r = _cmd.set(int(index),v,
                         string.strip(str(selection)),
                         int(state)-1,int(quiet),
                         int(updates))
         except:
            if(_feedback(fb_module.cmd,fb_mask.debugging)):
               traceback.print_exc()
               raise QuietException
            print "Error: unable to read setting value."
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
DESCRIPTION
  
   "meter_reset" resets the frames per secound counter
   
USAGE

   meter_reset
   '''
   try:
      lock()   
      r = _cmd.reset_rate()
   finally:
      unlock()
   return r

def delete(name):
   '''
DESCRIPTION
  
   "delete" removes an object or a selection. 
   
USAGE
 
   delete name  
   delete all   # deletes all objects

   name = name of object or selection
  
PYMOL API
 
   cmd.delete (string name = object-or-selection-name )

SEE ALSO

   remove
   '''
   try:
      lock()   
      r = _cmd.delete(str(name))
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

def full_screen(toggle=1):
   '''
DESCRIPTION

   "full_screen" enables or disables PyMOL's full_screen mode.  This
   is only functions well on PC's.

USAGE

   full_screen on
   full_screen off
   
'''
   toggle = toggle_dict[toggle_sc.auto_err(str(toggle),'toggle')]
   if thread.get_ident() == pymol.glutThread:
      try: 
         lock()
         r = _cmd.full_screen(int(toggle))
      finally:
         unlock()
   else:
      try:
         lock()
         r = _cmd.do("_cmd.full_screen(%d)"%int(toggle))
      finally:
         unlock()
   return r

def png(filename):
   '''
DESCRIPTION
  
   "png" writes a png format image file of the current image to disk.
   
USAGE
 
   png filename
 
PYMOL API
 
   cmd.png( string file )
   '''
   if thread.get_ident() ==pymol.glutThread:
      r = _png(str(filename))
   else:
      r = _cmd.do("cmd._png('"+str(filename)+"')")
   return r

def _png(a): # INTERNAL - can only be safely called by GLUT thread 
   try:
      lock()   
      fname = a
      if not re.search("\.png$",fname):
         fname = fname +".png"
      fname = os.path.expanduser(fname)
      fname = os.path.expandvars(fname)         
      r = _cmd.png(str(fname))
   finally:
      unlock()
   return r

def export_coords(obj,state): # experimental
   r = None
   try:
      lock()   
      r = _cmd.export_coords(str(obj),int(state)-1)
   finally:
      unlock()
   return r

def import_coords(obj,state,mechio): # experimental
   r = None
   try:
      lock()   
      r = _cmd.import_coords(str(obj),int(state)-1,mechio)
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

def _special(k,x,y): # INTERNAL (invoked when special key is pressed)
   k=int(k)
   if special.has_key(k):
      if special[k][1]:
         apply(special[k][1],special[k][2],special[k][3])
   return None

def _ctrl(k):
   if ctrl.has_key(k):
      if ctrl[k][0]!=None:
         apply(ctrl[k][0],ctrl[k][1],ctrl[k][2])
   return None


def set_key(key,fn,arg=(),kw={}):  
   '''
DESCRIPTION
  
   "set_key" binds a specific python function to a key press.
   
PYMOL API
 
   cmd.set_key( string key, function fn, tuple arg=(), dict kw={})
 
PYTHON EXAMPLE
 
   from pymol import cmd
 
   def color_blue(object):
      cmd.color("blue",object)
    
   cmd.set_key( 'F1' , make_it_blue, ( "object1" ) )
   cmd.set_key( 'F2' , make_it_blue, ( "object2" ) )
 
   // would turn object1 blue when the F1 key is pressed and
   // would turn object2 blue when the F2 key is pressed.

SEE ALSO

   button
   '''
   r = 0
   for a in special.keys():
      if special[a][0]==key:
         special[a][1]=fn
         special[a][2]=arg
         special[a][3]=kw
         r = 1
   if not r:
      print "Error: special '%s' key not found."%key
   return r
   
def mstop():
   '''
DESCRIPTION
  
   "mstop" stops the movie.
   
USAGE
 
   mstop
 
PYMOL API
 
   cmd.mstop()

SEE ALSO

   mplay, mset, mdo, mclear, mmatrix
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

SEE ALSO

   mstop, mset, mdo, mclear, mmatrix
   '''
   try:
      lock()   
      r = _cmd.mplay(1)
   finally:
      unlock()
   return r

def mray(): # deprecated
   try:
      lock()   
      r = _cmd.mplay(2)
   finally:
      unlock()
   return r

def viewport(width=-1,height=-1):
   '''
DESCRIPTION
  
   "viewport" changes the size of the viewing port (and thus the size
   of all png files subsequently output)
      
USAGE
 
   viewport width, height
 
PYMOL API
  
   cmd.viewport(int width, int height)
   '''
   r = None
   if not is_glut_thread():
      do("viewport %d,%d"%(int(width),int(height)))
   else:
      try:
         lock()
         r = _cmd.viewport(int(width),int(height))
      finally:
         unlock()
   return r

def mdo(frame,command):
   '''
DESCRIPTION
  
   "mdo" sets up a command to be executed upon entry into the
   specified frame of the movie.  These commands are usually created
   by a PyMOL utility program (such as util.mrock).  Command can
   actually contain several commands separated by semicolons ';'
 
USAGE
 
   mdo frame : command
 
PYMOL API
  
   cmd.mdo( int frame, string command )
 
EXAMPLE
 
   // Creates a single frame movie involving a rotation about X and Y
   
   load test.pdb
   mset 1
   mdo 1, turn x,5; turn y,5;
   mplay
   
NOTES
 
   The "mset" command must first be used to define the movie before
   "mdo" statements will have any effect.  Redefinition of the movie
   clears any existing mdo statements.

SEE ALSO

   mset, mplay, mstop
   '''
   try:
      lock()   
      r = _cmd.mdo(int(frame)-1,str(command),0)
   finally:
      unlock()
   return r

def mappend(frame,command):
   '''
DESCRIPTION
  
USAGE
 
   mappend frame : command
 
PYMOL API
  
EXAMPLE
   
NOTES
 
SEE ALSO

   mset, mplay, mstop
   '''
   try:
      lock()   
      r = _cmd.mdo(int(frame)-1,str(";"+command),1)
   finally:
      unlock()
   return r

def dummy(*arg):
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

SEE ALSO

   mset, backward, rewind
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

SEE ALSO

   mset, forward, rewind
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

def test(winid): # generic test routine for development
   try:
      lock()   
      r=_cmd.test(winid)
   finally:
      unlock()
   return r

def dump(fnam,obj):
   try:
      lock()
      r = _cmd.dump(str(fnam),obj)
   finally:
      unlock()
   return r

def multisave(filename,object,state=0):
   r = 1
   try:
      lock()
      _cmd.multisave(str(filename),str(object),int(state)-1,0)
   finally:
      unlock()
   return r

def save(filename,selection='(all)',state=0,format=''):
   '''
DESCRIPTION
  
   "save" writes selected atoms to a file.  The file format is
   autodetected if the extesion is ".pdb" or ".pkl"
 
USAGE
 
   save file [,(selection) [,state [,format]] ]
 
PYMOL API
  
   cmd.save(file, selection, state, format)

SEE ALSO

   load, get_model
   '''
   # preprocess selection
   selection = selector.process(selection)
   #   
   r = 1
   if format=='':
      format = 'pdb'
      if re.search("\.pdb$|\.ent$",filename):
         format = 'pdb'
      elif re.search("\.mol$",filename):
         format = 'mol'
      elif re.search("\.pkl$",filename):
         format = 'pkl'
      elif re.search("\.pkl$",filename):
         format = 'pkla'
      elif re.search("\.mmd$",filename):
         format = 'mmod'
      elif re.search("\.mmod$",filename):
         format = 'mmod'
      elif re.search("\.pmo$",filename):
         format = 'pmo'
   else:
      format = str(format)
   filename = os.path.expanduser(filename)
   filename = os.path.expandvars(filename)
   if format=='pdb':
      f=open(filename,"w")
      if f:
         try:
            lock()
            st = _cmd.get_pdb(str(selection),int(state)-1)
         finally:
            unlock()
            f.write(st)
            f.close()
         r = None
         print " Save: wrote \""+filename+"\"."
   elif format=='pkl': # python binary
      io.pkl.toFile(get_model(selection,state),filename)
      print " Save: wrote \""+filename+"\"."
   elif format=='pkla': # ascii override
      io.pkl.toFile(get_model(selection),filename,bin=0)
      print " Save: wrote \""+filename+"\"."
   elif format=='mmod': # macromodel
      io.mmd.toFile(get_model(selection),filename)
      print " Save: wrote \""+filename+"\"."
   elif format=='mol': 
      io.mol.toFile(get_model(selection),filename)
      print " Save: wrote \""+filename+"\"."
   return r

def get_model(selection="(all)",state=1):
   '''
DESCRIPTION
  
   "get_model" returns a Chempy "Indexed" format model from a selection.
 
PYMOL API
 
   cmd.get_model(string selection [,int state] )
 
   '''
   # preprocess selection
   selection = selector.process(selection)
   #   
   r = 1
   try:
      lock()
      r = _cmd.get_model(str(selection),int(state)-1)
   finally:
      unlock()
   return r

def get_area(selection="(all)",state=1,load_b=0):
   '''
   PRE-RELEASE functionality - API will change
   '''
   # preprocess selection
   selection = selector.process(selection)
   #      
   try:
      lock()
      r = _cmd.get_area(str(selection),int(state)-1,int(load_b))
   finally:
      unlock()
   return r

def get_names(type='objects'):
   '''
DESCRIPTION
  
   "get_names" returns a list of object and/or selection names.
 
PYMOL API
 
   cmd.get_names( [string: "objects"|"selections"|"all"] )
 
NOTES
 
   The default behavior is to return only object names.

SEE ALSO

   get_type, count_atoms, count_states
   '''
   mode = 1
   if type=='objects':
      mode = 1
   elif type=='selections':
      mode = 2
   elif type=='all':
      mode = 0
   elif type=='public':
      mode = 3
   try:
      lock()
      r = _cmd.get_names(int(mode))
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

SEE ALSO

   get_names
   '''
   try:
      lock()
      r = _cmd.get_type(str(name))
   finally:
      unlock()
   return r

def get_state():
   '''
DESCRIPTION
  
   "get_state" returns the current state index (1-based)
 
PYMOL API
 
   cmd.get_state()
 
NOTES
 
   States refer to different geometric configurations which an object
   can above.  By default, states and movie frames have a one-to-one
   relationship.  States can be visited in an arbitrary order to
   create frames.  The "mset" command allows you to build a
   relationship between states and frames.

SEE ALSO

   get_frame
   '''
   # NOTE: NO LOCKS...this is/can be called from cmd.refresh()
   r = _cmd.get_state()+1
   return r

def get_frame():
   '''
DESCRIPTION
  
   "get_frame" returns the current frame index (1-based)
 
PYMOL API
 
   Frames refers to sequences of images in a movie.  Sequential frames
   may contain identical molecular states, they may have one-to-one
   correspondance to molecular states (default), or they may have an
   arbitrary relationship, specific using the "mset" command.

SEE ALSO

   get_state
 
   '''
   # NOTE: NO LOCKS...this is/can be be called from cmd.refresh()
   r = _cmd.get_frame()
   return r


def id_atom(selection):
   '''
DESCRIPTION
  
   "id_atom" returns the original source id of a single atom, or
   raises and exception if the atom does not exist or if the selection
   corresponds to multiple atoms.
 
PYMOL API
 
   list = cmd.id_atom(string selection)
   '''
   r = -1
   selection = string(selection)
   l = apply(identify,(selection,))
   ll = len(l)
   if not ll:
      print "Error: atom %s not found by id_atom." % selection
      raise QuietException
   elif ll>1:
      print "Error: multiple atoms %s found by id_atom." % selection
      raise QuietException
   else:
      r = l[0]
   return r

def identify(selection="(all)"):
   '''
DESCRIPTION
  
   "identify" returns a list of atom IDs corresponding to the ID code
   of atoms in the selection.
 
PYMOL API
 
   list = cmd.identify(string selection="(all)")
 
   '''
   # preprocess selection
   selection = selector.process(selection)
   #      
   r = []
   try:
      lock()
      r = _cmd.identify(str(selection),0) # 0 = default mode
   finally:
      unlock()
   return r

def index(selection="(all)"):
   '''
DESCRIPTION
  
   "index" returns a list of tuples corresponding to the
   object name and index of the atoms in the selection.
 
PYMOL API
 
   list = cmd.index(string selection="(all)")

NOTE

  Atom indices are fragile and will change as atoms are added
  or deleted.  Whenever possible, use integral atom identifiers
  instead of indices.
  
   '''
   # preprocess selection
   selection = selector.process(selection)
   #      
   r = []
   try:
      lock()
      r = _cmd.index(str(selection),0) # 0 = default mode
   finally:
      unlock()
   return r

def find_pairs(selection1,selection2,state1=1,state2=1,cutoff=3.5,mode=0,angle=45):
   '''
DESCRIPTION
  
   "find_pairs" is currently undocumented.
 
   '''
   # preprocess selection
   selection1 = selector.process(selection1)
   selection2 = selector.process(selection2)
   #      
   r = []
   try:
      lock()
      r = _cmd.find_pairs(str(selection1),str(selection2),
                          int(state1)-1,int(state2)-1,int(mode),float(cutoff),float(angle))
      # 0 = default mode
   finally:
      unlock()
   return r

def get_extent(selection="(all)",state=0):
   '''
DESCRIPTION
  
   "get_extent" returns the minimum and maximum XYZ coordinates of a
   selection as an array:
    [ [ min-X , min-Y , min-Z ],[ max-X, max-Y , max-Z ]]
 
PYMOL API
 
   cmd.get_extent(string selection="(all)", state=0 )
 
   '''
   # preprocess selection
   selection = selector.process(selection)
   #      
   r = 1
   try:
      lock()
      r = _cmd.get_min_max(str(selection),int(state)-1)
   finally:
      unlock()
   return r

def create(name,selection,source_state=0,target_state=0):
   '''
DESCRIPTION
  
   "create" creates a new molecule object from a selection.  It can
   also be used to create states in an existing object.
 
   NOTE: this command has not yet been throughly tested.
 
USAGE

   create name, (selection) [,source_state [,target_state ] ]
   
   create name = (selection) [,source_state [,target_state ] ]
     # (DEPRECATED)

   name = object to create (or modify)
   selection = atoms to include in the new object
   source_state (default: 0 - copy all states)
   target_state (default: 0)
   
PYMOL API
  
   cmd.create(string name, string selection, int state, int target_state)

NOTES

   If the source and target states are zero (default), all states will
   be copied.  Otherwise, only the indicated states will be copied.

SEE ALSO

   load, copy
   '''
   # preprocess selection
   selection = selector.process(selection)
   #      
   try:
      lock()
      if name==None:
         sel_cnt = _cmd.get("sel_counter") + 1.0
         _cmd.legacy_set("sel_counter","%1.0f" % sel_cnt)
         name = "obj%02.0f" % sel_cnt
      _cmd.create(str(name),str(selection),
                  int(source_state)-1,int(target_state)-1)
   finally:
      unlock()
   return None

def get_feedback(): # INTERNAL
   l = []
   if lock_attempt():
      try:
         r = _cmd.get_feedback()
         while r:
            l.append(r)
            r = _cmd.get_feedback()
      finally:
         unlock()
#   try:
#      lock()
#      r = _cmd.get_feedback()
#      while r:
#         l.append(r)
#         r = _cmd.get_feedback()
#   finally:
#      unlock()
   return l

def load_coords(*arg): # UNSUPPORTED
   r = 1
   try:
      lock()
      ok = 1
      ftype = loadable.model
      state = 0
      model = arg[0];
      if len(arg)<2:
         ok=0
      if len(arg)>=2:
         oname = string.strip(arg[1])
      if len(arg)>=3:
         state = int(arg[2])-1
      if ok:
         r = _cmd.load_coords(str(oname),model,
                              int(state)-1,int(ftype))
      else:
         print "Error: invalid arguments."
   finally:
      unlock()
   return r

def finish_object(name):
   '''
DESCRIPTION

   "finish_object" is used in cases where many individual states are
   being loaded and it is advantageos to avoid processing them until
   all states have been loaded into RAM.  This function should always
   be called after loading an object with the finish flag set to zero.

PYMOL API

   cmd.finish(string name)

   "name" should be the name of the object
   '''
   r = 1
   try:
      lock()   
      r = _cmd.finish_object(name)
   finally:
      unlock()
   return r

def load_object(type,object,name,state=0,finish=1,discrete=0):
      # assume first argument is the object type (numeric)
   '''
DESCRIPTION

   "load_object" is a general developer function for loading Python objects
   into PyMOL.

PYMOL API

   cmd.load_object(type,object,name,state=0,finish=1,discrete=0)

   NOTE type is one one of the numberic cmd.loadable types
   '''
   r = 1
   try:
      lock()   
      r = _cmd.load_object(str(name),object,int(state)-1,
                              int(type),int(finish),int(discrete))
   finally:
      unlock()
   return r
   
def load_brick(*arg):
   '''
Temporary routine for GAMESS-UK project.
'''
   lst = [loadable.brick]
   lst.extend(list(arg))
   return apply(load_object,lst)

def load_map(*arg):
   '''
Temporary routine for the Phenix project.
'''
   
   lst = [loadable.map]
   lst.extend(list(arg))
   return apply(load_object,lst)

def load_callback(*arg):
   '''
DESCRIPTION

   "load_callback" is used to load a generic Python callback object.
   These objects are called every time the screen is updated and can be used
   to trigger OpenGL rendering calls (such as with PyOpenGL).

PYMOL API

   cmd.load_callback(object,name,state,finish,discrete)
   
'''
   
   lst = [loadable.callback]
   lst.extend(list(arg))
   return apply(load_object,lst)

def load_cgo(*arg):
   '''
DESCRIPTION

   "load_cgo" is used to load a compiled graphics object, which is
   actually a list of floating point numbers built using the constants
   in the $PYMOL_PATH/modules/pymol/cgo.py file.

PYMOL API

   cmd.load_cgo(object,name,state,finish,discrete)
   
'''
   
   lst = [loadable.cgo]
   lst.extend(list(arg))
   return apply(load_object,lst)

def load_model(*arg,**kw):
   '''
DESCRIPTION
  
   "load_model" reads a ChemPy model into an object

PYMOL API
  
   cmd.load_model(model, object [,state [,finish [,discrete ]]])
   '''
   lst = [loadable.model]
   lst.extend(list(arg))
   return apply(load_object,lst,kw)

def _load(oname,finfo,state,ftype,finish,discrete):
   # caller must already hold API lock
   # NOTE: state index assumes 1-based state
   r = 1
   if ftype not in (loadable.model,loadable.brick):
      if ftype == loadable.r3d:
         import cgo
         obj = cgo.from_r3d(finfo)
         if obj:
            _cmd.load_object(str(oname),obj,int(state)-1,loadable.cgo,
                             int(finish),int(discrete))
         else:
            print " load: couldn't load raster3d file."
      elif ftype == loadable.cc1: # ChemDraw 3D
         obj = io.cc1.fromFile(finfo)
         if obj:
            _cmd.load_object(str(oname),obj,int(state)-1,loadable.model,
                             int(finish),int(discrete))            
      else:
         r = _cmd.load(str(oname),finfo,int(state)-1,int(ftype),
                       int(finish),int(discrete))
   else:
      try:
         x = io.pkl.fromFile(finfo)
         if isinstance(x,types.ListType) or isinstance(x,types.TupleType):
            for a in x:
               r = _cmd.load_object(str(oname),a,int(state)-1,
                                    int(ftype),0,int(discrete))
               if(state>0):
                  state = state + 1
            _cmd.finish_object(str(oname))
         else:
            r = _cmd.load_object(str(oname),x,
                                 int(state)-1,int(ftype),
                                 int(finish),int(discrete))
      except:
         print 'Error: can not load file "%s"' % finfo
   return r

def load(filename,object='',state=0,format='',finish=1,discrete=0):
   '''
DESCRIPTION
  
   "load" reads several file formats.  The file extension is used to
   determine the format.  PDB files must end in ".pdb", MOL files must
   end in ".mol", Macromodel files must end in ".mmod", XPLOR
   maps must end in ".xplor", CCP4 maps must end in ".ccp4",
   Raster3D input (Molscript output) must end in ".r3d".

   Pickled ChemPy models with a ".pkl" can also be directly read.
 
   If an object is specified, then the file is load into that object.
   Otherwise, an object is created with the same name as the file
   prefix.
 
USAGE
 
   load filename [,object [,state [,format [,finish [,discrete ]]]]]
 
PYMOL API
  
   cmd.load( filename [,object [,state [,format [,finish [,discrete ]]]]]

NOTES

   You can override the file extension by giving a format string:

   'pdb' : PDB,  'mmod' : Macromodel, 'xyz' : Tinker, 'cc1' : ChemDraw3D  
   'mol' : MDL MOL-file, 'sdf' : MDL SD-file
   'xplor' : X-PLOR/CNS map, 'ccp4' : CCP4 map,
   'callback' : PyMOL Callback object (PyOpenGL)
   'cgo' : compressed graphics object (list of floats)

SEE ALSO

   save
   '''
   r = 1
   try:
      lock()
      type = format
      ftype = 0
      state = int(state)
      finish = int(finish)
      discrete = int(discrete)
      fname = filename
      fname = os.path.expanduser(fname)
      fname = os.path.expandvars(fname)

      if not len(str(type)):
         # determine file type if possible
         if re.search("\.pdb$|\.ent$",filename,re.I):
            ftype = loadable.pdb
         elif re.search("\.mol$",filename,re.I):
            ftype = loadable.mol
         elif re.search("\.mmod$|\.mmd$|\.dat$|\.out$",filename,re.I):
            ftype = loadable.mmod
         elif re.search("\.xplor$",filename,re.I):
            ftype = loadable.xplor
         elif re.search("\.ccp4$",filename,re.I):
            ftype = loadable.ccp4
         elif re.search("\.pkl$",filename,re.I):
            ftype = loadable.model
         elif re.search("\.r3d$",filename,re.I):
            ftype = loadable.r3d
         elif re.search("\.xyz$",filename,re.I):
            ftype = loadable.xyz
         elif re.search("\.cc1$|\.cc2$",filename,re.I):
            ftype = loadable.cc1
         elif re.search("\.xyz_[0-9]*$",filename,re.I):
            ftype = loadable.xyz
         elif re.search("\.sdf$",filename,re.I):
            ftype = loadable.sdf
         elif re.search("\.pmo$",filename,re.I):
            ftype = loadable.pmo
         else:
            ftype = loadable.pdb # default is PDB
      elif is_string(type):
         try:
            ftype = int(type)
         except:
            type = loadable_sc.auto_err(type,'file type')
            if hasattr(loadable,type):
               ftype = getattr(loadable,type)
            else:
               print "Error: unknown type '%s'",type
               raise QuietException
      else:
         ftype = int(type)
         
# get object name
      if len(str(object))==0:
         oname = re.sub(r".*\/|.*\\","",filename) # strip path
         oname = file_ext_re.sub("",oname) # strip extension
         oname = safe_oname_re.sub("_",oname)
         if not len(oname): # safety
            oname = 'obj01'
      else:
         oname = string.strip(object)

# special handling of sdf files
      if ftype == loadable.sdf:
         ftype = loadable.molstr
         sdf = SDF(filename)
         while 1:
            rec = sdf.read()
            if not rec: break
            r = _load(oname,string.join(rec.get('MOL'),''),state,ftype,0,1)
         del sdf
         _cmd.finish_object(str(oname))
         _cmd.do("zoom (%s)"%oname) 
         ftype = -1
         
# standard file handling

      if ftype>=0:
         r = _load(oname,fname,state,ftype,finish,discrete)
   finally:
      unlock()
   return r

def read_molstr(molstr,name,state=0,finish=1,discrete=1):
   '''
DESCRIPTION
  
   "read_molstr" reads an MDL MOL format file as a string
   
PYMOL API ONLY
 
   cmd.read_molstr( string molstr, string name, int state=0,
      int finish=1, int discrete=1 )

NOTES

   "state" is a 1-based state index for the object, or 0 to append.

   "finish" is a flag (0 or 1) which can be set to zero to improve
   performance when loading large numbers of objects, but you must
   call "finish_object" when you are done.

   "discrete" is a flag (0 or 1) which tells PyMOL that there will be
   no overlapping atoms in the file being loaded.  "discrete"
   objects save memory but can not be edited.
   '''
   r = 1
   try:
      lock()
      r = _cmd.load(str(name),str(molstr),int(state)-1,
                    loadable.molstr,int(finish),int(discrete))
   finally:
      unlock()
   return r

def read_mmodstr(*arg):
   '''
DESCRIPTION

   "read_mmodstr" reads a macromodel format structure from a Python
   string.

'''
   r = 1
   try:
      lock()   
      ftype = 6
      if len(arg)==2:
         oname = string.strip(arg[1])
         r = _cmd.load(str(oname),arg[0],-1,int(ftype),1,1)
      elif len(arg)==3:
         oname = string.strip(arg[1])
         r = _cmd.load(str(oname),arg[0],int(arg[2])-1,int(ftype),1,1)
      else:
         print "argument error."
   finally:
      unlock()
   return r

def read_pdbstr(pdb,name,state=0,finish=1,discrete=0):
   '''
DESCRIPTION
  
   "read_pdbstr" in an API-only function which reads a pdb file from a
   Python string.  This feature can be used to load or update
   structures into PyMOL without involving any temporary files.
   
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
   objects save memory but can not be edited.
'''
   r = 1
   try:
      lock()   
      ftype = loadable.pdbstr
      oname = string.strip(str(name))
      r = _cmd.load(str(oname),pdb,int(state)-1,int(ftype),
                       int(finish),int(discrete))
   finally:
      unlock()
   return r
   
def select(name,selection="",quiet=0,show=0):
   '''
DESCRIPTION
  
   "select" creates a named selection from an atom selection.
 
USAGE
 
   select (selection)
   select name, (selection)
   select name = (selection)            # (DEPRECATED)
 
PYMOL API
  
   cmd.select(string name, string selection)
 
EXAMPLES 

   select near , (ll expand 8)
   select near , (ll expand 8)
   select bb, (name ca,n,c,o )

NOTES

   'help selections' for more information about selections.
   '''   
   try:
      lock()
      if selection=="":
         sel_cnt = _cmd.get("sel_counter") + 1.0
         _cmd.legacy_set("sel_counter","%1.0f" % sel_cnt)
         selection = name
         name = "sel%02.0f" % sel_cnt
      else:
         name = name
      # preprocess selection (note: inside TRY)
      selection = selector.process(selection)
      #
      r = _cmd.select(str(name),str(selection),int(quiet))
      if r and show:
         r = _cmd.onoff(str(name),1);
   finally:
      unlock()
   return r

def phi_psi(selection="(byres pk1)"):
   result = get_phipsi(selection)
   if result!=None:
      kees = result.keys()
      kees.sort()
      feedback('push')
      feedback('disable','executive','actions')
      for a in kees:
         iterate("(%s`%d)"%a,"print ' %-9s "+("( %6.1f, %6.1f )"%result[a])+"'%(resn+'-'+resi+':')")
      feedback('pop')
   else:
      print "Error: can't compute phi_psi"
   return result

def count_atoms(selection="(all)",quiet=0):
   '''
DESCRIPTION
  
   "count_atoms" returns a count of atoms in a selection.
 
USAGE
 
   count (selection)
 
PYMOL API
  
   cmd.count(string selection)
 
   '''
   # preprocess selection
   selection = selector.process(selection)
   #      
   try:
      lock()   
      r = _cmd.select("_count_tmp",str(selection),1)
      _cmd.delete("_count_tmp")
   finally:
      unlock()
   if not int(quiet):
      print " count_atoms: %d atoms"%r
   return r
   
def color(color,selection="(all)"):
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
   # preprocess selection
   selection = selector.process(selection)
   #
   try:
      lock()
      r = _cmd.color(str(color),str(selection),0)
   finally:
      unlock()
   return r

def flag(number,selection):
   '''
DESCRIPTION
  
   "flag" sets the indicated flag for atoms in the selection and
    clears the indicated flag for atoms not in the selection.  This
    is primarily useful for passing selection information into
    Chempy models.
   
USAGE

   flag number, selection
   
   flag number = selection     # (DEPRECATED)

PYMOL API
  
   cmd.flag( int number, string selection )
 
EXAMPLES  
 
   flag 0, (name ca)
   flag 1, (resi 45 x; 6)
 
   '''
   # preprocess selection
   selection = selector.process(selection)
   #      
   try:
      lock()   
      r = _cmd.flag(int(number),str(selection))
   finally:
      unlock()
   return r

def set_color(name,rgb):
   '''
DESCRIPTION
  
   "set_color" defines a new color with color indices (0.0-1.0)
   
USAGE

   set_color name, [ red-float, green-float, blue-float ]
   
   set_color name = [ red-float, green-float, blue-float ]
     # (DEPRECATED)
     
PYMOL API
  
   cmd.set_color( string name, float-list rgb )
 
EXAMPLES 
 
   set_color red = [ 1.0, 0.0, 0.0 ]
   '''
   r = 1
   if is_string(rgb):
      rgb = eval(rgb)
   if not (isinstance(rgb,types.ListType) or isinstance(rgb,types.TupleType)):
      print "Error: color specification must be a list such as [ 1.0, 0.0, 0.0 ]"
   elif len(rgb)!=3:
      print "Error: color specification must be a list such as [ 1.0, 0.0, 0.0 ]"
   else:
      try:
         lock()

         if len(rgb)==3:
            r = _cmd.colordef(str(name),float(rgb[0]),float(rgb[1]),float(rgb[2]))
         else:
            print "Error: invalid color."
      finally:
         unlock()
   return r

def rebuild(selection='all',representation='everything'):
   '''
DESCRIPTION

   "rebuild" forces PyMOL to recreate geometric objects in
   case any of them have gone out of sync.

USAGE
   
   rebuild [selection [, representation ]]
   
PYMOL API

   cmd.rebuild(string selection = 'all', string representation = 'everything')

SEE ALSO

   refresh
'''
   r = 1
   selection = selector.process(selection)
   representation = repres_sc.auto_err(representation,'representation')
   repn = repres[representation];
   try:
      lock()
      r = _cmd.rebuild(selection,repn)
   finally:
      unlock()


def set_title(object,state,text):
   '''
DESCRIPTION

   "set_title" attaches a text string to the state of a particular
   object which can be displayed when the state is active.  This is
   useful for display the energies of a set of conformers.

USAGE
   
   set_title object,state,text

PYMOL API

   cmd.set_title(string object,int state,string text)

'''
   r = 1
   try:
      lock()
      r = _cmd.set_title(str(object),int(state)-1,str(text))
   finally:
      unlock()
      
def mpng(prefix):
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
      r = _mpng(prefix)
   else:
      r = _cmd.do("cmd._mpng('"+prefix+"')")
   return r

def _mpng(*arg): # INTERNAL
   try:
      lock()   
      fname = arg[0]
      if re.search("\.png$",fname):
         fname = re.sub("\.png$","",fname)
      fname = os.path.expanduser(fname)
      fname = os.path.expandvars(fname)
      r = _cmd.mpng_(str(fname))
   finally:
      unlock()
   return r

def show(representation="",selection=""):
   '''
DESCRIPTION
  
   "show" turns on atom and bond representations.
 
   The available representations are:
    
      lines     spheres   mesh      ribbon     cartoon
      sticks    dots      surface   labels
      nonbonded nb_spheres 
   
USAGE
 
   show
   show reprentation [,object]
   show reprentation [,(selection)]
   show (selection)
   
PYMOL API
 
   cmd.show( string representation="", string selection="" )
 
EXAMPLES
 
   show lines,(name ca or name c or name n)
   show ribbon
 
NOTES

   "selection" can be an object name
   "show" alone will turn on lines for all bonds.

SEE ALSO

   hide, enable, disable
   '''
   r=1
   try:
      lock()
      if (representation=="") and (selection==""):
         r = _cmd.showhide("(all)",repres['lines'],1); # show lines by default       
      elif (representation!="") and (selection!=""):
         rep = representation
         rep = repres_sc.auto_err(rep,'representation')
         repn = repres[rep];
         # preprocess selection 
         selection = selector.process(selection)
         #   
         r = _cmd.showhide(str(selection),int(repn),1);
      elif representation=='all':
         r = _cmd.showhide("all",repres['lines'],1); # show lines by default 
      elif (representation[0:1]=='(') or (string.find(representation,'/')>=0):
         # preprocess selection
         selection = selector.process(representation)
         #                  
         r = _cmd.showhide(str(selection),repres['lines'],1);
      else: # selection==""
         rep = representation
         rep = repres_sc.auto_err(rep,'representation')
         repn = repres[rep];
         r = _cmd.showhide("all",int(repn),1);
   finally:
      unlock()
   return r

def hide(representation="",selection=""):
   '''
DESCRIPTION
  
   "hide" turns of atom and bond representations.
 
   The available representations are:
    
      lines     spheres   mesh      ribbon     cartoon
      sticks    dots      surface   labels
      nonbonded nb_spheres
   
USAGE
 
   hide reprentation [,object]
   hide reprentation [,(selection)]
   hide (selection)
   
PYMOL API
 
   cmd.hide( string representation="", string selection="")
 
EXAMPLES
 
   hide lines,all
   hide ribbon

SEE ALSO

   show, enable, disable
   '''
   r = 1
   try:
      lock()
      if (representation=="") and (selection==""):
         r = _cmd.showhide("@",0,0);      
      elif (representation!="") and (selection!=""):
         rep = representation
         rep = repres_sc.auto_err(rep,'representation')
         repn = repres[rep];
         selection = selector.process(selection)
         r = _cmd.showhide(str(selection),int(repn),0);
      elif (representation=='all'):
         r = _cmd.showhide("@",0,0);
      elif (representation[0:1]=='(') or (string.find(representation,'/')>=0):
         selection = selector.process(representation)
         r = _cmd.showhide(str(selection),-1,0);
      else: # selection == ""
         rep = representation
         rep = repres_sc.auto_err(rep,'representation')
         repn = repres[rep];
         r = _cmd.showhide("all",int(repn),0);
   finally:
      unlock()
   return r

def paste(): # INTERNAL
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

   
def button(button,modifier,action):
   '''
DESCRIPTION
  
   "button" can be used to redefine what the mouse buttons do.
   
USAGE
 
   button <button>,<modifier>,<action>
 
PYMOL API
 
   cmd.button( string button, string modifier, string action )
 
NOTES

   button:      L, M, R
   modifers:    None, Shft, Ctrl, CtSh
   actions:     Rota, Move, MovZ, Clip, RotZ, ClpN, ClpF
                lb,   mb,   rb,   +lb,  +lbX, -lbX, +mb,  +rb, 
                PkAt, PkBd, RotF, TorF, MovF, Orig

   Switching from visualization to editing mode will redefine the
   buttons, so do not use the built-in switch if you want to preserve
   your custom configuration.

'''
   r=1
   try:
      lock()
      button = string.lower(button)
      button = button_sc.auto_err(button,'button')
      modifier = string.lower(modifier)
      modifier = but_mod_sc.auto_err(modifier,'modifier')
      action = string.lower(action)
      action = but_act_sc.auto_err(action,'action')
      but_code = button_code[button] + 3*but_mod_code[modifier]
      act_code = but_act_code[action]
      r = _cmd.button(but_code,act_code)
   finally:
      unlock()
   return r
   
def mmatrix(action):
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
      if action=="clear":
         r = _cmd.mmatrix(0)
      elif action=="store":
         r = _cmd.mmatrix(1)
      elif action=="recall":
         r = _cmd.mmatrix(2)
   finally:
      unlock()
   return r

def enable(name='all'):
   '''
DESCRIPTION
  
   "enable" enable display of an object and all currently visible representations.
   
USAGE
 
   enable name
   enable all

   name = object or selection name
   
PYMOL API
 
   cmd.enable( string object-name )
 
EXAMPLE
 
   enable my_object

SEE ALSO

   show, hide, disable
   '''
   try:
      lock()   
      r = _cmd.onoff(str(name),1);
   finally:
      unlock()
   return r

def disable(name='all'):
   '''
DESCRIPTION
  
   "disable" disables display of an object and all currently visible
   representations.
   
USAGE
 
   disable name
   disable all 

   "name" is the name of an object or a named selection
   
PYMOL API
 
   cmd.disable( string name ) 
 
EXAMPLE
 
   disable my_object

SEE ALSO

   show, hide, enable   
   '''
   try:
      lock()   
      r = _cmd.onoff(str(name),0);
   finally:
      unlock()
   return r

def check(selection=None,preserve=0):
# UNSUPPORTED
# This function relies on code that is not currently part of PyMOL/ChemPy
   # NOTE: the realtime module relies on code that is not yet part of PyMOL/ChemPy
   from chempy.tinker import realtime
   if selection==None:
      arg = get_names("objects")
      arg = arg[0:1]
      if arg:
         if len(arg):
            selection = arg
   if selection!=None:
      selection = selector.process(selection)
      realtime.assign("("+selection+")",int(preserve))
      realtime.setup("("+selection+")")

def fast_minimize(*arg):
# OBSOLETE, TO BE REMOVED
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
# OBSOLETE, TO BE REMOVED
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

   "cd" changes the current working directory.

USAGE
   
   cd <path>

SEE ALSO

   pwd, ls, system
   '''
   dir = os.path.expanduser(dir)
   dir = os.path.expandvars(dir)
   os.chdir(dir)


def pwd():
   '''
DESCRIPTION

   Print current working directory.

USAGE
   
   pwd

SEE ALSO

   cd, ls, system
   '''
   print os.getcwd()
   

def ls(pattern=None):
   '''
DESCRIPTION

   List contents of the current working directory.

USAGE
   
   ls [pattern]
   dir [pattern]

EXAMPLES

   ls
   ls *.pml

SEE ALSO

   cd, pwd, system   
   '''
   if pattern==None:
      pattern = "*"
   else:
      pattern = os.path.expanduser(pattern)
      pattern = os.path.expandvars(pattern)
   if string.find("*",pattern)<0:
      lst = glob(pattern+"/*")
   else:
      lst = []
   if not len(lst):
      lst = glob(pattern)
   lst = parsing.list_to_str_list(lst)
   for a in lst:
      print a
      
def mset(specification=""):
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

SEE ALSO

   mdo, mplay, mclear
   '''
   try:
      lock()
      output=[]
      input = string.split(string.strip(specification))
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

def extend(name,function):
   '''
DESCRIPTION

   "extend" is an API-only function which binds a new external
   function as a command into the PyMOL scripting language.

PYMOL API

   cmd.extend(string name,function function)
   
PYTHON EXAMPLE

   def foo(moo=2): print moo
   cmd.extend('foo',foo)

   The following would now be valid within PyMOL:

   foo
   foo 3
   foo moo=5

SEE ALSO

   alias, api
   '''
   keyword[name] = [function, 0,0,',',parsing.STRICT]
   kwhash.append(name)
   
class fb_action:
   set = 0
   enable = 1
   disable = 2
   push = 3
   pop = 4
   
fb_action_sc = Shortcut(fb_action.__dict__.keys())

class fb_module:

# This first set represents internal C systems

   all                       =0
   isosurf                   =1
   map                       =2
   matrix                    =3
   mypng                     =4
   triangle                  =5
   match                     =6
   raw                       =7
   
   feedback                  =12
   scene                     =13
   threads                   =14  
   symmetry                  =15
   ray                       =16
   setting                   =17
   object                    =18
   ortho                     =19
   movie                     =20
   python                    =21
   extrude                   =22
   rep                       =23
   
   coordset                  =25
   distset                   =26

   objectmolecule            =30
   objectmap                 =31
   objectmesh                =32
   objectdist                =33 
   objectcgo                 =34
   objectcallback            =35

   repwirebond               =45
   repcylbond                =46
   replabel                  =47
   repsphere                 =49
   repsurface                =50
   repmesh                   =51
   repdot                    =52
   repnonbonded              =53
   repnonbondedsphere        =54
   repdistdash               =55
   repdistlabel              =56
   repribbon                 =57
   repcartoon                =58
   
   executive                 =70
   selector                  =71
   editor                    =72

   export                    =75
   ccmd                      =76
   api                       =77   
   
   main                      =80  

# This second set, with negative indices
# represent python level systems

   parser                    =-1
   cmd                       =-2
   
fb_module_sc = Shortcut(fb_module.__dict__.keys())

class fb_mask:
   results =             0x01
   errors =              0x02
   actions =             0x04
   warnings =            0x08
   details =             0x10
   blather =             0x20
   debugging =           0x80
   everything =          0xFF
   
fb_mask_sc = Shortcut(fb_mask.__dict__.keys())

fb_dict ={}

for a in fb_module.__dict__.keys():
   vl = getattr(fb_module,a)
   if vl<0:
      fb_dict[vl] = 0x1F # default mask

fb_debug = sys.stderr # can redirect python debugging output elsewhere if desred...

def _feedback(module,mask): # feedback test routine
   r = 0
   module = int(module)
   mask = int(mask)
   if module>0:
      try:
         lock()
         r = _cmd.feedback(module,mask)
      finally:
         unlock()
   else:
      if fb_dict.has_key(module):
         r = fb_dict[module]&mask
   return r
   
def feedback(action="?",module="?",mask="?"):
   '''
DESCRIPTION

   "feedback" allows you to control what and how much text is output
   from PyMOL.

USAGE
   
   feedback action,module,mask

   action is one of ['set','enable','disable']
   module is a space-separated list of strings or simply "all"
   mask is a space-separated list of strings or simply "everything"

NOTES:

   "feedback" alone will print a list of the available choices
   
PYMOL API

   cmd.feedback(string action,string module,string mask)
   
EXAMPLES

   feedback enable, all , debugging
   feedback disable, selector, warnings actions
   feedback enable, main, blather
'''
   r = None

   # validate action

   if action=="?":
      print " feedback: available actions: set, enable, disable"
      act_kee = 0
   else:
      act_kee = fb_action_sc.interpret(action)
      if act_kee == None:
         print "Error: invalid feedback action '%s'."%action
         raise QuietException
      elif not is_string(act_kee):
         print "Error: ambiguous feedback action '%s'."%action
         print action_amb
         raise QuietException
      act_int = int(getattr(fb_action,act_kee))

   if (act_int<3) and ("?" in [action,module,mask]):
      if module=="?":
         print " feedback: available modules:"
         for a in fb_module.__dict__.keys():
            if a[0]!='_':
               print "   ",a
      if mask=="?":
         print " feedback: available masks:"
         for a in fb_mask.__dict__.keys():
            if a[0]!='_':
               print "   ",a
   else:
      if (act_int>=3):
         module='all'
         mask='everything'
            
      # validate and combine masks
      
      mask_int = 0
      mask_lst = string.split(mask)
      for mask in mask_lst:
         mask_kee = fb_mask_sc.interpret(mask)
         if mask_kee == None:
            print "Error: invalid feedback mask '%s'."%mask
            raise QuietException
         elif not is_string(mask_kee):
            print "Error: ambiguous feedback mask '%s'."%mask
            raise QuietException
         mask_int = int(getattr(fb_mask,mask_kee))

      # validate and iterate modules
      
      mod_lst = string.split(module)
      for module in mod_lst:
         mod_kee = fb_module_sc.interpret(module)
         if mod_kee == None:
            print "Error: invalid feedback module '%s'."%module
            raise QuietException         
         elif not is_string(mod_kee):
            print "Error: ambiguous feedback module '%s'."%module
            raise QuietException         
         mod_int = int(getattr(fb_module,mod_kee))
         if mod_int>=0:
            try:
               lock()
               r = _cmd.set_feedback(act_int,mod_int,mask_int)
            finally:
               unlock()
         if mod_int<=0:
            if mod_int:
               if act_int==0:
                  fb_dict[mod_int] = mask_int
               elif act_int==1:
                  fb_dict[mod_int] = fb_dict[mod_int] | mask_int
               elif act_int==2:
                  fb_dict[mod_int] = fb_dict[mod_int] & ( 0xFF - mask_int )
            else:
               for mod_int in fb_dict.keys():
                  if act_int==0:
                     fb_dict[mod_int] = mask_int
                  elif act_int==1:
                     fb_dict[mod_int] = fb_dict[mod_int] | mask_int
                  elif act_int==2:
                     fb_dict[mod_int] = fb_dict[mod_int] & ( 0xFF - mask_int )
            if _feedback(fb_module.feedback,fb_mask.debugging):
                sys.stderr.write(" feedback: mode %d on %d mask %d\n"%(
                   act_int,mod_int,mask_int))
   return r

def alias(name,command):
   '''
DESCRIPTION

   "alias" allows you to bind a commonly used command to a single word

USAGE
   
   alias name, command-sequence

PYMOL API

   cmd.alias(string name,string command)
   
EXAMPLES

   alias go,load "test.pdb"; zoom (i;500); show sticks,(i;500 a;4)
   go

SEE ALSO

   extend
   '''
   keyword[name] = [eval("lambda :do('''%s ''')"%command), 0,0,',',parsing.STRICT]
   kwhash.append(name)
   
#####################################################################
import util
import movie

keyword = {

   # keyword : [ command, # min_arg, max_arg, separator, mode ]

   # NOTE: min_arg, max_arg, and separator, are hold-overs from the
   #       original PyMOL parser which will eventually be removed.
   #       all new commands should use NO_CHECK or STRICT modes
   #       which make much better use of built-in python features.
   
   'abort'         : [dummy        , 0 , 0 , ',' , parsing.ABORT  ],
   'alias'         : [alias        , 0 , 0 , ''  , parsing.LITERAL1 ],
   'align'         : [align        , 0 , 0 , ''  , parsing.STRICT ],
   'alter'         : [alter        , 0 , 0 , ''  , parsing.LITERAL1 ],
   'alter_state'   : [alter_state  , 0 , 0 , ''  , parsing.LITERAL2 ],
   'api'           : [api          , 0 , 0 , ''  , parsing.STRICT ],
   'backward'      : [backward     , 0 , 0 , ''  , parsing.STRICT ],
   'bg_color'      : [bg_color     , 0 , 0 , ''  , parsing.STRICT ],
   'bond'          : [bond         , 0 , 0 , ''  , parsing.STRICT ],
   'button'        : [button       , 0 , 0 , ''  , parsing.STRICT ],
   'cartoon'       : [cartoon      , 0 , 0 , ''  , parsing.STRICT ],
   'cd'            : [cd           , 0 , 0 , ''  , parsing.STRICT ],  
   'check'         : [check        , 0 , 0 , ''  , parsing.STRICT ],
   'clip'          : [clip         , 0 , 0 , ''  , parsing.STRICT ],
   'cls'           : [cls          , 0 , 0 , ''  , parsing.STRICT ],
   'color'         : [color        , 0 , 0 , ''  , parsing.STRICT ],
   'commands'      : [commands     , 0 , 0 , ''  , parsing.STRICT ],
   'copy'          : [copy         , 0 , 0 , ''  , parsing.LEGACY ],
   'count_atoms'   : [count_atoms  , 0 , 0 , ''  , parsing.STRICT ],
   'count_states'  : [count_states , 0 , 0 , ''  , parsing.STRICT ],
   'cycle_valence' : [cycle_valence, 0 , 0 , ''  , parsing.STRICT ],
   'create'        : [create       , 0 , 0 , ''  , parsing.LEGACY ],   
   'delete'        : [delete       , 0 , 0 , ''  , parsing.STRICT ],
   'deprotect'     : [deprotect    , 0 , 0 , ''  , parsing.STRICT ],
   'dir'           : [ls           , 0 , 0 , ''  , parsing.STRICT ],  
   'disable'       : [disable      , 0 , 0 , ''  , parsing.STRICT ],
   'distance'      : [distance     , 0 , 0 , ''  , parsing.LEGACY ],   
   'dump'          : [dump         , 0 , 0 , ''  , parsing.STRICT ],
   'edit'          : [edit         , 0 , 0 , ''  , parsing.STRICT ],
   'edit_mode'     : [edit_mode    , 0 , 0 , ''  , parsing.STRICT ],
   'enable'        : [enable       , 0 , 0 , ''  , parsing.STRICT ],
   'ending'        : [ending       , 0 , 0 , ''  , parsing.STRICT ],
   'export_dots'   : [export_dots  , 0 , 0 , ''  , parsing.STRICT ],
   'extend'        : [extend       , 0 , 0 , ''  , parsing.STRICT ],
   'fast_minimize' : [fast_minimize, 1,  4 , ',' , parsing.SIMPLE ], # TO REMOVE
   'feedback'      : [feedback     , 0,  0 , ''  , parsing.STRICT ],
   'fit'           : [fit          , 0 , 0 , ''  , parsing.STRICT ],
   'flag'          : [flag         , 0 , 0 , ''  , parsing.LEGACY ],
   'fork'          : [spawn        , 1 , 2 , ',' , parsing.SPAWN  ],
   'forward'       : [forward      , 0 , 0 , ''  , parsing.STRICT ],
   'fragment'      : [fragment     , 0 , 0 , ''  , parsing.STRICT ],
   'full_screen'   : [full_screen  , 0 , 0 , ''  , parsing.STRICT ],
   'fuse'          : [fuse         , 0 , 0 , ''  , parsing.STRICT ],
   'frame'         : [frame        , 0 , 0 , ''  , parsing.STRICT ],
   'get_view'      : [get_view     , 0 , 0 , ''  , parsing.STRICT ],
   'h_add'         : [h_add        , 0 , 0 , ''  , parsing.STRICT ],
   'h_fill'        : [h_fill       , 0 , 0 , ''  , parsing.STRICT ],
   'help'          : [help         , 0 , 0 , ''  , parsing.STRICT ],
   'hide'          : [hide         , 0 , 0 , ''  , parsing.STRICT ],
   'intra_fit'     : [intra_fit    , 0 , 0 , ''  , parsing.STRICT ],
   'intra_rms'     : [intra_rms    , 0 , 0 , ''  , parsing.STRICT ],
   'intra_rms_cur' : [intra_rms_cur, 0 , 0 , ''  , parsing.STRICT ],
   'invert'        : [invert       , 0 , 0 , ''  , parsing.STRICT ],
   'isodot'        : [isodot       , 0 , 0 , ''  , parsing.LEGACY ],   
   'isomesh'       : [isomesh      , 0 , 0 , ''  , parsing.LEGACY ],
   'iterate'       : [iterate      , 0 , 0 , ''  , parsing.LITERAL1 ],
   'iterate_state' : [iterate_state, 0 , 0 , ''  , parsing.LITERAL2 ],
   'label'         : [label        , 0 , 0 , ''  , parsing.LITERAL1 ],
   'load'          : [load         , 0 , 0 , ''  , parsing.STRICT ],
   'log'           : [log          , 0 , 0 , ''  , parsing.STRICT ],
   'log_close'     : [log_close    , 0 , 0 , ''  , parsing.STRICT ],
   'log_open'      : [log_open     , 0 , 0 , ''  , parsing.STRICT ],
   'ls'            : [ls           , 0 , 0 , ''  , parsing.STRICT ],  
   'mask'          : [mask         , 0 , 0 , ''  , parsing.STRICT ],
   'map_set_border': [map_set_border,0 , 0 , ''  , parsing.STRICT ],    
   'mappend'       : [mappend      , 2 , 2 , ':' , parsing.SINGLE ], 
   'mem'           : [mem          , 0 , 0 , ''  , parsing.STRICT ],
   'meter_reset'   : [meter_reset  , 0 , 0 , ''  , parsing.STRICT ],
   'move'          : [move         , 0 , 0 , ''  , parsing.STRICT ],
   'mset'          : [mset         , 0 , 0 , ''  , parsing.STRICT ],
   'mdo'           : [mdo          , 2 , 2 , ':' , parsing.SINGLE ],
   'mpng'          : [mpng         , 0 , 0 , ''  , parsing.STRICT ],
   'mplay'         : [mplay        , 0 , 0 , ''  , parsing.STRICT ],
   'mray'          : [mray         , 0 , 0 , ''  , parsing.STRICT ],
   'mstop'         : [mstop        , 0 , 0 , ''  , parsing.STRICT ],
   'mclear'        : [mclear       , 0 , 0 , ''  , parsing.STRICT ],
   'middle'        : [middle       , 0 , 0 , ''  , parsing.STRICT ],
   'minimize'      : [minimize     , 0 , 4 , ',' , parsing.SIMPLE ], # TO REMOVE
   'mmatrix'       : [mmatrix      , 0 , 0 , ''  , parsing.STRICT ],
   'multisave'     : [multisave    , 0 , 0 , ''  , parsing.STRICT ],   
   'origin'        : [origin       , 0 , 0 , ''  , parsing.STRICT ],
   'orient'        : [orient       , 0 , 0 , ''  , parsing.STRICT ],
   'overlap'       : [overlap      , 0 , 0 , ''  , parsing.STRICT ],
   'pair_fit'      : [pair_fit     , 2 ,98 , ',' , parsing.SIMPLE ],
   'phi_psi'       : [phi_psi      , 0 , 0 , ''  , parsing.STRICT ],
   'protect'       : [protect      , 0 , 0 , ''  , parsing.STRICT ],
   'pwd'           : [pwd          , 0 , 0 , ''  , parsing.STRICT ],
   'ray'           : [ray          , 0 , 0 , ''  , parsing.STRICT ],
   'rebuild'       : [rebuild      , 0 , 0 , ''  , parsing.STRICT ],
   'redo'          : [redo         , 0 , 0 , ''  , parsing.STRICT ],
   'refresh'       : [refresh      , 0 , 0 , ''  , parsing.STRICT ],
   'remove'        : [remove       , 0 , 0 , ''  , parsing.STRICT ],
   'remove_picked' : [remove_picked, 0 , 0 , ''  , parsing.STRICT ],
   'rename'        : [rename       , 0 , 0 , ''  , parsing.STRICT ],
   'replace'       : [replace      , 0 , 0 , ''  , parsing.STRICT ],
   'reset'         : [reset        , 0 , 0 , ''  , parsing.STRICT ],
   'resume'        : [resume       , 0 , 0 , ''  , parsing.STRICT ],
   'rewind'        : [rewind       , 0 , 0 , ''  , parsing.STRICT ],
   'rock'          : [rock         , 0 , 0 , ''  , parsing.STRICT ],
   'run'           : [run          , 1 , 2 , ',' , parsing.RUN    ],
   'rms'           : [rms          , 0 , 0 , ''  , parsing.STRICT ],
   'rms_cur'       : [rms_cur      , 0 , 0 , ''  , parsing.STRICT ],
   'save'          : [save         , 0 , 0 , ''  , parsing.STRICT ],
   'select'        : [select       , 0 , 0 , ''  , parsing.LEGACY ],
   'set'           : [set          , 0 , 0 , ''  , parsing.LEGACY ],
   'set_color'     : [set_color    , 0 , 0 , ''  , parsing.LEGACY ],
   'set_title'     : [set_title    , 0 , 0 , ''  , parsing.STRICT ],   
   'set_key'       : [set_key      , 0 , 0 , ''  , parsing.STRICT ], # API only
   'set_view'      : [set_view     , 0 , 0 , ''  , parsing.STRICT ],   
   'show'          : [show         , 0 , 0 , ''  , parsing.STRICT ],
   'sort'          : [sort         , 0 , 0 , ''  , parsing.STRICT ],
   'spawn'         : [spawn        , 1 , 2 , ',' , parsing.SPAWN  ],
   'spheroid'      : [spheroid     , 0 , 0 , ''  , parsing.STRICT ],
   'splash'        : [splash       , 0 , 0 , ''  , parsing.STRICT ],
   '_special'      : [_special     , 0 , 0 , ''  , parsing.STRICT ],
   'stereo'        : [stereo       , 0 , 0 , ''  , parsing.STRICT ],
   'symexp'        : [symexp       , 0 , 0 , ''  , parsing.LEGACY ],
   'system'        : [system       , 0 , 0 , ''  , parsing.STRICT ],
   'test'          : [test         , 0 , 0 , ''  , parsing.STRICT ],
   'torsion'       : [torsion      , 0 , 0 , ''  , parsing.STRICT ],
   'turn'          : [turn         , 0 , 0 , ''  , parsing.STRICT ],
   'quit'          : [quit         , 0 , 0 , ''  , parsing.STRICT ],
   '_quit'         : [_quit        , 0 , 0 , ''  , parsing.STRICT ],
   'png'           : [png          , 0 , 0 , ''  , parsing.STRICT ],
   'unbond'        : [unbond       , 0 , 0 , ''  , parsing.STRICT ],
   'unpick'        : [unpick       , 0 , 0 , ''  , parsing.STRICT ],
   'undo'          : [undo         , 0 , 0 , ''  , parsing.STRICT ],
   'unmask'        : [unmask       , 0 , 0 , ''  , parsing.STRICT ],
   'unprotect'     : [deprotect    , 0 , 0 , ''  , parsing.STRICT ],
# utility programs 
   'util.cbag'     : [util.cbag    , 0 , 0 , ''  , parsing.STRICT ],
   'util.cbac'     : [util.cbac    , 0 , 0 , ''  , parsing.STRICT ],
   'util.cbay'     : [util.cbay    , 0 , 0 , ''  , parsing.STRICT ],
   'util.cbas'     : [util.cbas    , 0 , 0 , ''  , parsing.STRICT ],
   'util.cbap'     : [util.cbap    , 0 , 0 , ''  , parsing.STRICT ],
   'util.cbak'     : [util.cbak    , 0 , 0 , ''  , parsing.STRICT ],
   'util.cbaw'     : [util.cbaw    , 0 , 0 , ''  , parsing.STRICT ],
   'util.cbab'     : [util.cbab    , 0 , 0 , ''  , parsing.STRICT ],
   'util.cbc'      : [util.cbc     , 0 , 0 , ''  , parsing.STRICT ],
   'util.mrock'    : [util.mrock   , 0 , 0 , ''  , parsing.STRICT ], # LEGACY
   'util.mroll'    : [util.mroll   , 0 , 0 , ''  , parsing.STRICT ], # LEGACY
   'util.ss'       : [util.ss      , 0 , 0 , ''  , parsing.STRICT ],# secondary structure
   'util.rainbow'  : [util.rainbow , 0 , 0 , ''  , parsing.STRICT ],# secondary structure
# movie programs
   'movie.rock'    : [movie.rock   , 0 , 0 , ''  , parsing.STRICT ],
   'movie.roll'    : [movie.roll   , 0 , 0 , ''  , parsing.STRICT ],
   'movie.load'    : [movie.load   , 0 , 0 , ''  , parsing.STRICT ],
   'movie.zoom'    : [movie.zoom   , 0 , 0 , ''  , parsing.STRICT ],
   'movie.screw'   : [movie.screw  , 0 , 0 , ''  , parsing.STRICT ],
#
   'update'        : [update       , 0 , 0 , ''  , parsing.STRICT ],
   'view'          : [view         , 0 , 0 , ''  , parsing.STRICT ],   
   'viewport'      : [viewport     , 0 , 0 , ''  , parsing.STRICT ],
   'wizard'        : [wizard       , 0 , 0 , ''  , parsing.STRICT ],
   'zoom'          : [zoom         , 0 , 0 , ''  , parsing.STRICT ],
   }

kwhash = Shortcut(keyword.keys())

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
   'faster'        : [faster       , 0 , 0 , ',' , 0 ],  
   '@'             : [at_sign      , 0 , 0 , ',' , 0 ],  
}

help_sc = Shortcut(keyword.keys()+help_only.keys())

repres = {
   'everything'    : -1,
   'sticks'        : 0,
   'spheres'       : 1,
   'surface'       : 2,
   'labels'        : 3,
   'nb_spheres'    : 4,
   'cartoon'       : 5,
   'ribbon'        : 6,
   'lines'         : 7,
   'mesh'          : 8,
   'dots'          : 9,
   'dashes'        :10,
   'nonbonded'     :11,
   'cell'          :12,
   'cgo'           :13,
   'callback'      :14,
   'extent'        :15,   
}

repres_sc = Shortcut(repres.keys())

button_code = {
   'left' : 0,
   'middle' : 1,
   'right' : 2,
   }

button_sc = Shortcut(button_code.keys())

but_mod_code = {
   'none'  : 0,
   'shft'  : 1,
   'ctrl'  : 2,
   'ctsh'  : 3
   }

but_mod_sc = Shortcut(but_mod_code.keys())

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
   '+lbx' : 19 ,
   '-lbx' : 20 ,
   'lbbx' : 21 ,
   }

but_act_sc = Shortcut(but_act_code.keys())

special = {
   1        : [ 'F1'        , None                   , () , {} ],
   2        : [ 'F2'        , None                   , () , {} ],
   3        : [ 'F3'        , None                   , () , {} ],
   4        : [ 'F4'        , None                   , () , {} ],
   5        : [ 'F5'        , None                   , () , {} ],
   6        : [ 'F6'        , None                   , () , {} ],
   7        : [ 'F7'        , None                   , () , {} ],
   8        : [ 'F8'        , None                   , () , {} ],
   9        : [ 'F9'        , None                   , () , {} ],
   10       : [ 'F10'       , None                   , () , {} ],
   11       : [ 'F11'       , redo                   , () , {} ],
   12       : [ 'F12'       , undo                   , () , {} ],
   100      : [ 'left'      , backward               , () , {} ],
   101      : [ 'up'        , None                   , () , {} ],
   102      : [ 'right'     , forward                , () , {} ],
   103      : [ 'down'      , None                   , () , {} ],
   104      : [ 'pgup'      , None                   , () , {} ],
   105      : [ 'pgdown'    , None                   , () , {} ],
   106      : [ 'home'      , rewind                 , () , {} ],
   107      : [ 'end'       , ending                 , () , {} ],
   108      : [ 'insert'    , rock                   , () , {} ]   
}

ctrl = {
   'A' : [ redo                   , () , {}],
   'B' : [ replace                , ('Br',1,1), {} ],
   'C' : [ replace                , ('C',4,4), {} ],
   'D' : [ remove_picked          , () , {} ],
   'E' : [ invert                 , () , {} ],   
   'F' : [ replace                , ('F',1,1), {} ],   
   'G' : [ replace                , ('H',1,1), {} ],
   'I' : [ replace                , ('I',1,1), {} ],
   'J' : [ alter                  , ('pk1','formal_charge=-1.0'), {} ],
   'K' : [ alter                  , ('pk1','formal_charge =1.0'), {} ],
   'L' : [ replace                , ('Cl',1,1) , {}],   
   'N' : [ replace                , ('N',4,3) , {}],
   'O' : [ replace                , ('O',4,2) , {}],   
   'P' : [ replace                , ('P',4,1) , {}],
   'Q' : [ h_add                  , ("pk1",) , {}],   
   'R' : [ h_fill                 , () , {} ],   
   'S' : [ replace                , ('S',4,2) , {}],
   'T' : [ bond                   , () , {} ],   
   'U' : [ alter                  , ('pk1','formal_charge =0.0') , {}],
   'W' : [ cycle_valence          , () , {}],   
   'X' : [ None                   , () , {} ],
   'Y' : [ attach                 , ('H',1,1) , {} ],
   'Z' : [ undo                   , () , {} ],   
   }

def get_names_of_type(type):
   obj = get_names('objects')
   types = map(get_type,obj)
   mix = map(None,obj,types)
   lst = []
   for a in mix:
      if a[1]==type:
         lst.append(a[0])
   return lst

selection_sc = lambda sc=Shortcut,gn=get_names:sc(gn('public')+['all'])
object_sc = lambda sc=Shortcut,gn=get_names:sc(gn('objects'))
map_sc = lambda sc=Shortcut,gnot=get_names_of_type:sc(gnot('object:map'))

# Table for argument autocompletion

auto_arg =[
   {
   'set' : [ setting.setting_sc, 'settings', '=' ],
   'show' : [ repres_sc , 'representations',', ' ],
   'hide' : [ repres_sc , 'representations',', ' ],
   'stereo' : [ toggle_sc , 'options','' ],
   'full_screen' : [ toggle_sc , 'options','' ],
   'clip' : [ clip_action_sc , 'clipping actions',', ' ],
   'feedback' : [ fb_action_sc , 'actions',', ' ],
   'button' : [ button_sc , 'buttons',', ' ],
   'zoom' : [ selection_sc , 'selections','' ],
   'origin' : [ selection_sc , 'selections','' ],
   'protect' : [ selection_sc , 'selections','' ],
   'deprotect' : [ selection_sc , 'selections','' ],   
   'mask' : [ selection_sc , 'selections','' ],
   'unmask' : [ selection_sc , 'selections','' ],
   'delete' : [ selection_sc , 'selections','' ],
   'alter' : [ selection_sc , 'selections','' ],
   'iterate' : [ selection_sc , 'selections','' ],
   'iterate_state' : [ selection_sc , 'selections','' ],
   'help' : [ help_sc , 'selections','' ],         
   },
   {
   'feedback' : [ fb_module_sc , 'modules',', ' ],
   'button' : [ but_mod_sc , 'modifiers',', ' ],
   'show' : [ selection_sc , 'selections','' ],
   'hide' : [ selection_sc , 'selections','' ],
   'color' : [ selection_sc , 'selections','' ],
   'select' : [ selection_sc , 'selections','' ],
   'save' : [ selection_sc , 'selections',', ' ],
   'load' : [ selection_sc , 'selections',', ' ],
   'create' : [ selection_sc , 'selections',', ' ],
   'symexp' : [ object_sc , 'objects',', ' ],   
   'isomesh' : [ map_sc , 'map objects',', ' ],            
   },
   {
   'feedback' : [ fb_mask_sc , 'mask','' ],
   'button' : [ but_act_sc , 'button actions','' ],
   }
   ]
   
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
   xyz = 15      # xyz, tinker format
   sdf = 16      # sdf, only used within cmd.py
   cc1 = 17      # cc1 and cc2, only used within cmd.py
   ccp4 = 18     # CCP4 map, under development
   pmo = 19      # pmo, experimental molecular object format
   
loadable_sc = Shortcut(loadable.__dict__.keys()) 

