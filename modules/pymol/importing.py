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

import re
import string
import os
from pymol import cmd
from cmd import _cmd,lock,unlock,Shortcut,QuietException, \
     _feedback,fb_module,fb_mask, \
     file_ext_re,safe_oname_re, \
     _load

from chempy.sdf import SDF,SDFRec

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
      elif cmd.is_string(type):
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
