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

import cmd
from cmd import _cmd,lock,unlock
import selector

def align(source,target): # EXPERIMENTAL, BUGGY
   r = None
   try:
      lock()
      r = _cmd.align(str(source),str(target))
   finally:
      unlock()
   return r

def intra_fit(selection,state=1):
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
      r = _cmd.intrafit("("+str(selection)+")",int(state)-1,2)
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
      r = _cmd.intrafit("("+str(selection)+")",int(state)-1,1)
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
      r = _cmd.intrafit("("+str(selection)+")",int(state)-1,0)
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
   try:
      lock()   
      r = _cmd.fit("((%s) in (%s))" % (str(a),str(b)),
                  "((%s) in (%s))" % (str(b),str(a)),2)
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
   try:
      lock()   
      r = _cmd.fit("((%s) in (%s))" % (str(a),str(b)),
                  "((%s) in (%s))" % (str(b),str(a)),1)
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
   try:
      lock()   
      r = _cmd.fit("((%s) in (%s))" % (str(a),str(b)),
                  "((%s) in (%s))" % (str(b),str(a)),0)
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




