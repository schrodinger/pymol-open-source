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

# pmgui_tk 
# TKinter based gui for PyMol
# NOTE: must have treads support compiled into python to use this module
#
# **This is the only module which should be/need be imported by 
# **PyMol Programs

from PMGApp import *
import sys, os, threading

if sys.platform=='win32':
    if sys.version[0:4]=='2.1 ':
        if not os.environ.has_key('TCL_LIBRARY'):
            os.environ['TCL_LIBRARY']='c:\\python21\\tcl\\tcl8.3'
                
def run(pymol_instance,poll=0):
    if not hasattr(sys,"argv"):
        sys.argv=["pymol"]
    PMGApp(pymol_instance).run(poll)

def __init__(pymol_instance,poll=0):

    t = threading.Thread(target=run,args=(pymol_instance,poll))
    t.setDaemon(1)
    t.start()





