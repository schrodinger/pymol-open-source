#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* Copyright (c) Schrodinger, LLC. 
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

# abstract (external or internal) gui control interface

if __name__=='pymol.gui':
    
    import pymol
    import cmd

# external gui control 

def ext_hide(_self=cmd):
    pymol = _self._pymol
    if pymol._ext_gui != None:
        pymol._ext_gui.fifo.put('self.root.withdraw()')
    else:
        pass
    
def ext_show(_self=cmd):
    pymol = _self._pymol
    if pymol._ext_gui != None:
        pymol._ext_gui.fifo.put('self.root.deiconify()')
    else:
        pass

# common actions

def save_as(_self=cmd):
    pymol = _self._pymol
    if pymol._ext_gui != None:
        pymol._ext_gui.fifo.put('self.skin.session_save_as()')
    else:
        pass

def save_image(_self=cmd):
    pymol = _self._pymol
    if pymol._ext_gui != None:
        pymol._ext_gui.fifo.put('self.skin.file_save_png()')
    else:
        pass
    
