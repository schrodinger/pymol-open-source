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

import string
import selector

from cmd import _cmd,lock,unlock,Shortcut,QuietException

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

def config_mouse(quiet=0):# INTERNAL
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
      r = _cmd.mask("("+str(selection)+")",1)
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
      r = _cmd.mask("("+str(selection)+")",0)
   finally:
      unlock()
   return r
