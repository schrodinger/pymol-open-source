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

if __name__=='pymol.controlling':
   
   import string
   import selector
   import cmd
   from cmd import _cmd,lock,unlock,Shortcut,QuietException

   button_code = {
      'left' : 0,
      'middle' : 1,
      'right' : 2,
      'scroll' : 3,
      'double_left' : 4,
      'double_middle' : 5,
      'double_right' : 6,
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
      'none' : 22 ,
      'cent' : 23 ,
      'pktb' : 24 ,
      'slab' : 25 ,
      'movs' : 26 ,
      
      'pk1'  : 27 ,
      'mova' : 28 ,
      'menu' : 29 ,
      }
   but_act_sc = Shortcut(but_act_code.keys())

   ring_dict = {
      'three_button' : [   'three_button_viewing',
                           'three_button_editing' ],
      'two_button' : [ 'two_button_viewing',
                       'two_button_selecting',
                               ],
      'two_button_editing' : [ 'two_button_viewing',
                               'two_button_selecting',
                               'two_button_editing',
                               ]
      }

   def config_mouse(mode='three_button',quiet=1):
      global mouse_ring
      if ring_dict.has_key(mode):
         mouse_ring = ring_dict[mode]
         if not quiet:
            print " config_mouse: %s"%mode
         mouse(quiet=1)

      else:
         print " Error: unrecognized mouse ring: '%s'"%mode

   mouse_ring = ring_dict['three_button']

   mode_dict = {
      'three_button_viewing' : [ ('l','none','rota'),
                                 ('m','none','move'),
                                 ('r','none','movz'),
                                 ('l','shft','+lbx'),
                                 ('m','shft','-lbx'),
                                 ('r','shft','clip'),                 
                                 ('l','ctrl','+lb'),
                                 ('m','ctrl','pkat'),
                                 ('r','ctrl','pkbd'),                 
                                 ('l','ctsh','lb'),
                                 ('m','ctsh','cent'),
                                 ('r','ctsh','rb'),
                                 ('s','none','slab'),
                                 ('s','shft','movs'),
                                 ('s','ctrl','none'),
                                 ('s','ctsh','none'),
                                 ('double_left','none','menu'),
                                 ('double_middle','none','cent'),
                                 ('double_right','none', 'none'),
                                 ],

      'three_button_editing': [ ('l','none','rota'),
                                ('m','none','move'),
                                ('r','none','movz'),
                                ('l','shft','rotf'),
                                ('m','shft','movf'),
                                ('r','shft','clip') ,                 
                                ('l','ctrl','torf'),
                                ('m','ctrl','pkat'),
                                ('r','ctrl','pktb'),                  
                                ('l','ctsh','mova'),
                                ('m','ctsh','cent'),
                                ('r','ctsh','menu'),
                                ('s','none','slab'),
                                ('s','shft','movs'),
                                ('s','ctrl','none'),
                                ('s','ctsh','none'),
                                ],

      'two_button_viewing' : [ ('l','none','rota'),
                               ('m','none','none'),
                               ('r','none','movz'),
                               ('l','shft','pkat'),
                               ('m','shft','none'),
                               ('r','shft','clip'),                 
                               ('l','ctrl','move'),
                               ('m','ctrl','none'),
                               ('r','ctrl','pkbd'),                 
                               ('l','ctsh','lb'),
                               ('m','ctsh','none'),
                               ('r','ctsh','cent'),
                               ('s','none','none'),
                               ('s','shft','none'),
                               ],
      'two_button_selecting' : [ ('l','none','rota'),
                                 ('m','none','none'),
                                 ('r','none','movz'),
                                 ('l','shft','+lbx'),
                                 ('m','shft','none'),
                                 ('r','shft','-lbx'),                 
                                 ('l','ctrl','+lb'),
                                 ('m','ctrl','none'),
                                 ('r','ctrl','+rb'),                 
                                 ('l','ctsh','lb'),
                                 ('m','ctsh','none'),
                                 ('r','ctsh','rb'),
                                 ('s','none','none'),
                                 ('s','shft','none'),
                                 ('s','ctrl','none'),
                                 ('s','ctsh','none'),
                               ],
      'two_button_editing' : [ ('l','none','rota'),
                               ('m','none','none'),
                               ('r','none','movz'),
                               ('l','shft','pkma'),
                               ('m','shft','none'),
                               ('r','shft','clip'),                 
                               ('l','ctrl','torf'),
                               ('m','ctrl','none'),
                               ('r','ctrl','pktb'),                 
                               ('l','ctsh','rotf'),
                               ('m','ctsh','none'),
                               ('r','ctsh','movf'),
                               ('s','none','none'),
                               ('s','shft','none'),
                               ('s','ctrl','none'),
                               ('s','ctsh','none'),
                               ],
      }

   def mouse(action=None,quiet=1):# INTERNAL
      # NOTE: PyMOL automatically runs this routine upon start-up
      try:
         lock()

         if action=='forward':
            bm = _cmd.get_setting("button_mode")
            bm = (int(bm) + 1) % len(mouse_ring)
            cmd.set("button_mode",str(bm),quiet=1)
            action=None
         elif action=='backward':
            bm = _cmd.get_setting("button_mode")
            bm = (int(bm) - 1) % len(mouse_ring)
            cmd.set("button_mode",str(bm),quiet=1)
            action=None

         mode_list = None
         if action==None:
            bm = _cmd.get_setting("button_mode")
            bm = int(bm) % len(mouse_ring)
            mode = mouse_ring[bm]
            mode_list = mode_dict[mode]
         elif action in mode_dict.keys():
            mode = action
            mode_list = mode_dict[mode]
         if mode_list!=None:
            for a in mode_list:
               apply(button,a)
            if not quiet:
               print " mouse: %s"%mode
      finally:
         unlock()

   def edit_mode(mode=None):
      # legacy function
      if len(mouse_ring):
         mouse_mode = mouse_ring[0]
         if mouse_mode[0:10]=='two_button':
            mouse(action='two_button_editing',quiet=0)
         elif mouse_mode[0:12] == 'three_button':
            mouse(action='three_button_editing',quiet=0)

   #   try:
   #      lock()
   #      r = _cmd.get_setting("button_mode")
   #      r = int(r)
   #      if mode==None:
   #         if r:
   #            _cmd.legacy_set("button_mode","0")
   #         else:
   #            _cmd.legacy_set("button_mode","1")            
   #      else:
   #         if mode=='on':
   #            _cmd.legacy_set("button_mode","1")
   #         if mode=='off':
   #            _cmd.legacy_set("button_mode","0")
   #      config_mouse()
   #   finally:
   #      unlock()
   #   pass

   def set_key(key,fn,arg=(),kw={}):  
      '''
DESCRIPTION

   "set_key" binds a specific python function to a key press.

PYMOL API (ONLY)

   cmd.set_key( string key, function fn, tuple arg=(), dict kw={})

PYTHON EXAMPLE

   from pymol import cmd

   def color_blue(object): cmd.color("blue",object)

   cmd.set_key( 'F1' , make_it_blue, ( "object1" ) )
   // would turn object1 blue when the F1 key is pressed and

   cmd.set_key( 'F2' , make_it_blue, ( "object2" ) )
   // would turn object2 blue when the F2 key is pressed.

   cmd.set_key( 'CTRL-C' , cmd.zoom )   
   cmd.set_key( 'ALT-A' , cmd.turn, ('x',90) )

KEYS WHICH CAN BE REDEFINED

   F1 to F12
   left, right, pgup, pgdn, home, insert
   CTRL-A to CTRL-Z 
   ALT-0 to ALT-9, ALT-A to ALT-Z

SEE ALSO

   button
      '''
      r = 0
      if key[0:5]=='CTRL-':
         pat=key[5:]
         if len(pat)>1: # ctrl-special key
            if pat[0]!='F':
               pat=string.lower(pat)
            for a in cmd.ctrl_special.keys():
               if cmd.ctrl_special[a][0]==pat:
                  cmd.ctrl_special[a][1]=fn
                  cmd.ctrl_special[a][2]=arg
                  cmd.ctrl_special[a][3]=kw
                  r = 1
         else: # std. ctrl key
            for a in cmd.ctrl.keys():
               if a==pat:
                  cmd.ctrl[a][0]=fn
                  cmd.ctrl[a][1]=arg
                  cmd.ctrl[a][2]=kw
                  r = 1
      elif key[0:4]=='ALT-':
         pat=key[4:]
         if len(pat)>1: # alt-special key
            if pat[0]!='F':
               pat=string.lower(pat)
            for a in cmd.alt_special.keys():
               if cmd.alt_special[a][0]==pat:
                  cmd.alt_special[a][1]=fn
                  cmd.alt_special[a][2]=arg
                  cmd.alt_special[a][3]=kw
                  r = 1
         else: # std. alt key
            pat=string.lower(pat)
            for a in cmd.alt.keys():
               if a==pat:
                  cmd.alt[a][0]=fn
                  cmd.alt[a][1]=arg
                  cmd.alt[a][2]=kw
                  r = 1
      elif key[0:5]=='SHFT-':
         pat=key[5:]
         if len(pat)>1: # shft-special key
            if pat[0]!='F':
               pat=string.lower(pat)
            for a in cmd.shft_special.keys():
               if cmd.shft_special[a][0]==pat:
                  cmd.shft_special[a][1]=fn
                  cmd.shft_special[a][2]=arg
                  cmd.shft_special[a][3]=kw
                  r = 1
      else:
         if key[0]!='F':
            pat=string.lower(key)
         else:
            pat=key
         for a in cmd.special.keys():
            if cmd.special[a][0]==pat:
               cmd.special[a][1]=fn
               cmd.special[a][2]=arg
               cmd.special[a][3]=kw
               r = 1

      if not r:
         print "Error: special '%s' key not found."%key
         if cmd._raising(): raise QuietException
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

   button:      L, M, R, S
   modifers:    None, Shft, Ctrl, CtSh
   actions:     Rota, Move, MovZ, Clip, RotZ, ClpN, ClpF
                lb,   mb,   rb,   +lb,  +lbX, -lbX, +mb,  +rb, 
                PkAt, PkBd, RotF, TorF, MovF, Orig, Cent

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
         button_num = button_code[button]
         if button_num<3:
            but_code = button_num + 3*but_mod_code[modifier]
         elif button_num<4:
            but_code = 12 + but_mod_code[modifier]
         else:
            but_code = 15 + button_num - 3
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

