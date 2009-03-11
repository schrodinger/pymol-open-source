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
    import pymol
    
    from cmd import _cmd, Shortcut, QuietException, is_string, \
          boolean_dict, boolean_sc, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, is_ok, is_error, \
          location_code, location_sc 

    button_code = {
        'left' : 0,
        'middle' : 1,
        'right' : 2,
        'wheel' : 3,
        'double_left' : 4,
        'double_middle' : 5,
        'double_right' : 6,
        'single_left' : 7,
        'single_middle' : 8,
        'single_right' : 9
        }
    button_sc = Shortcut(button_code.keys())

    but_mod_code = {
        'none'  : 0,
        'shft'  : 1,
        'ctrl'  : 2,
        'ctsh'  : 3,
	'alt'   : 4,
	'alsh'  : 5,
	'ctal'  : 6,
	'ctas'  : 7,
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

        'sele' : 30 ,
        '+/-'  : 31 ,
        '+box' : 32 ,
        '-box' : 33 ,

        'mvsz' : 34 ,
        'dgrt' : 36 ,
        'dgmv' : 37 ,
        'dgmz' : 38 ,

        'roto' : 39 ,
        'movo' : 40 ,
        'mvoz' : 41 ,
        'mvfz' : 42 ,
        'mvaz' : 43 ,
        'drgm' : 44 ,

        'rotv' : 45 ,
        'movv' : 46 ,
        'mvvz' : 47 ,

        # 48 for internal use only
        
        'drgo' : 49 ,

        'imsz' : 50 ,
        'imvz' : 51 ,
        'box'  : 52 ,
        'irtz' : 53
        }
    
    but_act_sc = Shortcut(but_act_code.keys())

    ring_dict = {
        'maestro' : [   'three_button_maestro' ],
        
        'three_button' : [   'three_button_viewing',
                             'three_button_editing' ],
        'two_button' : [ 'two_button_viewing',
                         'two_button_selecting',
                         ],

        'two_button_editing' : [ 'two_button_viewing',
                                 'two_button_selecting',
                                 'two_button_editing',
                                 ],
        'three_button_motions' : [   'three_button_viewing',
                                     'three_button_editing',
                                     'three_button_motions' ],
        'one_button' : [   'one_button_viewing' ],
        }

    def config_mouse(ring='three_button', quiet=1, _self=cmd):

        '''
DESCRIPTION

    "config_mouse" sets the current mouse configuration ring.

USAGE

    config_mouse ring

EXAMPLES

    config_mouse three_button
    config_mouse two_button
    config_mouse one_button
    
PYMOL API

    cmd.config_mouse(string ring, int quiet)

SEE ALSO

    mouse, button
    '''
        
        global mouse_ring
        if ring_dict.has_key(ring):
            mouse_ring = ring_dict[ring]
            if not quiet:
                print " config_mouse: %s"%ring
            _self.mouse(quiet=1,_self=_self)

        else:
            print " Error: unrecognized mouse ring: '%s'"%ring

    mouse_ring = ring_dict['three_button']

    mode_name_dict = {
        'three_button_maestro' : '3-Button Maestro',        
        'three_button_viewing' : '3-Button Viewing',
        'three_button_editing' : '3-Button Editing',
        'three_button_motions':  '3-Button Motions',
        'two_button_viewing'   : '2-Button Viewing',
        'two_button_selecting' : '2-Btn. Selecting',
        'two_button_editing'   : '2-Button Editing',
        'one_button_viewing'   : '1-Button Viewing',
        }

    mode_dict = {
        'three_button_maestro' : [ ('l','none','box'),
				   ('m','none','rota'),
				   ('r','none','move'),
				   ('l','shft','+Box'),
				   ('m','shft','-Box'),
				   ('r','shft','clip'),                 
				   ('l','ctrl','+/-'),
				   ('m','ctrl','irtz'),
				   ('r','ctrl','pk1'),                 
				   ('l','ctsh','Sele'),
				   ('m','ctsh','orig'),
				   ('r','ctsh','clip'),
				   ('l','alt' ,'none'),
				   ('m','alt' ,'none'),
				   ('r','alt' ,'none'),
				   ('w','none','imvz'),
				   ('w','shft','movs'),
				   ('w','ctrl','none'), # disable since ctrl-middle is irtz
				   ('w','ctsh','slab'),
				   ('double_left','none','menu'),
				   ('double_middle','none','none'),
				   ('double_right','none', 'pkat'),
				   ('single_left','none','sele'),
				   ('single_middle','none','cent'),
				   ('single_right','none', 'menu'),
                                   ('single_left','shft','+/-'),
				   ],
          'three_button_viewing' : [ ('l','none','rota'),
				   ('m','none','move'),
				   ('r','none','movz'),
				   ('l','shft','+Box'),
				   ('m','shft','-Box'),
				   ('r','shft','clip'),                 
				   ('l','ctrl','+/-'),
				   ('m','ctrl','pkat'),
				   ('r','ctrl','pk1'),                 
				   ('l','ctsh','Sele'),
				   ('m','ctsh','orig'),
				   ('r','ctsh','clip'),
				   ('l','alt' ,'none'),
				   ('m','alt' ,'none'),
				   ('r','alt' ,'none'),
				   ('w','none','slab'),
				   ('w','shft','movs'),
				   ('w','ctrl','mvsz'),
				   ('w','ctsh','movz'),
				   ('double_left','none','menu'),
				   ('double_middle','none','none'),
				   ('double_right','none', 'pkat'),
				   ('single_left','none','+/-'),
				   ('single_middle','none','cent'),
				   ('single_right','none', 'menu'),
				   ],
        'three_button_editing': [ ('l','none','rota'),
				  ('m','none','move'),
				  ('r','none','movz'),
				  ('l','shft','roto'),
				  ('m','shft','movo'),
				  ('r','shft','mvoz') ,                 
				  ('l','ctrl','torf'),
				  ('m','ctrl','+/-'),
				  ('r','ctrl','pktb'),                  
				  ('l','ctsh','mova'),
				  ('m','ctsh','orig'),
				  ('r','ctsh','clip'),
				  ('l','alt' ,'none'),
				  ('m','alt' ,'none'),
				  ('r','alt' ,'none'),
				  ('w','none','slab'),
				  ('w','shft','movs'),
				  ('w','ctrl','mvsz'),
				  ('w','ctsh','movz'),
				  ('double_left','none','torf'),
				  ('double_middle','none','drgm'),
				  ('double_right','none', 'pktb'),
				  ('single_left','none','pkat'),
				  ('single_middle','none','cent'),
				  ('single_right','none', 'menu'),
				  ],
	
        'three_button_motions': [ ('l','none','rota'),
				  ('m','none','move'),
				  ('r','none','movz'),
				  ('l','shft','rotv'),
				  ('m','shft','movv'),
				  ('r','shft','mvvz') ,                 
				  ('l','ctrl','torf'),
				  ('m','ctrl','pkat'),
				  ('r','ctrl','pktb'),                  
				  ('l','ctsh','mova'),
				  ('m','ctsh','orig'),
				  ('r','ctsh','clip'),
				  ('l','alt' ,'none'),
				  ('m','alt' ,'none'),
				  ('r','alt' ,'none'),
				  ('w','none','slab'),
				  ('w','shft','movs'),
				  ('w','ctrl','mvsz'),
				  ('w','ctsh','movz'),
				  ('double_left','none','menu'),
				  ('double_left','none','torf'),
				  ('double_middle','none','drgm'),
				  ('double_right','none', 'pktb'),
				  ('single_left','none','pkat'),
				  ('single_middle','none','cent'),
				  ('single_right','none', 'menu'),
				  ],

        'two_button_viewing' : [ ('l','none','rota'),
				 ('m','none','none'),
				 ('r','none','movz'),
				 ('l','shft','pk1'),
				 ('m','shft','none'),
				 ('r','shft','clip'),                 
				 ('l','ctrl','move'),
				 ('m','ctrl','none'),
				 ('r','ctrl','pkat'),                 
				 ('l','ctsh','sele'),
				 ('m','ctsh','none'),
				 ('r','ctsh','cent'),
				 ('l','alt' ,'none'),
				 ('m','alt' ,'none'),
				 ('r','alt' ,'none'),
				 ('w','none','none'),
				 ('w','shft','none'),
				 ('w','ctrl','none'),
				 ('w','ctsh','none'),
				 ('double_left','none','menu'),
				 ('double_middle','none','none'),                               
				 ('double_right','none','cent'),
				 ('single_left','none','pkat'),
				 ('single_middle','none','none'),
				 ('single_right','none', 'menu'),
				 ],
        'two_button_selecting' : [ ('l','none','rota'),
				   ('m','none','none'),
				   ('r','none','movz'),
				   ('l','shft','+Box'),
				   ('m','shft','none'),
				   ('r','shft','-Box'),                 
				   ('l','ctrl','+/-'),
				   ('m','ctrl','none'),
				   ('r','ctrl','pkat'),                 
				   ('l','ctsh','sele'),
				   ('m','ctsh','none'),
				   ('r','ctsh','cent'),
				   ('l','alt' ,'none'),
				   ('m','alt' ,'none'),
				   ('r','alt' ,'none'),
				   ('w','none','none'),
				   ('w','shft','none'),
				   ('w','ctrl','none'),
				   ('w','ctsh','none'),
				   ('double_left','none','menu'),
				   ('double_left','none','menu'),
				   ('double_middle','none','none'),                               
				   ('double_right','none','cent'),
				   ('single_left','none','+/-'),
				   ('single_right','none', 'menu'),
				   ],
        'two_button_editing' : [ ('l','none','rota'),
				 ('m','none','none'),
				 ('r','none','movz'),
				 ('l','shft','pkat'),
				 ('m','shft','none'),
				 ('r','shft','clip'),                 
				 ('l','ctrl','torf'),
				 ('m','ctrl','none'),
				 ('r','ctrl','pktb'),                 
				 ('l','ctsh','rotf'),
				 ('m','ctsh','none'),
				 ('r','ctsh','movf'),
				 ('l','alt' ,'none'),
				 ('m','alt' ,'none'),
				 ('r','alt' ,'none'),
				 ('w','none','none'),
				 ('w','shft','none'),
				 ('w','ctrl','none'),
				 ('w','ctsh','none'),
				 ('double_left','none','menu'),
				 ('double_middle','none','none'),                               
				 ('double_right','none','cent'),
				 ('single_left','none','pkat'),
				 ('single_middle','none','none'),                               
				 ('single_right','none','menu'),                               
				 ],
         'one_button_viewing' : [ ('l','none','rota'), # approximate the old GLUT behavior on Mac
				  ('m','none','none'),
				  ('r','none','none'),
				  ('l','shft','+Box'),
				  ('m','shft','none'),
				  ('r','shft','none'),
				  ('l','ctrl','movZ'),
				  ('m','ctrl','none'),
				  ('r','ctrl','none'),
				  ('l','ctsh','clip'),
				  ('m','ctsh','none'),
				  ('r','ctsh','none'),
				  ('l','alt' ,'move'),
				  ('m','alt' ,'none'),
				  ('r','alt' ,'none'),
				  ('l','alsh' ,'-Box'),
				  ('m','alsh' ,'none'),
				  ('r','alsh' ,'none'),
				  ('l','ctal' ,'none'),
				  ('m','ctal' ,'none'),
				  ('r','ctal' ,'none'),
				  ('l','ctas' ,'none'),
				  ('m','ctas' ,'none'),
				  ('r','ctas' ,'none'),
				  ('w','none','none'),
				  ('w','shft','none'),
				  ('w','ctrl','none'),
				  ('w','ctsh','none'),
				  ('double_left','none','menu'),
				  ('double_middle','none','none'),
				  ('double_right','none', 'none'),
				  ('single_left','none','+/-'),
				  ('single_middle','none','none'),
				  ('single_right','none', 'none'),

				  ('single_left','shft','none'),
				  ('single_left','ctrl','menu'),
				  ('single_left','ctsh','pkat'), 
				  ('single_left','alt', 'cent'),
				  ('single_left','alsh','none'),
				  ('single_left','ctal','none'),
				  ('single_left','ctas','none'),

				  ],
        }

    def order(names,sort=0,location='current',_self=cmd):
        '''
DESCRIPTION

    "order" changes the ordering of names in the control panel.

USAGE

    order names, sort, location

ARGUMENTS

    names = string: a space-separated list of names

    sort = yes or no {default: no}

    location = top, current, or bottom {default: current}

EXAMPLES

    order 1dn2 1fgh 1rnd       # sets the order of these three objects
    order *,yes                # sorts all names
    order 1dn2_*, yes          # sorts all names beginning with 1dn2_
    order 1frg, location=top   # puts 1frg at the top of the list

PYMOL API

    cmd.order(string names, string sort, string location)

NOTES

    "order" can also be used to reorder objects within a group.
    
SEE ALSO

    set_name, group
        '''

        r = DEFAULT_ERROR
        location=location_code[location_sc.auto_err(location,'location')]
        if is_string(sort):
            sort=boolean_dict[boolean_sc.auto_err(sort,'sort option')]
        try:
            _self.lock(_self)
            r = _cmd.order(_self._COb,str(names),int(sort),int(location))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException         
        return r

    def mouse(action=None, quiet=1, _self=cmd):# INTERNAL

        '''
DESCRIPTION

    "mouse" cycles through the mouse modes defined in the current
    mouse configuration ring.

USAGE

    mouse 

'''
        # NOTE: PyMOL automatically runs this routine upon start-up
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)

            if action=='forward':
                bm = _cmd.get_setting(_self._COb,"button_mode")
                bm = (int(bm) + 1) % len(mouse_ring)
                _self.set("button_mode",str(bm),quiet=1)
                action=None
            elif action=='backward':
                bm = _cmd.get_setting(_self._COb,"button_mode")
                bm = (int(bm) - 1) % len(mouse_ring)
                _self.set("button_mode",str(bm),quiet=1)
                action=None
            elif action=='select_forward':
                sm = _cmd.get_setting(_self._COb,"mouse_selection_mode")
                sm = sm + 1
                if sm>6: sm = 0
                _self.set("mouse_selection_mode",sm,quiet=1)
            elif action=='select_backward':
                sm = _cmd.get_setting(_self._COb,"mouse_selection_mode")
                sm = sm - 1
                if sm<0: sm = 6
                _self.set("mouse_selection_mode",sm,quiet=1)
            
            mode_list = None
            if action==None:
                bm = _cmd.get_setting(_self._COb,"button_mode")
                bm = int(bm) % len(mouse_ring)
                mode = mouse_ring[bm]
                _self.set("button_mode_name",mode_name_dict.get(mode,mode))
                mode_list = mode_dict[mode]
            elif action in mode_dict.keys():
                mode = action
                _self.set("button_mode_name",mode_name_dict.get(mode,mode))
                if mode in mouse_ring:
                    bm = mouse_ring.index(mode)
                    _self.set("button_mode",bm)
                mode_list = mode_dict[mode]
            if mode_list!=None:
                kw_dict = {'_self':_self}
                for a in mode_list:
                    apply(button,a,kw_dict)
                if not quiet:
                    print " mouse: %s"%mode
                if mode[-7:]!='editing': _self.unpick()
                if mode[-7:]=='editing': _self.deselect()
            r = DEFAULT_SUCCESS
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                     
        return r
            

    def edit_mode(active=1,quiet=1,_self=cmd):
        '''
DESCRIPTION

    "edit_mode" switches the mouse into editing mode, if such a mode
    is available in the current mouse ring.
    
    '''
        # legacy function
        if is_string(active):
            active=boolean_dict[boolean_sc.auto_err(active,'active')]
        active = int(active)
        if len(mouse_ring):
            bm = int(_cmd.get_setting(_self._COb,"button_mode"))
            mouse_mode = mouse_ring[bm]
            if active:
                if mouse_mode[0:10]=='two_button':
                    if mouse_mode!='two_button_editing':
                        mouse(action='two_button_editing',quiet=quiet,_self=_self)
                elif mouse_mode[0:12] == 'three_button':
                    if mouse_mode!='three_button_editing':
                        mouse(action='three_button_editing',quiet=quiet,_self=_self)
            else:
                if mouse_mode[0:10]=='two_button':
                    if mouse_mode!='two_button_viewing':               
                        mouse(action='two_button_viewing',quiet=quiet,_self=_self)
                elif mouse_mode[0:12] == 'three_button':
                    if mouse_mode!='three_button_viewing':
                        mouse(action='three_button_viewing',quiet=quiet,_self=_self)
        return DEFAULT_SUCCESS
    
    def set_key(key,fn,arg=(),kw={},_self=cmd):  
        '''
DESCRIPTION

    "set_key" binds a specific python function to a key press.

PYMOL API (ONLY)

    cmd.set_key( string key, function fn, tuple arg=(), dict kw={})

PYTHON EXAMPLE

    from pymol import cmd

    def color_blue(object): cmd.color("blue",object)

    cmd.set_key( 'F1' , color_blue, ( "object1" ) )
    // would turn object1 blue when the F1 key is pressed and

    cmd.set_key( 'F2' , color_blue, ( "object2" ) )
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
        r = DEFAULT_ERROR
        if key[0:5]=='CTRL-':
            pat=key[5:]
            if len(pat)>1: # ctrl-special key
                if pat[0]!='F':
                    pat=string.lower(pat)
                for a in _self.ctrl_special.keys():
                    if _self.ctrl_special[a][0]==pat:
                        _self.ctrl_special[a][1]=fn
                        _self.ctrl_special[a][2]=arg
                        _self.ctrl_special[a][3]=kw
                        r = DEFAULT_SUCCESS
            else: # std. ctrl key
                for a in _self.ctrl.keys():
                    if a==pat:
                        _self.ctrl[a][0]=fn
                        _self.ctrl[a][1]=arg
                        _self.ctrl[a][2]=kw
                        r = DEFAULT_SUCCESS
        elif key[0:4]=='ALT-':
            pat=key[4:]
            if len(pat)>1: # alt-special key
                if pat[0]!='F':
                    pat=string.lower(pat)
                for a in _self.alt_special.keys():
                    if _self.alt_special[a][0]==pat:
                        _self.alt_special[a][1]=fn
                        _self.alt_special[a][2]=arg
                        _self.alt_special[a][3]=kw
                        r = DEFAULT_SUCCESS
            else: # std. alt key
                pat=string.lower(pat)
                for a in _self.alt.keys():
                    if a==pat:
                        _self.alt[a][0]=fn
                        _self.alt[a][1]=arg
                        _self.alt[a][2]=kw
                        r = DEFAULT_SUCCESS
        elif key[0:5]=='SHFT-':
            pat=key[5:]
            if len(pat)>1: # shft-special key
                if pat[0]!='F':
                    pat=string.lower(pat)
                for a in _self.shft_special.keys():
                    if _self.shft_special[a][0]==pat:
                        _self.shft_special[a][1]=fn
                        _self.shft_special[a][2]=arg
                        _self.shft_special[a][3]=kw
                        r = DEFAULT_SUCCESS
        else:
            if key[0]!='F':
                pat=string.lower(key)
            else:
                pat=key
            for a in _self.special.keys():
                if _self.special[a][0]==pat:
                    _self.special[a][1]=fn
                    _self.special[a][2]=arg
                    _self.special[a][3]=kw
                    r = DEFAULT_SUCCESS
        if is_error(r):
            print "Error: special '%s' key not found."%key
        if _self._raising(r,_self): raise pymol.CmdException         
        return r

    def button(button, modifier, action, _self=cmd):
        '''
DESCRIPTION

    "button" can be used to redefine what the mouse buttons do.

USAGE

    button button, modifier, action

ARGUMENTS

    button = left, middle, right, wheel, double_left, double_middle,
        double_right, single_left, single_middle, or single_right
       
    modifiers = None, Shft, Ctrl, CtSh, CtAl, CtAl, CtAS, 

    actions = None, Rota, Move, MovZ, Slab, +Box, -Box, Clip, MovS,
        +/-, PkAt, Pk1, MvSZ, Sele, Orig, Menu, PkAt, Pk1 RotO, MovO,
        MvOZ, MovA, PkAt, PkTB, MvSZ MvAZ, DrgM, RotZ, PkBd, ClpN,
        ClpF

NOTES

   Changes made using the button command are easily overridden when
   the user iterates through the mouse modes.  This behavior needs to
   be changed.

   Obsolete actions: lb, mb, rb, +lb, +mb, +rb, +lbX, -lbX,
      
   Unsupported, Internal, or Future Actions: RotD, MovD, MvDZ, RotF,
    MovF, MvFZ, TorF, RotV, MovV, MvVZ, DgMZ, DgRT

PYMOL API

    cmd.button(string button, string modifier, string action)

SEE ALSO

    config_mouse
    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            button = string.lower(button)
            button = button_sc.auto_err(button,'button')
            modifier = string.lower(modifier)
            modifier = but_mod_sc.auto_err(modifier,'modifier')
            action = string.lower(action)
            action = but_act_sc.auto_err(action,'action')
            button_num = button_code[button]
	    but_mod_num = but_mod_code[modifier]
            if button_num<3: # normal button (L,M,R)
		if but_mod_num < 4: # none, shft, ctrl, ctsh
		    but_code = button_num + 3*but_mod_num
		else: # alt, alsh, alct, alcs
		    but_code = button_num + 68 + 3*(but_mod_num-4)
            elif button_num < 4: # wheel
		if but_mod_num < 4: # none, shft, ctrl, ctsh
		    but_code = 12 + but_mod_num
		else:
		    but_code = 64 + but_mod_num - 4
            else: # single and double clicks
                but_code = (16 + button_num - 4) + but_mod_num * 6
            act_code = but_act_code[action]
            r = _cmd.button(_self._COb,but_code,act_code)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException         
        return r


    def mask(selection="(all)",quiet=1,_self=cmd):
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
        r = DEFAULT_ERROR
        # preprocess selection
        selection = selector.process(selection)
        #
        try:
            _self.lock(_self)   
            r = _cmd.mask(_self._COb,"("+str(selection)+")",1,int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException         
        return r

    def unmask(selection="(all)",quiet=1,_self=cmd):
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
        r = DEFAULT_ERROR
        # preprocess selection
        selection = selector.process(selection)
        #   
        try:
            _self.lock(_self)   
            r = _cmd.mask(_self._COb,"("+str(selection)+")",0,int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException         
        return r

