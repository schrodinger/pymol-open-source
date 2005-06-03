import _cmd
import threading
import pymol

P_GLUT_IDLE_EVENT          =  0
P_GLUT_DISPLAY_EVENT       =  1
P_GLUT_RESHAPE_EVENT       =  2
P_GLUT_MOUSE_EVENT         =  3
P_GLUT_MOTION_EVENT        =  4
P_GLUT_CHAR_EVENT          =  5
P_GLUT_SPECIAL_EVENT       =  6
P_GLUT_PASSIVE_MOTION_EVENT=  7
P_GLUT_ACTIVE_ALT          =  4
P_GLUT_ACTIVE_CTRL         =  2
P_GLUT_ACTIVE_SHIFT        =  1

P_GLUT_LEFT_BUTTON         =  0
P_GLUT_MIDDLE_BUTTON       =  1
P_GLUT_RIGHT_BUTTON        =  2
P_GLUT_BUTTON_SCROLL_FORWARD= 3
P_GLUT_BUTTON_SCROLL_BACKWARD=4

P_GLUT_DOWN                =  0
P_GLUT_UP                  =  1

P_GLUT_KEY_F1         =  1
P_GLUT_KEY_F2         =  2
P_GLUT_KEY_F3         =  3
P_GLUT_KEY_F4         =  4
P_GLUT_KEY_F5         =  5
P_GLUT_KEY_F6         =  6
P_GLUT_KEY_F7         =  7
P_GLUT_KEY_F8         =  8
P_GLUT_KEY_F9         =  9
P_GLUT_KEY_F10        =  10
P_GLUT_KEY_F11        =  11
P_GLUT_KEY_F12        =  12
P_GLUT_KEY_LEFT       =  100
P_GLUT_KEY_UP         =  101
P_GLUT_KEY_RIGHT      =  102
P_GLUT_KEY_DOWN       =  103
P_GLUT_KEY_PAGE_UP    =  104
P_GLUT_KEY_PAGE_DOWN  =  105
P_GLUT_KEY_HOME       =  106
P_GLUT_KEY_END        =  107
P_GLUT_KEY_INSERT     =  108

special_dict = {
    'F1'   :        P_GLUT_KEY_F1,
    'F2'   :        P_GLUT_KEY_F2,
    'F3'   :        P_GLUT_KEY_F3,
    'F4'   :        P_GLUT_KEY_F4,
    'F5'   :        P_GLUT_KEY_F5,
    'F6'   :        P_GLUT_KEY_F6,
    'F7'   :        P_GLUT_KEY_F7,
    'F8'   :        P_GLUT_KEY_F8,
    'F9'   :        P_GLUT_KEY_F9,
    'F10'  :        P_GLUT_KEY_F10,
    'F11'  :        P_GLUT_KEY_F11,
    'F12'  :        P_GLUT_KEY_F12,
    'LEFT' :        P_GLUT_KEY_LEFT,
    'UP'   :        P_GLUT_KEY_UP,
    'RIGHT'     :   P_GLUT_KEY_RIGHT,
    'DOWN'      :   P_GLUT_KEY_DOWN,
    'PAGE_UP'   :   P_GLUT_KEY_PAGE_UP,
    'PAGE_DOWN' :   P_GLUT_KEY_PAGE_DOWN,
    'HOME'      :   P_GLUT_KEY_HOME,
    'END'       :   P_GLUT_KEY_END,
    'INSERT'    :   P_GLUT_KEY_INSERT,
}

mod_dict = {}

def get_mod_value(shift,control,meta):
    global mod_dict
    mod = 0
    if shift:
        mod = mod + P_GLUT_ACTIVE_SHIFT
    if control:
        mod = mod + P_GLUT_ACTIVE_CTRL
    if meta:
        mod = mod + P_GLUT_ACTIVE_ALT
    mod_dict[(shift,control,meta)]=mod
    return mod

class EmbeddedPyMOL:

    def ep_get_pymol(self):
        return pymol
    
    def ep_swap_dummy(self):
        pass
     
    def ep_init(self):
        _cmd.runwxpymol()
        pymol._swap_buffers = lambda s=self: s.ep_swap_dummy() # dummy swap function
        self.ep_mod = 0
        self.ep_button = None
        self.ep_swap = None
        
    def ep_reshape(self, width, height):
        _cmd.runwxpymol() 
        _cmd.p_glut_event(P_GLUT_RESHAPE_EVENT,width,height,0,0,0)

    def ep_char(self, x, y, code, shift, control, meta):
        self.ep_mod = mod_dict.get((shift,control,meta))
        if self.ep_mod == None:
            self.ep_mod = get_mod_value(shift,control,meta)
        _cmd.runwxpymol() 
        _cmd.p_glut_event(P_GLUT_CHAR_EVENT,x,y,code,0,self.ep_mod)

    def ep_special(self,x,y,code,shift,control,meta):
        self.ep_mod = mod_dict.get((shift,control,meta))
        if self.ep_mod == None:
            self.ep_mod = get_mod_value(shift,control,meta)
        code = special_dict.get(code)
        if code!=None:
            _cmd.runwxpymol() 
            _cmd.p_glut_event(P_GLUT_SPECIAL_EVENT,x,y,code,0,self.ep_mod)
    
    def ep_mouse_down(self,x,y,left,middle,right,shift,control,meta):
        self.ep_mod = mod_dict.get((shift,control,meta))
        if self.ep_mod == None:
            self.ep_mod = get_mod_value(shift,control,meta)
        if left:
            self.ep_button = P_GLUT_LEFT_BUTTON
        elif middle:
            self.ep_button = P_GLUT_MIDDLE_BUTTON
        elif right:
            self.ep_button = P_GLUT_RIGHT_BUTTON
        _cmd.runwxpymol() 
        _cmd.p_glut_event(P_GLUT_MOUSE_EVENT,x,y,self.ep_button,P_GLUT_DOWN,self.ep_mod)
        
    def ep_mouse_up(self, x, y):
        if self.ep_button != None:
            _cmd.runwxpymol() 
            _cmd.p_glut_event(P_GLUT_MOUSE_EVENT,x,y,self.ep_button,P_GLUT_UP,self.ep_mod)
        self.ep_button = None
        
    def ep_motion(self, x, y, left, middle, right, shift, control, meta):
        self.ep_mod = mod_dict.get((shift,control,meta))
        if self.ep_mod == None:
            self.ep_mod = get_mod_value(shift,control,meta)
        if left:
            self.ep_button = P_GLUT_LEFT_BUTTON
        elif middle:
            self.ep_button = P_GLUT_MIDDLE_BUTTON
        elif right:
            self.ep_button = P_GLUT_RIGHT_BUTTON
        _cmd.runwxpymol() 
        _cmd.p_glut_event(P_GLUT_MOTION_EVENT,x,y,self.ep_button,0,self.ep_mod)

    def ep_passive_motion(self, x, y, shift, control, meta):
        self.ep_mod = mod_dict.get((shift,control,meta))
        if self.ep_mod == None:
            self.ep_mod = get_mod_value(shift,control,meta)
        _cmd.runwxpymol() 
        _cmd.p_glut_event(P_GLUT_PASSIVE_MOTION_EVENT,x,y,0,0,self.ep_mod)

    def ep_wheel(self, x, y, direction, shift, control, meta):
        _cmd.runwxpymol()
        if direction>0:
            _cmd.p_glut_event(P_GLUT_MOUSE_EVENT, x,y,P_GLUT_BUTTON_SCROLL_FORWARD,0,self.ep_mod)
        else:
            _cmd.p_glut_event(P_GLUT_MOUSE_EVENT, x,y,P_GLUT_BUTTON_SCROLL_BACKWARD,0,self.ep_mod)
        
    def ep_draw(self):
        _cmd.runwxpymol() 
        _cmd.p_glut_event(P_GLUT_DISPLAY_EVENT,0,0,0,0,0) # draw event

    def ep_idle(self):
        _cmd.runwxpymol() 
        _cmd.p_glut_event(0,0,0,0,0,0)

    def ep_get_redisplay(self):
        _cmd.runwxpymol()
        result = _cmd.p_glut_get_redisplay()
        return result
    
    def ep_set_swap_callback(self,swap):
        self.ep_swap = swap
        pymol._swap_buffers = lambda s=self: s.swap()
        

