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

if __name__=='pymol.moving':
    
    import thread
    import string

    import selector
    import pymol

    import cmd
    from cmd import _cmd,Shortcut, \
          toggle_dict,toggle_sc, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error        

    def accept(_self=cmd):
        '''

DESCRIPTION

    "accept" is an internal method for handling of session file
    security.

        '''
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)
            r = _cmd.accept(_self._COb)
            _self.set_wizard()
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def decline(_self=cmd):
        '''
DESCRIPTION

    "decline" is an internal method for handling of session file
    security.

        '''
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)
            r = _cmd.decline(_self._COb)
            _self.set_wizard()
        finally:
            _self.unlock(r,_self)

    def get_movie_playing(_self=cmd):
        '''
DECRIPTION

    "get_movie_playing" returns a boolean value informing the caller
    as to whether or not the movie is currently playing.
    
        '''
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)
            r = _cmd.get_movie_playing(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r
    
    def mdump(_self=cmd):
        '''
DESCRIPTION

    "mdump" dumps the current set of movie commands as text output.

USAGE

    mdump

SEE ALSO

    mplay, mset, mdo, mclear, mmatrix
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.mdump(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def mtoggle(_self=cmd):
        '''
DESCRIPTION

    "mtoggle" toggles playing of the movie.
    
    '''        
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)   
            r = _cmd.mplay(_self._COb,-1)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r
    
    def mstop(_self=cmd):
        '''
DESCRIPTION

    "mstop" stops playing of the movie.

USAGE

    mstop

SEE ALSO

    mplay, mset, mdo, mclear, mmatrix
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.mplay(_self._COb,0)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r


    mview_action_dict = {
        'store'       : 0,
        'clear'       : 1,
        'interpolate' : 2,
        'reinterpolate' : 3,
        'smooth'   : 4,
        'reset'    : 5,
        }

    mview_action_sc = Shortcut(mview_action_dict.keys())

    def mview(action='store', first=0, last=0, power=1.4,
              bias=1.0, simple=0, linear=0.0, object='',
              wrap=-1, hand=1, window=5, cycles=1, scene='',
              cut=0.5, quiet=1, _self=cmd):

        '''
DESCRIPTION

    "mview" stores camera and object matrices for use in movie
    interpolation.

USAGE

    mview [action [, first [, last [, power [, bias [, simple
       [, linear [, object [, wrap [, hand [, window [, cycles [,scene
       [, cut [, quiet ]]]]]]]]]]]]]]]

NOTES

    This command is not yet fully supported in PyMOL 1.x.
    
SEE ALSO

    mplay, mset, mdo, mclear, mmatrix
        '''
        
        r = DEFAULT_ERROR
        first = int(first)
        last = int(last)
        if first<0:
            first = _self.count_frames() + first + 1
            if last == 0:
                last = _self.count_frames()
        if last<0:
            last = _self.count_frames() + last + 1
        action = mview_action_dict[mview_action_sc.auto_err(action,'action')]
        if scene==None:
            scene = _self.get("scene_current_name")
        scene = str(scene)
        if (scene!=''):
            cmd.scene(scene,"recall","",animate=0,frame=0)
        try:
            _self.lock(_self)
            r = _cmd.mview(_self._COb,int(action),int(first)-1,int(last)-1,
                           float(power),float(bias),
                           int(simple), float(linear),str(object),
                           int(wrap),int(hand),int(window),int(cycles),
                           str(scene),float(cut),int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r
    
    def mplay(_self=cmd):
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
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.mplay(_self._COb,1)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def mray(_self=cmd): # deprecated
        '''
DESCRIPTION

    "mray" is an unsupported command of unknown function.
    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.mplay(_self._COb,2)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def mdo(frame,command,_self=cmd):
        '''
DESCRIPTION

    "mdo" defines (or redefines) the command-line operations
    associated with a particular movie frame.  These "generalized
    movie commands" will be executed every time the numbered frame is
    played.

USAGE

    mdo frame: command

PYMOL API

    cmd.mdo( int frame, string command )

EXAMPLE

    // Creates a single frame movie involving a rotation about X and Y

    load test.pdb
    mset 1
    mdo 1, turn x,5; turn y,5;
    mplay

NOTES

 These commands are usually created
    by a PyMOL utility program (such as util.mrock).  Command can
    actually contain several commands separated by semicolons ';'
    
    The "mset" command must first be used to define the movie before
    "mdo" statements will have any effect.  Redefinition of the movie
    clears any existing mdo statements.

SEE ALSO

    mset, mplay, mstop
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.mdo(_self._COb,int(frame)-1,str(command),0)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def mappend(frame,command,_self=cmd):
        '''
DESCRIPTION

    "mappend" associates additional command line operations with a
    particular movie frame.  These "generalized movie commands" will
    be executed every time the numbered frame is played.
    
USAGE

    mappend frame: command

ARGUMENTS

    frame = integer: the frame to modify

    command = literal command-line text
    
EXAMPLE

    mappend 1: hide everything; show sticks
    mappend 60: hide sticks; show spheres
    mappend 120: hide spheres; show surface
    
NOTES

    The "mset" command must first be used to define the movie before
    "mdo" statements will have any effect.  Redefinition of the movie
    clears any existing movie commands specified with mdo or mappend.

SEE ALSO

    mset, madd, mdo, mplay, mstop
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.mdo(_self._COb,int(frame)-1,str(";"+command),1)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def mpng(prefix,first=0,last=0,preserve=0,modal=0,
             mode=-1,quiet=1,_self=cmd):
        '''
DESCRIPTION

    "mpng" writes movie frames as a series of numbered png files.

USAGE

    mpng prefix [, first [, last ]]

ARGUMENTS

    prefix = string: filename prefix for saved images -- output files
    will be numbered and end in ".png"

    first = integer: starting frame {default: 0 (first frame)}

    last = integer: last frame {default: 0 (last frame)}

    modal = integer: will frames be rendered with a modal draw loop
    
NOTES

    If the "ray_trace_frames" variable is non-zero, then the frames
    will be ray-traced.  Note that this can take many hours for a long
    movie with complex content displayed.

    Also, be sure to avoid setting "cache_frames" when rendering a
    long movie to avoid running out of memory.
    
    Arguments "first" and "last" can be used to specify an inclusive
    interval over which to render frames.  Thus, you can write a smart
    Python program that will automatically distribute rendering over a
    cluster of workstations.  If these options are left at zero, then
    the entire movie will be rendered.

PYMOL API

    cmd.mpng(string prefix, int first, int last)

SEE ALSO

    png, save
    
        '''
        r = DEFAULT_ERROR
        if thread.get_ident() ==pymol.glutThread:
            r = _self._mpng(prefix,int(first)-1,int(last)-1,
                            int(preserve),int(modal),-1,int(mode),int(quiet))
        else:
            r = _self.do('cmd._mpng("""'+prefix+'""","'+
                         str(int(first)-1)+'","'+
                         str(int(last)-1)+'","'+
                         str(int(preserve))+'","'+
                         str(int(modal))+'","'+
                         str(int(-1))+'","'+
                         str(int(mode))+'","'+
                         str(int(quiet))+
                         '",_self=cmd)',0)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def mclear(_self=cmd):
        '''
DESCRIPTION

    "mclear" clears the movie frame image cache.

USAGE

    mclear

PYMOL API

    cmd.mclear()
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.mclear(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r


    def frame(frame,_self=cmd):
        '''
DESCRIPTION

    "frame" sets the viewer to the indicated movie frame.

USAGE

    frame frame

ARGUMENTS

    frame = integer: frame number to display

EXAMPLE

    frame 10
    
PYMOL API

    cmd.frame( int frame_number )

NOTES

    Frame numbers are 1-based.

SEE ALSO

    count_states
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.frame(_self._COb,int(frame))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def madd(specification=""):
        '''
DESCRIPTION

    "madd" extends the existing movie specification using the same
    syntax as mset.

SEE ALSO

    mset, mdo, mplay, mclear

    '''
        mset(specification,0)
        
    def mset(specification="",frame=1,_self=cmd):
        '''
DESCRIPTION

    "mset" sets up a relationship between molecular states and movie
    frames.  This makes it possible to control which states are shown
    in which frame.

USAGE

    mset specification [ ,frame ]

PYMOL API

    cmd.mset( string specification [, int frame] )

EXAMPLES

    # simplest case, one state -> one frame

    mset 1

    # ten frames, all corresponding to state 1
    
    mset 1 x10     

    # the first thirty frames are state 1
    # the next 15 frames pass through states 1-15
    # the next 30 frames are of state 15
    # the next 15 frames iterate back to state 1

    mset 1 x30 1 -15 15 x30 15 -1

SEE ALSO

    mdo, mplay, mclear
        '''
        r = DEFAULT_ERROR
        cur_state = _self.get_state()-1 # use the current state 
        try:
            _self.lock(_self)
            output=[]
            input = string.split(string.strip(specification))
            last = -1
            for x in input:
                if x[0]>"9" or x[0]<"0":
                    if x[0]=="x":
                        if last<0:
                            last = cur_state
                            cnt = int(x[1:])                     
                        else:
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
            r = _cmd.mset(_self._COb,string.join(output," "),int(frame)-1)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r


    def mmatrix(action,_self=cmd):
        '''
DESCRIPTION

    "mmatrix" sets up a matrix to be used for the first frame of the movie.

USAGE

    mmatrix action

ARGUMENTS

    action = clear, store, or recall

NOTES

    This command ensures that the movie always starts from the same
    camera view.

    "mmatrix" should not be used when controlling the camera using
    "mview".
    
PYMOL API

    cmd.mmatrix( string action )

EXAMPLES

    mmatrix store
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            if action=="clear":
                r = _cmd.mmatrix(_self._COb,0)
            elif action=="store":
                r = _cmd.mmatrix(_self._COb,1)
            elif action=="recall":
                r = _cmd.mmatrix(_self._COb,2)
            elif action=="check":
                r = _cmd.mmatrix(_self._COb,3)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r


    def forward(_self=cmd):
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
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.setframe(_self._COb,5,1)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def backward(_self=cmd):
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
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.setframe(_self._COb,5,-1)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r


    def rewind(_self=cmd):
        '''
DESCRIPTION

    "rewind" goes to the beginning of the movie.

USAGE

    rewind

PYMOL API

    cmd.rewind()
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.setframe(_self._COb,4,0)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def ending(_self=cmd):
        '''
DESCRIPTION

    "ending" goes to the end of the movie.

USAGE

    ending

PYMOL API

    cmd.ending()
        '''
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)   
            r=_cmd.setframe(_self._COb,6,0)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def middle(_self=cmd):
        '''
DESCRIPTION

    "middle" goes to the middle of the movie.

USAGE

    middle

PYMOL API

    cmd.middle()
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _cmd.setframe(_self._COb,3,0)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r


    def get_state(_self=cmd):
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
        r = _cmd.get_state(_self._COb)+1
        return r

    def get_frame(_self=cmd):
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
        r = _cmd.get_frame(_self._COb)
        return r
