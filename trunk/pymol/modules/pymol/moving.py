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
    from cmd import _cmd,lock,unlock,Shortcut, \
          toggle_dict,toggle_sc, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error        

    def accept():
        '''

DESCRIPTION

    "accept" is an internal method for handling of session file
    security.

        '''
        r = DEFAULT_ERROR      
        try:
            lock()
            r = _cmd.accept()
            cmd.set_wizard()
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException
        return r

    def decline():
        '''
DESCRIPTION

    "decline" is an internal method for handling of session file
    security.

        '''
        r = DEFAULT_ERROR      
        try:
            lock()
            r = _cmd.decline()
            cmd.set_wizard()
        finally:
            unlock(r)

    def get_movie_playing():
        '''
DECRIPTION

    "get_movie_playing" returns a boolean value informing the caller
    as to whether or not the movie is currently playing.
    
        '''
        r = DEFAULT_ERROR      
        try:
            lock()
            r = _cmd.get_movie_playing()
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException
        return r
    
    def mdump():
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
            lock()   
            r = _cmd.mdump(0)
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException
        return r

    def mtoggle():
        '''
DESCRIPTION

    "mtoggle" toggles playing of the movie.
    
    '''        
        r = DEFAULT_ERROR      
        try:
            lock()   
            r = _cmd.mplay(-1)
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException
        return r
    
    def mstop():
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
            lock()   
            r = _cmd.mplay(0)
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException
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

    def mview(action='store',first=0,last=0,power=1.4,
              bias=1.0,simple=0,linear=0.0,object='',
              wrap=-1,hand=1,window=5,cycles=1,scene='',
              cut=0.5,quiet=1):
        r = DEFAULT_ERROR
        first = int(first)
        last = int(last)
        if first<0:
            first = cmd.count_frames() + first + 1
            if last == 0:
                last = cmd.count_frames()
        if last<0:
            last = cmd.count_frames() + last + 1
        action = mview_action_dict[mview_action_sc.auto_err(action,'action')]
        if scene==None:
            scene = ''
        else:
            scene = cmd.get("scene_current_name")
        try:
            lock()
            r = _cmd.mview(int(action),int(first)-1,int(last)-1,
                           float(power),float(bias),
                           int(simple), float(linear),str(object),
                           int(wrap),int(hand),int(window),int(cycles),
                           str(scene),float(cut),int(quiet))
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException
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
        r = DEFAULT_ERROR
        try:
            lock()   
            r = _cmd.mplay(1)
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException
        return r

    def mray(): # deprecated
        r = DEFAULT_ERROR
        try:
            lock()   
            r = _cmd.mplay(2)
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException
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
        r = DEFAULT_ERROR
        try:
            lock()   
            r = _cmd.mdo(int(frame)-1,str(command),0)
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException
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
        r = DEFAULT_ERROR
        try:
            lock()   
            r = _cmd.mdo(int(frame)-1,str(";"+command),1)
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException
        return r

    def mpng(prefix,first=0,last=0,preserve=0):
        '''
DESCRIPTION

    "mpng" writes a series of numbered movie frames to png files with
    the specified prefix.  If the "ray_trace_frames" variable is
    non-zero, these frames will be ray-traced.  This operation can take
    several hours for a long movie.

    Be sure to disable "cache_frames" when issuing this operation on a
    long movie (>100 frames) to avoid running out of memory.

USAGE

    mpng prefix [, first [, last]]

    Options "first" and "last" can be used to specify an inclusive
    interval over which to render frames.  Thus, you can write a smart
    Python program that will automatically distribute rendering over a
    cluster of workstations.  If these options are left at zero, then
    the entire movie will be rendered.

PYMOL API

    cmd.mpng( string prefix, int first=0, int last=0 )
        '''
        r = DEFAULT_ERROR
        if thread.get_ident() ==pymol.glutThread:
            r = cmd._mpng(prefix,int(first)-1,int(last)-1,int(preserve))
        else:
            r = cmd.do('cmd._mpng("'+prefix+'","'+
                            str(int(first)-1)+'","'+str(int(last)-1)+','+str(int(preserve))+'")',0)
        if _raising(r): raise pymol.CmdException
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
        r = DEFAULT_ERROR
        try:
            lock()   
            r = _cmd.mclear()
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException
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
        r = DEFAULT_ERROR
        try:
            lock()   
            r = _cmd.frame(int(frame))
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException
        return r

    def madd(specification=""):
        mset(specification,0)
        
    def mset(specification="",frame=1):
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
        r = DEFAULT_ERROR
        cur_state = cmd.get_state()-1 # use the current state 
        try:
            lock()
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
            r = _cmd.mset(string.join(output," "),int(frame)-1)
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException
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
        r = DEFAULT_ERROR
        try:
            lock()   
            if action=="clear":
                r = _cmd.mmatrix(0)
            elif action=="store":
                r = _cmd.mmatrix(1)
            elif action=="recall":
                r = _cmd.mmatrix(2)
            elif action=="check":
                r = _cmd.mmatrix(3)
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException
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
        r = DEFAULT_ERROR
        try:
            lock()   
            r = _cmd.setframe(5,1)
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException
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
        r = DEFAULT_ERROR
        try:
            lock()   
            r = _cmd.setframe(5,-1)
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException
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
        r = DEFAULT_ERROR
        try:
            lock()   
            r = _cmd.setframe(4,0)
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException
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
        r = DEFAULT_ERROR      
        try:
            lock()   
            r=_cmd.setframe(6,0)
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException
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
        r = DEFAULT_ERROR
        try:
            lock()   
            r = _cmd.setframe(3,0)
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException
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
