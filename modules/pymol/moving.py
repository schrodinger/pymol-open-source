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
   from cmd import _cmd,lock,unlock,Shortcut,QuietException
   from cmd import toggle_dict,toggle_sc


   def accept():
      '''
      SECURITY FEATURE
      '''
      try:
         lock()
         r = _cmd.accept()
         cmd.set_wizard()
      finally:
         unlock()

   def decline():
      '''
      SECURITY FEATURE
      '''
      try:
         lock()
         r = _cmd.decline()
         cmd.set_wizard()
      finally:
         unlock()

   def mdump():
      '''
DESCRIPTION

   "mdump" dumps the current set of movie commands

USAGE

   mdump

PYMOL API

   cmd.mdump()

SEE ALSO

   mplay, mset, mdo, mclear, mmatrix
      '''
      try:
         lock()   
         r = _cmd.mdump(0)
      finally:
         unlock()
      return r

   def mstop():
      '''
DESCRIPTION

   "mstop" stops the movie.

USAGE

   mstop

PYMOL API

   cmd.mstop()

SEE ALSO

   mplay, mset, mdo, mclear, mmatrix
      '''
      try:
         lock()   
         r = _cmd.mplay(0)
      finally:
         unlock()
      return r


   mview_action_dict = {
      'store'       : 0,
      'clear'       : 1,
      'interpolate' : 2,
      }

   mview_action_sc = Shortcut(mview_action_dict.keys())

   def mview(action='store',first=0,last=0,power=1.4):
      action = mview_action_dict[mview_action_sc.auto_err(action,'action')]
      try:
         lock()   
         r = _cmd.mview(int(action),int(first)-1,int(last)-1,float(power))
      finally:
         unlock()
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
      try:
         lock()   
         r = _cmd.mplay(1)
      finally:
         unlock()
      return r

   def mray(): # deprecated
      try:
         lock()   
         r = _cmd.mplay(2)
      finally:
         unlock()
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
      try:
         lock()   
         r = _cmd.mdo(int(frame)-1,str(command),0)
      finally:
         unlock()
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
      try:
         lock()   
         r = _cmd.mdo(int(frame)-1,str(";"+command),1)
      finally:
         unlock()
      return r

   def mpng(prefix,first=0,last=0):
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
      if thread.get_ident() ==pymol.glutThread:
         r = cmd._mpng(prefix,int(first)-1,int(last)-1)
      else:
         r = _cmd.do('cmd._mpng("'+prefix+'","'+
                     str(int(first)-1)+'","'+str(int(last)-1)+'")',0)
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
      try:
         lock()   
         r = _cmd.mclear()
      finally:
         unlock()
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
      try:
         lock()   
         r = _cmd.frame(int(frame))
      finally:
         unlock()
      return r

   def mset(specification=""):
      '''
DESCRIPTION

   "mset" sets up a relationship between molecular states and movie
   frames.  This makes it possible to control which states are shown
   in which frame.

USAGE

   mset specification

PYMOL API

   cmd.mset( string specification )

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
      try:
         lock()
         output=[]
         input = string.split(string.strip(specification))
         last = 0
         for x in input:
            if x[0]>"9" or x[0]<"0":
               if x[0]=="x":
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
         r = _cmd.mset(string.join(output," "))
      finally:
         unlock()
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
      r = 1
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
         unlock()
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
      try:
         lock()   
         r = _cmd.setframe(5,1)
      finally:
         unlock()
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
      try:
         lock()   
         r = _cmd.setframe(5,-1)
      finally:
         unlock()
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
      try:
         lock()   
         r = _cmd.setframe(4,0)
      finally:
         unlock()
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
      try:
         lock()   
         r=_cmd.setframe(6,0)
      finally:
         unlock()
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
      try:
         lock()   
         r = _cmd.setframe(3,0)
      finally:
         unlock()
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
