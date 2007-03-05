#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information. 
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-* Peter Haebel, Byron DeLaBarre
#-* 
#-*
#Z* -------------------------------------------------------------------

import cmd
import glob
import math
import string

def sweep(pause=0,cycles=1):
    pause = int(pause)
    cycles = int(cycles)
    n_state = cmd.count_states("all")
    if pause>0:
        pass_string = "1 x%d 1 -%d %d x%d %d -1"%(pause, n_state, n_state, pause, n_state)
    else:
        pass_string = "1 -%d %d -1"%(n_state, n_state)
    movie_list = [ pass_string ] * cycles
    movie_string = string.join(movie_list," ")
    cmd.mset(movie_string)

def pause(pause=15,cycles=1):
    pause = int(pause)
    cycles = int(cycles)
    n_state = cmd.count_states("all")
    if pause>0:
        pass_string = "1 x%d 1 -%d %d x%d"%(pause, n_state, n_state, pause)
    else:
        pass_string = "1 -%d %d -1"%(n_state, n_state)
    movie_list = [ pass_string ] * cycles
    movie_string = string.join(movie_list," ")
    cmd.mset(movie_string)
                
def load(*args,**kw):
    nam = "mov"
    if len(args)>1:
        nam = args[1]
    fils = glob.glob(args[0])
    fils.sort()
    if not len(fils):
        print "Error: no matching files"
    else:
        for a in fils:
            apply(cmd.load,(a,nam),kw)
#         cmd.load(a,nam)

def rock(first,last,angle=30,phase=0,loop=1,axis='y'):
    first=int(first)
    last=int(last)
    angle=float(angle)
    phase=float(phase)
    loop=int(loop)
    nstep = (last-first)+1
    if nstep<0:
        nstep = 1
    if loop:
        subdiv = nstep
    else:
        subdiv = nstep+1
    ang_cur = math.pi*phase/180
    ang_inc = 2*math.pi/subdiv
    ang_cur = ang_cur - ang_inc
    a = 0
    while a<nstep:
        last = angle*math.sin(ang_cur)/2
        ang_cur = ang_cur + ang_inc
        disp = angle*math.sin(ang_cur)/2
        diff = disp-last
        # com = "mdo %d:turn %s,%8.3f" % (first+a,axis,diff)
        # cmd.do(com)
        cmd.mdo("%d"%(first+a),"turn %s,%8.3f"% (axis,diff))      
        a = a + 1

def roll(first,last,loop=1,axis='y'):
    first=int(first)
    last=int(last)
    loop=int(loop)
    n = last - first
    if loop:
        step = 2*math.pi/(n+1)
    else:
        step = 2*math.pi/n   
    a = 0
    deg = (180*step/math.pi)
    while a<=n:
        # com = "mdo %d:turn %s,%8.3f" % (first+a,axis,deg)
        # cmd.do(com)
        cmd.mdo("%d" % (first+a), "turn %s,%8.3f" % (axis,deg))
        a = a + 1

def tdroll(first,rangex,rangey,rangez,skip=1):
    '''
AUTHOR

    Byron DeLaBarre

USAGE

    movie.tdroll(rangx,rangey,rangez,skip=1,mset=0)
    
    rangex/y/z = rotation range on respective axis
    enter 0 for no rotation.

    skip is angle increment in each frame
    
    Use skip to reduce final movie size or to speed up rotation.
    
EXAMPLE

    movie.tdroll 360,360,360,5
    
'''
    first = int(first)
    rangex = float(rangex)
    rangey = float(rangey)
    rangez = float(rangez)
    skip = int(skip)
    axis=['x','y','z']
    rangel=[rangex,rangey,rangez]
    axpos=0
    frpos=1
    for ax in axis:
        range = int(rangel[axpos])
        if range:
            leftover = divmod(range,skip)
            print leftover[1]
            if leftover[1]:
              range = range + int(leftover[1])
            a = 0
            while a<range:
                cmd.mdo("%d" % (first+frpos-1), "turn %s,%8.3f" % (ax,skip))
                a = a + skip
                frpos = frpos + 1
            axpos = axpos + 1
        else:
            axpos = axpos + 1
    print (" tdroll: defined rotations for", frpos - 1,
             "frames, starting at frame %d"%first)

def zoom(first,last,step=1,loop=1,axis='z'):
    # Author: Peter Haebel
    first=int(first)
    last=int(last)
    step=int(step)
    loop=int(loop)
    n = last - first
    a = 0
    while a<=n:
        if (loop and a>n/2):
            s = -step
        else:
            s = step
        # com = "mdo %d:move %s,%8.3f" % (first+a,axis,s)
        # cmd.do(com)
        cmd.mdo("%d" % (first+a),"move %s,%8.3f" % (axis,s))
        a = a + 1

def nutate(first,last,angle=30,phase=0,loop=1,shift=math.pi/2.0,factor=0.01):
    first=int(first)
    last=int(last)
    angle=float(angle)
    phase=float(phase)
    loop=int(loop)
    nstep = (last-first)+1
    if nstep<0:
        nstep = 1
    if loop:
        subdiv = nstep
    else:
        subdiv = nstep+1
    ang_cur = math.pi*phase/180
    ang_inc = 2*math.pi/subdiv
    ang_cur = ang_cur - ang_inc
    a = 0
    while a<nstep:
        lastx = angle*math.sin(ang_cur)/2
        lasty = angle*math.sin(ang_cur+shift)/2
        ang_cur = ang_cur + ang_inc
        nextx = angle*math.sin(ang_cur)/2
        nexty = angle*math.sin(ang_cur+shift)/2      
        # com = "mdo %d:turn %s,%8.3f" % (first+a,axis,diff)
        # cmd.do(com)
        cmd.mdo("%d"%(first+a),"turn x,%8.3f;turn y,%8.3f;turn y,%8.3f;turn x,%8.3f"%
                  (-lastx,-lasty,nexty,nextx))
        a = a + 1

def screw(first,last,step=1,angle=30,phase=0,loop=1,axis='y'):
    # Author: Peter Haebel
    first=int(first)
    last=int(last)
    step=int(step)
    angle=float(angle)
    phase=float(phase)
    loop=int(loop)
    nstep = (last-first)+1
    if nstep<0:
        nstep = 1
    if loop:
        subdiv = nstep
    else:
        subdiv = nstep+1
    ang_cur = math.pi*phase/180
    ang_inc = 2*math.pi/subdiv
    ang_cur = ang_cur - ang_inc
    a = 0
    while a<nstep:
        if (loop and a>=nstep/2):
            s = -step
        else:
            s = step
        last = angle*math.sin(ang_cur)/2
        ang_cur = ang_cur + ang_inc
        disp = angle*math.sin(ang_cur)/2
        diff = disp-last
        # com = "mdo %d:turn %s,%8.3f; move z,%8.3f" % (first+a,axis,diff,s)
        # cmd.do(com)
        cmd.mdo("%d" % (first+a), "turn %s,%8.3f; move z,%8.3f" % (axis,diff,s))
        a = a + 1

def timed_roll(period=12.0,cycles=1,axis='y'):
    frames_per_sec = float(cmd.get('movie_fps'))
    if frames_per_sec<1.0:
        frames_per_sec=30.0
    frames_per_cycle = int(period*frames_per_sec)
    total = frames_per_cycle * cycles
    
    cmd.mset("1 x%d"%total)
    step = 2*math.pi/(frames_per_cycle)
    deg = (180*step/math.pi)
    cmd.mview('reset')
    cmd.rewind()
    frame = 1
    for cycle in range(cycles):
        for cnt in range(frames_per_cycle):
            cmd.turn(axis,deg)
            cmd.mview('store',frame)
            frame = frame + 1
    
    
