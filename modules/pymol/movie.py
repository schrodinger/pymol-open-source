#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* Copyright (c) Schrodinger, LLC.
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

import sys
cmd = sys.modules["pymol.cmd"]
import math
import os
import glob
import threading
import time
from . import colorprinting

def get_movie_fps(_self):
    r = _self.get_setting_float('movie_fps')
    if r <= 0:
        return 30
    return r

def sweep(pause=0,cycles=1,_self=cmd):
    pause = int(pause)
    cycles = int(cycles)
    n_state = _self.count_states("all")
    if pause>0:
        pass_string = "1 x%d 1 -%d %d x%d %d -1"%(pause, n_state, n_state, pause, n_state)
    else:
        pass_string = "1 -%d %d -1"%(n_state, n_state)
    movie_list = [ pass_string ] * cycles
    movie_string = " ".join(movie_list)
    _self.mset(movie_string)

def pause(pause=15,cycles=1,_self=cmd):
    pause = int(pause)
    cycles = int(cycles)
    n_state = _self.count_states("all")
    if pause>0:
        pass_string = "1 x%d 1 -%d %d x%d"%(pause, n_state, n_state, pause)
    else:
        pass_string = "1 -%d %d -1"%(n_state, n_state)
    movie_list = [ pass_string ] * cycles
    movie_string = " ".join(movie_list)
    _self.mset(movie_string)

def load(pattern, nam = "mov", _self=cmd, **kw):
    fils = glob.glob(pattern)
    if not fils:
        print("Error: no matching files")
        return
    for a in sorted(fils):
        _self.load(a, nam, **kw)

def rock(first=1,last=-1,angle=30,phase=0,loop=1,axis='y',_self=cmd):
    first=int(first)
    last=int(last)
    if last<0:
        last = _self.count_frames()
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
        # _self.do(com)
        _self.mdo("%d"%(first+a),"turn %s,%8.3f"% (axis,diff))
        a = a + 1

def roll(first=1,last=-1,loop=1,axis='y',_self=cmd):
    first=int(first)
    last=int(last)
    loop=int(loop)
    axis=str(axis)
    if last<0:
        last = _self.count_frames()
    n = last - first
    if loop:
        step = 2*math.pi/(n+1)
    else:
        step = 2*math.pi/n
    a = 0
    invert = 1
    if axis[0:1]=='-':
        axis=axis[1:]
        invert = -1
    deg = invert * (180*step/math.pi)
    while a<=n:
        # com = "mdo %d:turn %s,%8.3f" % (first+a,axis,deg)
        # _self.do(com)
        _self.mdo("%d" % (first+a), "turn %s,%8.3f" % (axis,deg))
        a = a + 1


def tdroll(first,rangex,rangey,rangez,skip=1,_self=cmd):
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
            print(leftover[1])
            if leftover[1]:
              range = range + int(leftover[1])
            a = 0
            while a<range:
                _self.mdo("%d" % (first+frpos-1), "turn %s,%8.3f" % (ax,skip))
                a = a + skip
                frpos = frpos + 1
            axpos = axpos + 1
        else:
            axpos = axpos + 1
    print((" tdroll: defined rotations for", frpos - 1,
             "frames, starting at frame %d"%first))

def zoom(first,last,step=1,loop=1,axis='z',_self=cmd):
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
        # _self.do(com)
        _self.mdo("%d" % (first+a),"move %s,%8.3f" % (axis,s))
        a = a + 1

def nutate(first,last,angle=30,phase=0,loop=1,shift=math.pi/2.0,_self=cmd):
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
        # _self.do(com)
        _self.mdo("%d"%(first+a),"turn x,%8.3f;turn y,%8.3f;turn y,%8.3f;turn x,%8.3f"%
                  (-lastx,-lasty,nexty,nextx))
        a = a + 1

def screw(first,last,step=1,angle=30,phase=0,loop=1,axis='y',_self=cmd):
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
        # _self.do(com)
        _self.mdo("%d" % (first+a), "turn %s,%8.3f; move z,%8.3f" % (axis,diff,s))
        a = a + 1

def timed_roll(period=12.0,cycles=1,axis='y',_self=cmd):
    frames_per_sec = get_movie_fps(_self)
    if frames_per_sec<1.0:
        frames_per_sec=30.0
    frames_per_cycle = int(period*frames_per_sec)
    total = frames_per_cycle * cycles

    _self.mset("1 x%d"%total)
    step = 2*math.pi/(frames_per_cycle)
    deg = (180*step/math.pi)
    _self.mview('reset')
    _self.rewind()
    frame = 1
    for cycle in range(cycles):
        for cnt in range(frames_per_cycle):
            _self.turn(axis,deg)
            _self.mview('store',frame,freeze=1)
            frame = frame + 1

# blank frames to be specified by the user

def add_blank(duration=12.0,start=0,_self=cmd):
    '''
DESCRIPTION

    Add "blank" time to the movie (without any key frames)

ARGUMENTS

    duration = float: time to add to movie in seconds {default: 12}

    start = int: start frame (1 = first frame; 0 = end of movie) {default: 0}

SEE ALSO

    mset
    '''
    cmd = _self
    if not start:
        start = cmd.get_movie_length()+1
    duration = float(duration)
    fps = get_movie_fps(_self)
    n_frame = int(round(fps * duration))
    if n_frame > 0:
        cmd.mset("1 x%d"%n_frame,start)
        cmd.frame(start)

# new matrix-based camera interpolation routines

def add_roll(duration=12.0,loop=1,axis='y',start=0,_self=cmd):
    '''
DESCRIPTION

    Append a 360 degree camera rotation to the movie, using key frames.

ARGUMENTS

    duration = float: time to add to movie in seconds {default: 12}

    loop = 0/1: ???

    axis = x, y or z: rotation axis in camera space {default: y}

    start = int: start frame (1 = first frame; 0 = end of movie) {default: 0}
    '''
    cmd = _self
    if not start:
        start = cmd.get_movie_length()+1
    duration = float(duration)
    fps = get_movie_fps(_self)
    n_frame = int(round(fps * duration))
    if n_frame > 0:
        cmd.mset("1 x%d"%n_frame,start,freeze=1)
        cmd.mview("store",start,power=1,freeze=1)
        cmd.turn(axis,120)
        cmd.mview("store",start+n_frame/3,power=1,freeze=1)
        cmd.turn(axis,120)
        cmd.mview("store",start+(2*n_frame)/3,power=1,freeze=1)
        if loop:
            if (start == 1):
                cmd.mview("interpolate",wrap=1)
                cmd.turn(axis,120)
                cmd.mview("store",start+n_frame-1,power=1,freeze=1)
                cmd.turn(axis,120)
            else:
                adjustment = 360.0/n_frame
                cmd.turn(axis,120 - adjustment)
                cmd.mview("store",start+n_frame-1,power=1,freeze=1)
                cmd.mview("interpolate")
                cmd.turn(axis,adjustment)
        else:
            cmd.turn(axis,120)
            cmd.mview("store",start+n_frame-1,power=1,freeze=1)
            cmd.mview("interpolate")
        cmd.frame(start)
        # PYMOL-2881
        if cmd.get_setting_int('movie_auto_interpolate'):
            cmd.mview("reinterpolate")

def add_rock(duration=8.0,angle=30.0,loop=1,axis='y',start=0,_self=cmd):
    '''
DESCRIPTION

    Append a rocking camera motion to the movie, using key frames.

ARGUMENTS

    duration = float: time to add to movie in seconds {default: 8}

    angle = float: degrees {default: 30}

    loop = 0/1: ???

    axis = x, y or z: rotation axis in camera space {default: y}

    start = int: start frame (1 = first frame; 0 = end of movie) {default: 0}
    '''
    cmd = _self
    if not start:
        start = cmd.get_movie_length()+1
    duration = float(duration)
    angle=float(angle)
    loop=int(loop)
    fps = get_movie_fps(_self)
    n_frame = int(round(fps * duration))
    if n_frame > 0:
        cmd.mset("1 x%d"%n_frame,start,freeze=1)
        cmd.turn(axis,angle/2.0)
        cmd.mview("store",start+n_frame/4,power=-1,freeze=1)
        cmd.turn(axis,-angle)
        cmd.mview("store",start+(3*n_frame)/4,power=-1,freeze=1)
        if loop and (start == 1):
            cmd.mview("interpolate",wrap=1)
        else:
            cmd.mview("interpolate")
        cmd.frame(start)

def add_state_sweep(factor=1,pause=2.0,first=-1,last=-1,loop=1,start=0,_self=cmd):
    cmd = _self
    if not start:
        start = cmd.get_movie_length() + 1
    loop = int(loop)
    fps = get_movie_fps(_self)
    n_state = cmd.count_states()
    duration = (2 * pause) + (2 * factor * n_state) / fps
    n_frame = int(round(fps * duration))
    if n_frame > 0:
        cmd.mset("1 x%d"%n_frame, start, freeze=1)
        cmd.mview("store",start, state=1, freeze=1)
        cmd.mview("store",start + (n_frame * pause) / duration, state=1, freeze=1)
        cmd.mview("store",start + n_state * factor + (n_frame * pause) / duration - 1, state=n_state, freeze=1)
        cmd.mview("store",start + n_state * factor + (2 * n_frame * pause) / duration, state=n_state, freeze=1)
        cmd.mview("store",start + n_frame-1, state=1, freeze=1)
        if loop and (start == 1):
            cmd.mview("interpolate",wrap=1)
        else:
            cmd.mview("interpolate")
        cmd.frame(start)
        # PYMOL-2881
        if cmd.get_setting_int('movie_auto_interpolate'):
            cmd.mview("reinterpolate")

def add_state_loop(factor=1,pause=2.0,first=-1,last=-1,loop=1,start=0,_self=cmd):
    cmd = _self
    if not start:
        start = cmd.get_movie_length() + 1
    loop = int(loop)
    fps = get_movie_fps(_self)
    n_state = cmd.count_states()
    duration = (pause) + (factor * n_state) / fps
    n_frame = int(round(fps * duration))
    if n_frame > 0:
        cmd.mset("1 x%d"%n_frame, start, freeze=1)
        cmd.mview("store",start, state=1, freeze=1)
        cmd.mview("store",start + (n_frame * pause * 0.5) / duration, state=1, freeze=1)
        cmd.mview("store",start + n_state * factor + (n_frame * pause * 0.5) / duration, state=n_state, freeze=1)
        cmd.mview("store",start + n_frame-1, state=n_state, freeze=1)
        if loop and (start == 1):
            cmd.mview("interpolate",wrap=1)
        else:
            cmd.mview("interpolate")
        cmd.frame(start)
        # PYMOL-2881
        if cmd.get_setting_int('movie_auto_interpolate'):
            cmd.mview("reinterpolate")

def add_nutate(duration=8.0, angle=30.0, spiral=0, loop=1,
               offset=0, phase=0, shift=math.pi/2.0, start=0,
               _self=cmd):
    '''
DESCRIPTION

    Append a nutating camera motion to the movie, using key frames.

ARGUMENTS

    duration = float: time to add to movie in seconds {default: 8}

    angle = float: degrees {default: 30}

    spiral = -1/0/1: If zero, do a circular motion, otherwise do a spiral
    motion starting at the center (clockwise if 1, ccw if -1) {default: 0}

    loop: (unused)

    offset: (unused)

    phase = float: phase offset in degrees {default: 0}

    shift = float: x to y offset in radians {default: pi/2}

    start = int: start frame (1 = first frame; 0 = end of movie) {default: 0}
    '''
    cmd = _self
    if not start:
        start = cmd.get_movie_length()+1
    duration = float(duration)
    angle = float(angle)
    spiral = int(spiral)
    loop = int(loop)
    fps = get_movie_fps(_self)
    n_frame = int(round(fps * duration))
    if n_frame > 0:
        cmd.mset("1 x%d"%n_frame,start,freeze=1)
        for index in range(0,n_frame):
            if spiral>0:
                sp_angle = angle*float(index+1)/n_frame
            elif spiral<0:
                sp_angle = angle*(1.0 - float(index+1)/n_frame)
            else:
                sp_angle = angle
            ang_cur = math.pi*phase/180.0 + (2*math.pi*index)/n_frame
            x_rot = sp_angle * math.sin(ang_cur)/2
            y_rot = sp_angle * math.sin(ang_cur+shift)/2
            cmd.turn('x',x_rot)
            cmd.turn('y',y_rot)
            cmd.mview('store',start+index,freeze=1)
            cmd.turn('y',-y_rot)
            cmd.turn('x',-x_rot)
    # PYMOL-2881
    if cmd.get_setting_int('movie_auto_interpolate'):
        cmd.mview("reinterpolate")

def _rock(mode,axis,first,last,period,pause,_self=cmd):
    cmd = _self
    n_frame = last - first + 1
    angle = cmd.get_setting_float('sweep_angle')
    if (period * 1.5) < pause:
        n_cyc = int(round(pause / period))
        frame_list = []
        for cyc in range(n_cyc):
            frame_list.extend( [(first + ((1+4*cyc)*n_frame)/(4*n_cyc)),
                                (first + ((3+4*cyc)*n_frame)/(4*n_cyc))] )
    else:
        frame_list = [ first+n_frame/4, first+(3*n_frame)/4 ]
    if 1 or mode:
        direction = 0
        for frame in frame_list:
            if not direction:
                cmd.turn(axis,angle/2.0)
                cmd.mview("store",frame,power=-1,freeze=1)
                direction = -1
            else:
                cmd.turn(axis, direction * angle)
                cmd.mview("store",frame,power=-1,freeze=1)
                direction = -direction
        cmd.turn(axis,direction * angle/2.0)
        cmd.mview("store",last,power=-1,freeze=1)
        cmd.mview("interpolate",first,last)

def _nutate_sub(start_frame, stop_frame, angle=15.0, spiral=0, loop=1,
                offset=0, phase=0, shift=math.pi/2.0, _self=cmd):
    cmd = _self
    angle = float(angle)
    spiral = int(spiral)
    loop = int(loop)
    fps = get_movie_fps(_self)
    duration = (stop_frame - start_frame)/fps
    n_frame = int(round(fps * duration))
    if n_frame > 0:
        for index in range(0,n_frame):
            if spiral>0:
                sp_angle = angle*float(index+1)/n_frame
            elif spiral<0:
                sp_angle = angle*(1.0 - float(index+1)/n_frame)
            else:
                sp_angle = angle
            ang_cur = math.pi*phase/180.0 + (2*math.pi*index)/n_frame
            x_rot = sp_angle * math.sin(ang_cur)/2
            y_rot = sp_angle * math.sin(ang_cur+shift)/2
            cmd.turn('x',x_rot)
            cmd.turn('y',y_rot)
            cmd.mview('store',start_frame+index,freeze=1)
            cmd.turn('y',-y_rot)
            cmd.turn('x',-x_rot)

def _nutate(mode,first,last,period,pause,_self=cmd):
    cmd = _self
    n_frame = last - first + 1
    axis = 'y'
    angle = cmd.get_setting_float('sweep_angle')
    if (period * 1.5) < pause:
        n_cyc = int(round(pause / period))
        frame_list = []
        for cyc in range(n_cyc):
            frame_list.append( [(first + ((cyc)*n_frame)/(n_cyc)),
                                (first + ((cyc+1)*n_frame)/(n_cyc))] )
    else:
        frame_list = [ [first, first+n_frame] ]
    direction = 0
    spiral = 1
    for frame in frame_list:
        _nutate_sub(frame[0], frame[1], angle, spiral, _self=_self)
        spiral = 0

def add_scenes(names=None, pause=8.0, cut=0.0, loop=1,
               rock=-1, period=8.0, animate=-1, start=0,
               _self=cmd):
    '''
DESCRIPTION

    Append a sequence of scenes to the movie, with camera animation
    (rock/nutate) at each scene.

ARGUMENTS

    names = list: list of scenes names {default: all scenes}

    pause = float: display time per scene in seconds {default: 8}

    cut = float 0.0-1.0: scene switch moment (0.0: beginning of transition,
    1.0: end of transition) {default: 0.0}

    loop = 0/1: end of movie interpolates back to first frame {default: 1}

    rock = int: (sweep_mode + 1): 2=x-axis-rock, 3=y-axis-rock, 4=nutate
    {default: -1, use sweep_mode setting}

    period = float: rock/nutate time (e.g. "period=2, pause=6" will rock 3
    times per scene) {default: 8}

    animate = float: scene transition time in seconds
    {default: -1, use scene_animation_duration setting}

    start = int: start frame (1 = first frame; 0 = end of movie) {default: 0}
    '''
    cmd = _self
    if not start:
        start = cmd.get_movie_length()+1
    animate = float(animate)
    pause = float(pause)
    period = float(period)
    if period>pause:
        period = pause
    rock = int(rock)
    if animate<0:
        animate = float(cmd.get("scene_animation_duration"))
    if names is None:
        names = cmd.get_scene_list()
    elif cmd.is_string(names):
        names = cmd.safe_alpha_list_eval(names)
    n_scene = len(names)
    duration = n_scene*(pause+animate)
    fps = get_movie_fps(_self)
    n_frame = int(round(fps * duration))

    if not loop:
        act_n_frame = int(round(fps * (duration - animate)))
    else:
        act_n_frame = n_frame
    if rock<0:
        sweep_mode = int(cmd.get('sweep_mode'))
    elif rock:
        sweep_mode = rock - 1
    else:
        sweep_mode = 0
    if n_frame > 0:
        cmd.mset("1 x%d"%act_n_frame,start,freeze=1)
        cnt = 0
        for scene in names:
            frame = start+int((cnt*n_frame)/n_scene)
            cmd.mview("store",frame,scene=scene,freeze=1)
            if rock:
                cmd.mview("interpolate",cut=cut,wrap=0)
                sweep_first = frame + 1
                sweep_last = start + int(pause*fps+(cnt*n_frame)/n_scene) - 1
                if sweep_last > act_n_frame:
                    sweep_last = act_n_frame
                n_sweep_frame = sweep_last - sweep_first + 1
                if n_sweep_frame > 0:
                    if sweep_mode==1: # x-axis rock
                        _rock(sweep_mode, 'x', sweep_first, sweep_last,
                                period, pause, _self=_self)
                    elif sweep_mode<3: # y-axis rock
                        _rock(sweep_mode, 'y', sweep_first, sweep_last,
                                period, pause, _self=_self)
                    elif sweep_mode == 3:
                        _nutate(sweep_mode, sweep_first, sweep_last,
                                period, pause, _self=_self)
            frame = start+int(pause*fps+(cnt*n_frame)/n_scene)
            if frame <= act_n_frame:
                if sweep_mode!=3:
                    cmd.mview("store",frame,scene=scene,freeze=1)
            cnt = cnt + 1
        cmd.mview("interpolate",cut=cut,wrap=loop)
        if rock:
            cmd.mview("smooth")
        cmd.frame(start)

_prefix = "mov"

def _watch(filename,done_event):
    tries = 5
    size = 0
    import time,os
    while not os.path.exists(filename):
        if done_event.isSet():
            break
        tries = tries - 1
        if tries < 0:
            break
        time.sleep(1)
    if os.path.exists(filename):
        tries = 5
    while os.path.exists(filename):
        stat = os.stat(filename)
        if size != stat[6]:
            tries = 5
            size = stat[6]
            if done_event.isSet():
                break
            print(" produce: %d bytes written..."%size)
        else:
            tries = tries - 1
            if tries < 0:
                break
        time.sleep(2)
        if done_event.isSet():
            break

def _encode(filename,first,last,preserve,
            encoder,tmp_path,prefix,img_ext,quality,quiet,_self=cmd):
    import os
    tries = 10
    while 1: # loop until all of the files have been created...
        done = 1
        # check for the required output files
        for index in range(first,last+1):
            path = os.path.join(tmp_path, "%s%04d%s" % (prefix, index, img_ext))
            if not os.path.exists(path):
                done = 0
                break;
        if done:
            break
        elif _self.get_modal_draw(): # keep looping so long as we're rendering...
            tries = 10
        else:
            tries = tries - 1
            if tries < 0:
                done = 0
                break
        time.sleep(0.25)
    _self.sync()
    ok = 1
    result = None

    # reduce chance of passing non-ascii file paths to sub processes
    # by changing directory
    fn_rel = os.path.relpath(filename, tmp_path)
    old_cwd = os.getcwd()
    fps = get_movie_fps(_self)

    if done and ok and (encoder == 'mpeg_encode'):
        try:
            from pymol import mpeg_encode
        except:
            ok = 0
            print("produce-error: Unable to import module pymol.mpeg_encode.")
        if ok:
            if not mpeg_encode.validate():
                ok = 0
                print("produce-error: Unable to validate pymol.mpeg_encode.")
        if not ok:
            print("produce-error: Unable to create mpeg file.")
        else:
            mpeg_quality = 1+int(((100-quality)*29)/100) # 1 to 30
            input = mpeg_encode.input(fn_rel, '.',
                                      prefix,first,last,mpeg_quality);

            FPS_LEGAL_VALUES = [23.976, 24, 25, 29.97, 30, 50, 59.94, 60]
            fps_legal = min(FPS_LEGAL_VALUES, key=lambda v: abs(v - fps))
            if fps_legal != round(fps, 3):
                colorprinting.warning(
                    " Warning: Adjusting frame rate to {} fps (legal values are: {})"
                    .format(fps_legal, FPS_LEGAL_VALUES))
            input = input.replace('FRAME_RATE 30',
                                  'FRAME_RATE {:.3f}'.format(fps_legal))

            if not quiet:
                print(" produce: creating '%s' (in background)..."%(filename))

            os.chdir(tmp_path)

            done_event = None
            if not quiet:
                done_event = threading.Event()
                _self.async_(_watch, fn_rel, done_event, _self=_self)

            try:
                result = mpeg_encode.run(input)
            finally:
                os.chdir(old_cwd)
                if done_event is not None:
                    done_event.set()
    elif encoder == 'ffmpeg':
        import subprocess
        os.chdir(tmp_path)
        try:
            args = ['ffmpeg',
                '-f', 'image2',
                '-framerate', '{:.3f}'.format(fps),
                '-i', prefix + '%04d' + img_ext,
            ]
            if fn_rel.endswith('.webm'):
                args_crf = ['-crf', '{:.0f}'.format(65 - (quality / 2))]
                args += ['-c:v', 'libvpx-vp9', '-b:v', '0'] + args_crf
            elif not fn_rel.endswith('.gif'):
                args += [
                '-crf', '10' if quality > 90 else '15' if quality > 80 else '20',
                '-pix_fmt', 'yuv420p', # needed for Mac support
                ]
            process = subprocess.Popen(args + [fn_rel], stderr=subprocess.PIPE)
            stderr = process.communicate()[1]
            colorprinting.warning(stderr.strip().decode(errors='replace'))
            if process.returncode != 0:
                colorprinting.error('ffmpeg failed with '
                        'exit status {}'.format(process.returncode))
        finally:
            os.chdir(old_cwd)
    elif encoder == 'convert':
        import subprocess
        exe = find_exe(encoder)
        try:
            subprocess.check_call([exe,
                '-delay', '{:.3f}'.format(100. / fps), # framerate
                os.path.join(tmp_path, prefix) + '*' + img_ext,
                filename])
        finally:
            pass
    if not quiet:
                if not os.path.exists(filename):
                    if result is not None:
                        print(input, result[0], result[1])
                    print(" produce: compression failed")
                else:
                    print(" produce: finished.")
    _self.unset("keep_alive")
    if preserve<1:
        if os.path.isdir(tmp_path):
            for fil in glob.glob(os.path.join(tmp_path,prefix+"*")):
                os.unlink(fil)
            os.rmdir(tmp_path)

produce_mode_dict = {
    'normal'  : 0,
    'draw'    : 1,
    'ray'     : 2,
    }

produce_mode_sc = cmd.Shortcut(produce_mode_dict.keys())


def find_exe(exe):
    '''Return full path to executable or None.
    Excludes C:\Windows\System32\convert.exe
    Tests .exe extension on Unix (e.g. for legacy "mpeg_encode.exe" name).
    '''
    from shutil import which

    if exe.startswith('convert') and sys.platform == 'win32':
        # filter out C:\Windows\System32
        path = os.pathsep.join(p
                for p in os.getenv('PATH', '').split(os.pathsep)
                if r'\windows\system32' not in p.lower())
        return which(exe, path=path)

    if exe == 'mpeg_encode':
        legacy = which(exe + '.exe')
        if legacy:
            return legacy

    return which(exe)


def produce(filename, mode='', first=0, last=0, preserve=0,
            encoder='', quality=-1, quiet=1,
            width=0, height=0, _self=cmd):
    '''
DESCRIPTION

    Export a movie to an MPEG, WEBM, or GIF file.

    Which video formats and codecs are available depends on the availability
    and feature set of your ffmpeg and ImageMagick installations.

ARGUMENTS

    filename = str: filename of movie file to produce

    mode = draw or ray: {default: check "ray_trace_frames" setting}

    first = int: first frame to export {default: 1}

    last = int: last frame to export {default: last frame of movie}

    preserve = 0 or 1: don't delete temporary files {default: 0}

    encoder = ffmpeg|convert|mpeg_encode: Tool used for video encoding

    quality = 0-100: encoding quality {default: 90 (movie_quality setting)}

    width = int: Width in pixels {default: from viewport}

    height = int: Height in pixels {default: from viewport}

EXAMPLE

    movie.produce video.mp4, height=1080
    movie.produce video.webm, height=720
    '''
    from pymol import CmdException

    def has_exe(exe):
        return bool(find_exe(exe))

    prefix = _prefix

    if mode == '':
        mode = -1
    elif _self.is_string(mode):
        mode = produce_mode_sc.auto_err(mode,"mode")
        mode = produce_mode_dict[mode]
    else:
        mode = int(mode)
    first = int(first)
    last = int(last)
    quiet = int(quiet)
    preserve = int(preserve)
    quality = int(quality)
    if quality<0:
        quality = _self.get_setting_int('movie_quality')
    if quality>100:
        quality = 100
    ok = 1
    splitext = os.path.splitext(filename)
    tmp_path = splitext[0]+".tmp"
    if splitext[1]=='':
        splitext=(splitext[0],'.mpg')
    filename = splitext[0]+splitext[1]
    width, height = int(width), int(height)
    img_ext = '.png'

    # guess encoder
    if not encoder:
        if splitext[1] in ('.mpeg', '.mpg'):
            encoder = 'mpeg_encode'
        elif has_exe('ffmpeg'):
            encoder = 'ffmpeg'
        elif has_exe('convert'):
            encoder = 'convert'
        else:
            raise CmdException('neither "ffmpeg" nor "convert" available for '
                    'video encoding')
        print('using encoder "%s"' % encoder)

    # check encoder
    if encoder == 'mpeg_encode':
        img_ext = '.ppm'
    elif encoder not in ('ffmpeg', 'convert'):
        raise CmdException('unknown encoder "%s"' % encoder)
    elif not has_exe(encoder):
        raise CmdException('encoder "%s" not available' % (encoder))

    if img_ext == '.png':
        _self.set('opaque_background', quiet=quiet)

    # MP4 needs dimensions divisible by 2
    if splitext[1] in ('.mp4', '.mov', '.webm'):
        if width < 1 or height < 1:
            w, h = _self.get_viewport()
            if width > 0:
                height = width * h / w
            elif height > 0:
                width = height * w / h
            else:
                width, height = w, h
        if width % 2:
            width -= 1
        if height % 2:
            height -= 1

    # clean up old files if necessary
    if os.path.exists(filename):
        os.unlink(filename)
    if not os.path.exists(tmp_path):
        os.mkdir(tmp_path)
    elif preserve==0:
        # get rid of existing frames (if they exist)
        for fil in glob.glob(os.path.join(tmp_path,prefix+"*")):
            os.unlink(fil)
    if preserve<0:
        preserve = 0
    if os.path.isdir(tmp_path):
        if first <= 0:
            first = 1
        if last <= 0:
            last = _self.count_frames()
        if last <= 1:
            last = 1
        _self.set("keep_alive")
        _self.mpng(os.path.join(tmp_path,prefix + img_ext),first,last,
                   preserve,mode=mode,modal=-1,quiet=quiet,
                   width=width, height=height)
        # this may run asynchronously
    else:
        ok = 0
    if ok:
        args = ((filename, first, last, preserve,
                                   encoder,tmp_path,prefix,img_ext,
                                   quality,quiet,_self))
        if _self.get_modal_draw():
            t = threading.Thread(target=_encode, args=args)
            t.setDaemon(1)
            t.start()
        else:
            _encode(*args)
    if ok:
        return _self.DEFAULT_SUCCESS
    else:
        return _self.DEFAULT_ERROR
