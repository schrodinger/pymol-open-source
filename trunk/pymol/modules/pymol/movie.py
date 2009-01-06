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
import os
import glob
import threading
import time

def sweep(pause=0,cycles=1,_self=cmd):
    pause = int(pause)
    cycles = int(cycles)
    n_state = _self.count_states("all")
    if pause>0:
        pass_string = "1 x%d 1 -%d %d x%d %d -1"%(pause, n_state, n_state, pause, n_state)
    else:
        pass_string = "1 -%d %d -1"%(n_state, n_state)
    movie_list = [ pass_string ] * cycles
    movie_string = string.join(movie_list," ")
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
    movie_string = string.join(movie_list," ")
    _self.mset(movie_string)
                
def load(*args,**kw):
    _self = kw.get('_self',cmd)
    nam = "mov"
    if len(args)>1:
        nam = args[1]
    fils = glob.glob(args[0])
    fils.sort()
    if not len(fils):
        print "Error: no matching files"
    else:
        for a in fils:
            apply(_self.load,(a,nam),kw)
#         _self.load(a,nam)

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
            print leftover[1]
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
    print (" tdroll: defined rotations for", frpos - 1,
             "frames, starting at frame %d"%first)

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
    frames_per_sec = float(_self.get('movie_fps'))
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
            _self.mview('store',frame)
            frame = frame + 1
    
# new matrix-based camera interpolation routines

def add_roll(duration=12.0,loop=1,axis='y',_self=cmd):
    cmd = _self
    start_frame = cmd.get_movie_length()+1
    duration = float(duration)
    fps = float(cmd.get('movie_fps'))
    n_frame = int(round(fps * duration))
    if n_frame > 0:
        cmd.madd("1 x%d"%n_frame)
        cmd.frame(start_frame)
        cmd.mview("store")
        cmd.frame(start_frame+n_frame/3)
        cmd.turn(axis,120)
        cmd.mview("store")
        cmd.frame(start_frame+(2*n_frame)/3)
        cmd.turn(axis,120)
        cmd.mview("store")
        if loop:
            if (start_frame == 1):
                cmd.mview("interpolate",power=1,wrap=1)
                cmd.frame(start_frame+n_frame-1)
                cmd.mview("store")
                cmd.turn(axis,120)
            else:
                adjustment = 360.0/n_frame
                cmd.frame(start_frame+n_frame-1)
                cmd.turn(axis,120 - adjustment)
                cmd.mview("store")
                cmd.mview("interpolate",power=1)
                cmd.turn(axis,adjustment) 
        else:
            cmd.frame(start_frame+n_frame-1)
            cmd.turn(axis,120)
            cmd.mview("store")
            cmd.mview("interpolate",power=1)
        cmd.frame(start_frame)
        
def add_rock(duration=8.0,angle=30.0,loop=1,axis='y',_self=cmd):
    cmd = _self
    start_frame = cmd.get_movie_length()+1
    duration = float(duration)
    angle=float(angle)
    loop=int(loop)
    fps = float(cmd.get('movie_fps'))
    n_frame = int(round(fps * duration))
    if n_frame > 0:
        cmd.madd("1 x%d"%n_frame)
        cmd.frame(start_frame+n_frame/4)
        cmd.turn(axis,angle/2.0)
        cmd.mview("store")
        cmd.frame(start_frame+(3*n_frame)/4)
        cmd.turn(axis,-angle)
        cmd.mview("store")
        cmd.turn(axis,angle/2.0)
        if loop and (start_frame == 1):
            cmd.mview("interpolate",power=-1,wrap=1)
        else:
            cmd.mview("interpolate",power=-1)
        cmd.frame(start_frame)

def add_nutate(duration=8.0, angle=30.0, spiral=0, loop=1, 
               offset=0, phase=0, shift=math.pi/2.0,
               _self=cmd):
    cmd = _self
    start_frame = cmd.get_movie_length()+1
    duration = float(duration)
    angle = float(angle)
    spiral = int(spiral)
    loop = int(loop)
    fps = float(cmd.get('movie_fps'))
    n_frame = int(round(fps * duration))
    if n_frame > 0:
        cmd.madd("1 x%d"%n_frame)
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
            cmd.frame(start_frame+index)
            cmd.turn('x',x_rot)
            cmd.turn('y',y_rot)
            cmd.mview('store')
            cmd.turn('y',-y_rot)
            cmd.turn('x',-x_rot)
    
def _rock_y(mode,first,last,period,pause,_self=cmd):
    cmd = _self
    n_frame = last - first + 1
    axis = 'y'
    angle = 10
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
                cmd.frame(frame)
                cmd.turn(axis,angle/2.0)
                cmd.mview("store")
                direction = -1
            else:
                cmd.frame(frame)
                cmd.turn(axis, direction * angle)
                cmd.mview("store")
                direction = -direction

        cmd.frame(last)
        cmd.turn(axis,direction * angle/2.0)
        cmd.mview("store")

        cmd.mview("interpolate",first,last,power=-1)

def _nutate_sub(start_frame, stop_frame, angle=15.0, spiral=0, loop=1, 
                offset=0, phase=0, shift=math.pi/2.0, _self=cmd):
    cmd = _self
    angle = float(angle)
    spiral = int(spiral)
    loop = int(loop)
    fps = float(cmd.get('movie_fps'))
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
            cmd.frame(start_frame+index)
            cmd.turn('x',x_rot)
            cmd.turn('y',y_rot)
            cmd.mview('store')
            cmd.turn('y',-y_rot)
            cmd.turn('x',-x_rot)
    
def _nutate(mode,first,last,period,pause,_self=cmd):
    cmd = _self
    n_frame = last - first + 1
    axis = 'y'
    angle = 10
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
        _nutate_sub(frame[0], frame[1], spiral=spiral, _self=cmd)
        spiral = 0
        
def add_scenes(names=None, pause=8.0, cut=0.0, loop=1,
               rock=-1, period=8.0, animate=-1, _self=cmd):
    cmd = _self
    animate = float(animate)
    pause = float(pause)
    period = float(period)
    if period>pause:
        period = pause
    rock = int(rock)
    if animate<0:
        animate = float(cmd.get("scene_animation_duration"))
    if names == None:
        names = cmd.get_scene_list()
    elif cmd.is_str(names):
        names = cmd.safe_alpha_list_eval(names)
    n_scene = len(names)
    start_frame = cmd.get_movie_length()+1
    duration = n_scene*(pause+animate)
    fps = float(cmd.get('movie_fps'))
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
        cmd.madd("1 x%d"%act_n_frame)
        cnt = 0
        for scene in names:
            frame = start_frame+int((cnt*n_frame)/n_scene)
            cmd.frame(frame)
            cmd.mview("store",scene=scene)
            if rock:
                cmd.mview("interpolate",cut=cut,wrap=0)
                sweep_first = frame + 1
                sweep_last = start_frame + int(pause*fps+(cnt*n_frame)/n_scene) - 1
                if sweep_last > act_n_frame:
                    sweep_last = act_n_frame
                n_sweep_frame = sweep_last - sweep_first + 1
                if n_sweep_frame > 0:
                    if sweep_mode<3: # y-axis rock
                        _rock_y(sweep_mode, sweep_first, sweep_last,
                                period, pause, _self=_self)
                    elif sweep_mode == 3:
                        _nutate(sweep_mode, sweep_first, sweep_last,
                                period, pause, _self=_self)
            frame = start_frame+int(pause*fps+(cnt*n_frame)/n_scene)
            if frame <= act_n_frame:
                if sweep_mode!=3:
                    cmd.frame(frame)
                    cmd.mview("store",scene=scene)
            cnt = cnt + 1
        cmd.mview("interpolate",cut=cut,wrap=loop)
        if rock:
            cmd.mview("smooth")
        cmd.frame(start_frame)

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
            print " produce: %d bytes written..."%size
        else:
            tries = tries - 1
            if tries < 0:
                break
        time.sleep(2)
        if done_event.isSet():
            break
            
def _encode(filename,mode,first,last,preserve,
            encoder,tmp_path,prefix,quality,quiet,_self=cmd):
    import os
    tries = 10
    while 1: # loop until all of the files have been created...
        done = 1
        # check for the required output files
        for index in range(first,last+1):
            path = os.path.join(tmp_path,prefix+"%04d.ppm"%index)
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
    if done and ok and (encoder == 'mpeg_encode'):
        try:
            from freemol import mpeg_encode
        except:
            ok = 0
            print "produce-error: Unable to import module freemol.mpeg_encode."
        if ok:
            if not mpeg_encode.validate():
                ok = 0
                print "produce-error: Unable to validate freemol.mpeg_encode."
        if not ok:
            print "produce-error: Unable to create mpeg file."            
        else:
            mpeg_quality = 1+int(((100-quality)*29)/100) # 1 to 30
            input = mpeg_encode.input(filename,tmp_path,
                                      prefix,first,last,mpeg_quality);
            if not quiet:
                print " produce: creating '%s' (in background)..."%(filename)

            done_event = None
            if not quiet:
                done_event = threading.Event()
                t = threading.Thread(target=_watch,
                                     args=(filename,done_event))
                t.setDaemon(1)
                t.start()

            try:
                result = mpeg_encode.run(input)
            finally:
                if done_event != None:
                    done_event.set()
            if not quiet:
                if not os.path.exists(filename):
                    if result != None:
                        print input, result[0], result[1]
                    print " produce: compression failed"
                else:
                    print " produce: finished."
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


def produce(filename, mode='', first=0, last=0, preserve=0,
            encoder='mpeg_encode', quality=60, quiet=1, _self=cmd):
    prefix = _prefix

    if _self.is_string(mode):
        if mode == '':
            if int(cmd.get_setting_legacy('ray_trace_frames')):
                mode = 'ray'
            else:
                mode = 'draw'
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
        quality = 0
    if quality>100:
        quality = 100
    ok = 1
    splitext = os.path.splitext(filename)
    tmp_path = splitext[0]+".tmp"
    if splitext[1]=='':
        splitext=(splitext[0],'.mpg')
    filename = splitext[0]+splitext[1]
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
        _self.mpng(os.path.join(tmp_path,prefix+".ppm"),first,last,
                   preserve,mode=mode,modal=1,quiet=quiet) 
        # this may run asynchronously
    else:
        ok = 0
    if ok:
        t = threading.Thread(target=_encode,
                             args=(filename,mode,first,last,preserve,
                                   encoder,tmp_path,prefix,quality,quiet,_self))
        t.setDaemon(1)
        t.start()
    if ok:
        return _self.DEFAULT_SUCCESS
    else:
        return _self.DEFAULT_ERROR
    
