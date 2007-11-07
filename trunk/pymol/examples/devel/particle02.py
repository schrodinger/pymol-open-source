
# demonstration of a particle cloud (non-interacting)

# similar to particle01.py except that each new state is retained in
# memory and used to assemble a movie with the number of states
# provided below

from random import random
from pymol import cmd
from time import sleep

particle_count = 1000

box_size = 500.0

n_states = 200

# constants

half_box = box_size / 2

# create N particle system [x,y,z,r,vx,vy,vz]

particle = [] 
for resi in range(0,particle_count):
    particle.append([resi] + 
                     map(lambda x:(random()-0.5)*box_size/2,[0]*3) + # x,y,z
                    [random()+0.5] + # r
                    map(lambda x:(random()-0.5),[0]*3) # vx,vy,vz
                    )
        
# create cloud object

for part in particle:
    cmd.pseudoatom("cloud",
                   resi = part[0],
                   pos = part[1:4],
                   vdw = part[4])

# draw spheres efficiently
cmd.as("spheres")
cmd.unset("cull_spheres") 

# defer geometry generation until needed

cmd.set("defer_builds",1)

# position the camera

cmd.zoom()
cmd.zoom("center",box_size)

# let there be color

cmd.spectrum()

# this is the main loop

def simulation():
    state = 1 
    import traceback
    try:
        while state < n_states:
            state = state + 1
            for part in particle:
                # simplistic Euler intergration

                # p = p + v
                
                part[1] = (half_box + part[1] + part[5]) % box_size - half_box
                part[2] = (half_box + part[2] + part[6]) % box_size - half_box
                part[3] = (half_box + part[3] + part[7]) % box_size - half_box

                # v = v + pseudo-gravitational acceleration
                
                factor = max(0.1*box_size, 0.1*(part[1]**2+part[2]**2+part[3]**2)**1.5)
                
                part[5] = part[5] - part[1] / factor
                part[6] = part[6] - part[2] / factor
                part[7] = part[7] - part[3] / factor

            # copy initial coordinates to a new state
            
            cmd.create("cloud","cloud",1,state) 

            # update the new state coordinates
            cmd.alter_state(state,"cloud","(x,y,z) = particle[int(resi)][1:4]",space=globals())

            cmd.forward()
            cmd.refresh()

            # don't hog the CPU entirely
            sleep(0.01)
            
        cmd.mplay()
    except:
        traceback.print_exc()

# launch the main loop in a separate thread

import threading

thread = threading.Thread(target=simulation)
thread.setDaemon(1)
thread.start()
