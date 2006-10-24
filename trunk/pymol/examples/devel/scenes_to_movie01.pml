# SUMMARY
#

# This script demonstrates one way of creating a movie from scenes.
# It assumes that we have three scenes, each running for 10 seconds
# (300 frames apiece) including 2-second transitions.

# 1) Load or create content for three scenes (this could just as easily
#    come from a session file).

load $TUT/1hpv.pdb
util.cbc
turn x,180
orient
as cartoon
scene 001, store

show sticks, organic
orient organic
scene 002, store

hide cartoon
show lines, byres organic expand 5
turn x,45
turn y,45
scene 003, store

# 2) Specify a 30-second movie -- state 1, 900 frames at 30 frames per second.

mset 1 x900

# 3) Program scene matrices as movie views at appopriate frames
#    and also add y-axis rocking between scenes.

scene 001, animate=0
mview store, 1
mview store, 240

turn y,-30
mview store, 70
turn y,60
mview store, 170

scene 002, animate=0
mview store, 300
mview store, 540

turn y,-30
mview store, 370
turn y,60
mview store, 470

scene 003, animate=0
mview store, 600
mview store, 840

turn y,-30
mview store, 670
turn y,60
mview store, 770

# 4) Now interpolate the movie camera.

mview interpolate
mview smooth
mview smooth

# 5) Activate scene content at the appropriate movie frames.
 
mdo 1: scene 001, view=0, quiet=1
mdo 240: scene 002, view=0, quiet=1
mdo 540: scene 003, view=0, quiet=1
mdo 840: scene 001, view=0, quiet=1

# 6) Force frame 1 content to load.

rewind

# 6) And play the movie.

mplay



