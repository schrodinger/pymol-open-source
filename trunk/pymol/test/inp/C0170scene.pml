# -c

# for counting what is displayed
set ray_default_renderer=2

/print "BEGIN-LOG"

scene -10
scene 9999

cmd._dump_floats(cmd.get_view(0))
scene 1,store

load dat/pept.pdb
cmd._dump_floats(cmd.get_view(0))
scene 2,store

scene 1,recall
ray

scene 2
ray
dist tst,5/ca,8/ca
turn x,5
cmd._dump_floats(cmd.get_view(0))
scene 3,store

ray
hide dashes,tst
turn y,5
cmd._dump_floats(cmd.get_view(0))
scene 4,store
ray

scene 3
disable tst
turn z,20
scene 5,store

scene 1
cmd._dump_floats(cmd.get_view(0))
ray

scene 2
cmd._dump_floats(cmd.get_view(0))
ray

scene 3
cmd._dump_floats(cmd.get_view(0))
ray

scene 4
cmd._dump_floats(cmd.get_view(0))
ray

scene 5
cmd._dump_floats(cmd.get_view(0))
ray

scene 0
scene 6
scene 1000

scene 4,clear
scene 4

scene *
scene *,clear
scene *

scene 1

hide
show lines
color red
scene 1, store
count_atoms color red
count_atoms color blue
ray

hide
show sticks
color blue
scene 2, store

count_atoms color red
count_atoms color blue
ray

scene 1
count_atoms color red
count_atoms color blue
ray

scene 2
count_atoms color red
count_atoms color blue
ray

scene 2,rename,new_key=3

scene 1
count_atoms rep lines
count_atoms color red
count_atoms rep sticks
count_atoms color blue
ray

scene 3
count_atoms rep lines
count_atoms color red
count_atoms rep sticks
count_atoms color blue
ray

scene 2
count_atoms rep lines
count_atoms color red
count_atoms rep sticks
count_atoms color blue
ray

scene 3
hide
color green
show spheres
scene auto, update
count_atoms rep spheres
count_atoms color green
count_atoms rep lines
count_atoms color red
count_atoms rep sticks
count_atoms color blue

scene 1
count_atoms rep spheres
count_atoms color green
count_atoms rep lines
count_atoms color red
count_atoms rep sticks
count_atoms color blue

scene 3
count_atoms rep spheres
count_atoms color green
count_atoms rep lines
count_atoms color red
count_atoms rep sticks
count_atoms color blue


/print "END-LOG"



