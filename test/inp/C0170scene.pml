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

/print "END-LOG"



