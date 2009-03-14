# -c

# unix-specific externing tests

# Tests the commands:
#
# cd
# ls
# pwd
# system

set ray_default_renderer=2 
# for counting graphics primitives

feedback enable, repsurface, blather
feedback enable, isomesh, blather

/print "BEGIN-LOG"


# bg_color

bg white
bg_color black
bg_color 0
bg_color 1
bg_color blue
bg_color -28374283
bg_color 239283498
bg_color asdfasdf

# clip

dele all
reset
clip
clip near
clip near,-5
clip far,10
clip near,2
clip far,-3
cmd._dump_floats(cmd.get_view(0))

# cartoon

dele all
load dat/1tii.pdb
hide
show cartoon
ray
cartoon rect
ray
cartoon loop
ray
cartoon tube
ray
cartoon oval
ray
cartoon arrow
ray
cartoon dumbbell
ray
cartoon automatic
ray
cartoon loop,ss h
cartoon tube,ss s
ray
cartoon skip,resi 1-100
ray
dele all

# color

dele all
/from pymol.cgo import *
/cgo_list = [ BEGIN, LINES, VERTEX, 0.0, 0.0, 0.0, VERTEX, 1.2, -10.3, 14.0, END ]
/cmd.load_cgo(cgo_list,'cgo')
load dat/pept.pdb
color
color -123454
color 1023912
color blue
color red
color 99
color blue,pept
color green,cgo
color red,(all)
dele all

# cgo tuple input (convenience for users who make mistakes)

cgo_tuple = (SPHERE, 0.0, 0.0, 0.0, 1.0)
cmd.load_cgo(cgo_tuple,'cgo')
ray renderer=2
dele all

# disable & enable
dele all
load dat/pept.pdb
load dat/3al1.pdb
disable
enable pept
ray
disable
enable 3al1
ray
disable
enable
ray
disable
enable all
ray
disable
enable 3al1
enable pept
ray
dele all

# get_view

reset
cmd._dump_floats(cmd.get_view(0))

# show & hide

load dat/pept.pdb
hide
ray
show (name ca)
ray
show (name c)
ray
show (all)
ray
hide (name ca)
ray
hide
ray
show sticks
ray
show spheres,(name ca)
ray
show surface,(not name n)
hide (name o)
ray
hide surf
ray
hide sph
ray
hide sticks
ray
dele all


# label

dele all
load dat/pept.pdb
label (all),resi
label (name ca),"ca"
label (name n),"%s-%s-%s"%(chain,resi,resn)
label (all),''
label (all)
dele all

# move

dele all
reset
cmd._dump_floats(cmd.get_view(0))
move
move x
move y
move z
move x,10
move y,5
move z,2
move z,-1
move x,-10
move y,2
cmd._dump_floats(cmd.get_view(0))
dele all

# orient

dele all
load dat/pept.pdb
orient
orient (none)
orient (name ca)
orient (resi 10)
orient (name ca and i. 10)
dele all

# origin

dele all
load dat/pept.pdb
origin
origin (none)
origin (name ca)
origin (resi 10)
origin (name ca and i. 10)
dele all

# rebuild & refresh

dele all
load dat/pept.pdb
show dots
ray
alter (all),vdw=5.0
refresh
ray
rebuild
ray
dele all

# set_color & recolor

dele all
load dat/pept.pdb
color shoot
set_color shoot,[1.0,0.9,1.0]
color shoot
refresh
set_color shoot,[1.0,0.0,1.0]
recolor
refresh
dele all

# set_view

dele all
reset
/v1 = cmd.get_view(0)
cmd._dump_floats(v1)
load dat/pept.pdb
cmd._dump_floats(cmd.get_view(0))
/cmd.set_view(v1)
cmd._dump_floats(cmd.get_view(0))
dele all

# turn

dele all
reset
cmd._dump_floats(cmd.get_view(0))
turn
turn x
turn y
turn z
turn x,10
turn y,5
turn z,2
turn z,-1
turn x,-10
turn y,2
cmd._dump_floats(cmd.get_view(0))
dele all

# view

dele all
reset
view
view 1,store
view 2,store
move x,10
view 2,st
reset
move y,-10
cmd._dump_floats(cmd.get_view(0))
view 1,recall
cmd._dump_floats(cmd.get_view(0))
view 2,rec
cmd._dump_floats(cmd.get_view(0))

turn x,10
view test,store
turn x,20
view rest,store
view test
cmd._dump_floats(cmd.get_view(0))
view rest
cmd._dump_floats(cmd.get_view(0))
view rest,clear
view rest
cmd._dump_floats(cmd.get_view(0))
view test
cmd._dump_floats(cmd.get_view(0))

# viewport

viewport
viewport 300
viewport 100,100
viewport 300,300
viewport 640,480

# zoom

dele all
load dat/pept.pdb
cmd._dump_floats(cmd.get_view(0))
zoom
zoom resi 10
zoom resi 4,5,6
cmd._dump_floats(cmd.get_view(0))

# coordinate set based zoom, orient, origin

dele all
load dat/1tii.pdb,m1
load dat/il2.pdb,m1
load dat/pept.pdb,m1

zoom m1
cmd._dump_floats(cmd.get_view(0))
zoom m1,state=0
cmd._dump_floats(cmd.get_view(0))
zoom m1,state=1
cmd._dump_floats(cmd.get_view(0))
zoom m1,state=2
cmd._dump_floats(cmd.get_view(0))
zoom m1,state=3
cmd._dump_floats(cmd.get_view(0))
zoom m1,state=4
cmd._dump_floats(cmd.get_view(0))

origin m1
cmd._dump_floats(cmd.get_view(0))
origin m1,state=0
cmd._dump_floats(cmd.get_view(0))
origin m1,state=1
cmd._dump_floats(cmd.get_view(0))
origin m1,state=2
cmd._dump_floats(cmd.get_view(0))
origin m1,state=3
cmd._dump_floats(cmd.get_view(0))
origin m1,state=4
cmd._dump_floats(cmd.get_view(0))

orient m1
cmd._dump_ufloats(cmd.get_view(0))
orient m1,state=0
cmd._dump_ufloats(cmd.get_view(0))
orient m1,state=1
cmd._dump_ufloats(cmd.get_view(0))
orient m1,state=2
cmd._dump_ufloats(cmd.get_view(0))
orient m1,state=3
cmd._dump_ufloats(cmd.get_view(0))
orient m1,state=4
cmd._dump_ufloats(cmd.get_view(0))



# safety
zoom m1,state=-5000 
zoom (none),state=5
zoom (all),state=-1
zoom (all),state=0
cmd._dump_ufloats(cmd.get_view(0))

dele all
load dat/1tii.pdb
cmd._dump_floats(cmd.get_extent("1tii")[0])
cmd._dump_floats(cmd.get_extent("1tii")[1])
cmd._dump_floats(cmd.get_extent("1TII")[0])
cmd._dump_floats(cmd.get_extent("1TII")[1])
cmd._dump_floats(cmd.get_extent("1Tii")[0])
cmd._dump_floats(cmd.get_extent("1Tii")[1])

cmd._dump_floats(cmd.get_extent("1tii////cA")[0])
cmd._dump_floats(cmd.get_extent("1tii////ca")[1])
cmd._dump_floats(cmd.get_extent("1TII////cA")[0])
cmd._dump_floats(cmd.get_extent("1TII////ca")[1])
cmd._dump_floats(cmd.get_extent("1Tii////cA")[0])
cmd._dump_floats(cmd.get_extent("1Tii////ca")[1])

dist TSt,/1tii//A/TYR`3/CD1, /1tii//F/ARG`12/NE
cmd._dump_floats(cmd.get_extent("tst")[0])
cmd._dump_floats(cmd.get_extent("tst")[1])
cmd._dump_floats(cmd.get_extent("TST")[0])
cmd._dump_floats(cmd.get_extent("TST")[1])
cmd._dump_floats(cmd.get_extent("tSt")[0])
cmd._dump_floats(cmd.get_extent("tSt")[1])

/print "END-LOG"



