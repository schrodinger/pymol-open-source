# -c

# massive test of miscellanous functionality at the command level
# this script isn't designed so much to test results as to make sure
# that PyMOL is able to parser the commands given.
# 
# This ugly test was created in order to provide a means of rapidly 
# testing broad range of features while I muck around with the parser.
# Better focused and more thorough tests will come later.

/print "BEGIN-LOG"


load dat/pept.pdb

# don't change these, they're used later on

select lb=(i;11&n;n)
select rb=(i;11&n;c) 

# trivial stuff

cls

# nearly trivial

print cmd.count_states()
print cmd.count_states("pept")

# do

cmd.do("print 'it be done'")

# display stuff which shouldn't have any effect in command mode

viewport 200,150
stereo on
turn x,90
origin 
origin (i;11)
origin pept
zoom (i;11)
zoom pept
zoom
frame 1
frame 2
move x,2
clip near,2
clip far,2
refresh
reset
meter_reset
mclear
mset 1 x10
mdo 1:turn x,5
mdo 2:turn y,10
rock
rock
forward
backward
rewind
ending
middle
mplay
mstop
color red
color blue,(name n)
flag 3=(resn his)
set_color dude = [ 1.0,0.6,0.3]
color dude,(name o)
show sticks,(all)
hide 
show lines
mmatrix store
turn x,45
mmatrix recall
disable
enable
mmatrix clear
refresh
rebuild
refresh

set line_smooth=1

# labels

label (all),chain
label (all),''


# iterate, alter, etc.

stored.cnt = 0
iterate (all),stored.cnt=stored.cnt+1
print stored.cnt

stored.x = 0
iterate_state 1,(all),stored.x=stored.x+x
print stored.x

alter (all),segi='PEPT'

alter_state 1,(all),x = x + 1.0

# stuff which relies on explicit selections

print cmd.get_dihedral("(i;11&n;n)","(i;11&n;ca)","(i;11&n;cb)","(i;11&n;cg)")

# stuff which relies on mouse button selections

dist 
dist (lb),(rb)
dist tst 
dist tst2 = (lb),(rb)

bond 
unbond 
bond (lb),(rb),2
unbond (lb),(rb)

# copying, and some 2 object stuff

copy cpy=pept
select cpy_sel = (cpy)

alter_state 1,(cpy),x=x+math.cos(y+z)

rms cpy,pept
rms (cpy and name ca),pept
rms (cpy and name ca),(pept &n;ca)
rms_cur cpy,pept
rms_cur (cpy and name ca),pept
rms_cur (cpy and name ca),(pept &n;ca)
fit cpy,pept
fit (cpy),(pept)
fit (cpy and name ca),(pept &n;ca)

pair_fit (cpy&i;11&n;ca),(pept&i;11&n;ca),\
         (cpy&i;5&n;ca),(pept&i;5&n;ca),\
         (cpy&i;1&n;ca),(pept&i;1&n;ca),\
         (cpy&i;8&n;ca),(pept&i;8&n;ca)
update (cpy),(pept)

remove (cpy and name c)

edit (cpy and i;11 & n;ca)
remove_picked
edit (cpy and i;5 & n;ca),(cpy and i;5 and n;cb)

alter (cpy),name=''
rename cpy
rename cpy,1

cycle_valence
cycle_valence
cycle_valence
remove_picked

edit (cpy and i;11 &n;n)
cmd.attach("H",1,1)
create tmp = (cpy&i;5)
fuse (cpy&i;11&n;cb),(tmp&n;cb)

edit (cpy and i;10 & n;n)
h_fill
h_fill

h_add (cpy)

edit (cpy and i;8 &n;o)
cmd.replace("S",3,1)

delete tmp
delete cpy
delete cpy_sel

edit (i;11&n;ca)
unpick

#
protect (pept and n;ca)
protect pept
protect
deprotect (pept and n;ca)
deprotect 

mask (pept and n;ca)
mask pept
mask
unmask (pept and n;ca)
unmask


# stuff which relies on an active editing selection

select lb=(i;11&n;n)
select rb=(i;11&n;c) 
edit (i;11&n;ca)
invert 
invert (lb),(rb)

edit (i;11&n;ca),(i;11&n;cb)
torsion 10
torsion -10
torsion 180
torsion 180

# API only

print cmd.get_model().__class__
print cmd.get_area()
print cmd.get_area()
print cmd.get_area("(name ca)")
print cmd.get_names()
print cmd.get_type('pept')
print cmd.identify("(i;7)")
print cmd.get_extent()

# unsupported features

sort
sort pept

#feedback ena,objectmol,debug
spheroid
spheroid pept

# internal functions

cmd.config_mouse(quiet=0)


/print "END-LOG"

