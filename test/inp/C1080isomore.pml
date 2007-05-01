# -c

feedback enable, isosurface, blather
feedback enable, isomesh, blather
 
print "BEGIN-LOG"

load dat/pept.pdb,obj
load dat/il2.pdb,obj
load dat/3al1.pdb,obj,state=4

# all states one giant map
map_new tst1,gaussian,1.0,state=-3

# selected states
map_new tst2,gaussian,1.0,state=2
map_new tst2,gaussian,1.0,state=4

# independent states in one extent
map_new tst3,gaussian,1.0,state=-4 

# selected empty state
map_new tst4,gaussian,1.0,state=3

# selected single state
map_new tst5,gaussian,1.0,state=4

# indepenent states in independent extents
map_new tst6,gaussian,1.0

frame 4
isomesh m8,tst2,state=-1

frame 3
isomesh m9,tst1,state=-1,source_state=1

isomesh m7,tst1,state=4,source_state=1
refresh
isomesh m7,tst6,state=3,source_state=1
refresh
isomesh m7,tst5,state=2,source_state=4
refresh
isomesh m7,tst2,state=1,source_state=2
refresh

isomesh m1,tst1
refresh
isomesh m2,tst2
refresh
isomesh m3,tst3
refresh

isomesh m4,tst4
refresh
isomesh m5a,tst5,source_state=1
refresh
isomesh m5b,tst5

isomesh m6,tst6

isosurface s1,tst1
refresh
isosurface s2,tst2
refresh
isosurface s3,tst3
refresh
isosurface s4,tst4
refresh
isosurface s4,tst5
refresh
isosurface s5a,tst5,source_state=1
refresh
isosurface s5b,tst5

isosurface s6,tst6

isosurface s7,tst1,state=4,source_state=1
refresh
isosurface s7,tst6,state=3,source_state=1
refresh
isosurface s7,tst5,state=2,source_state=4
refresh
isosurface s7,tst2,state=1,source_state=2
refresh

frame 4
isosurface s8,tst2,state=-1

frame 3
isosurface s9,tst1,state=-1,source_state=1

print "END-LOG"