# -c

/print "BEGIN-LOG"

load lrg/map.ccp4

isomesh m1,map,1.0
ray renderer=2
isomesh m1,map,0.5
ray renderer=2
isomesh m1,map,500
ray renderer=2
isomesh m1,map,-3000
ray renderer=2

dele m1

/for a in range(1,10): cmd.isomesh("m2","map",1.0+a/10.0,state=a)

ray renderer=2

frame 3
ray renderer=2

frame 8
ray renderer=2

dele all

load lrg/map.xplor

isomesh m1,map,1.0
ray renderer=2
isomesh m1,map,0.5
ray renderer=2
isomesh m1,map,500
ray renderer=2
isomesh m1,map,-3000
ray renderer=2

dele m1

/for a in range(1,10): cmd.isomesh("m2","map",1.0+a/10.0,state=a)

ray renderer=2

frame 3
ray renderer=2

frame 8
ray renderer=2

/print "END-LOG"



