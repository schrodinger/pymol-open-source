# -c

/print "BEGIN-LOG"

load lrg/map.ccp4

load lrg/ref.pdb

isomesh m1,map,1.0,(ref)
ray renderer=2

isomesh m1,map,1.0,A/10:20/
ray renderer=2

isomesh m1,map,1.0,A/50:60/
ray renderer=2

isomesh m1,map,1.0,A/100/
ray renderer=2

isomesh m1,map,1.0,A/100/,2.0
ray renderer=2

isomesh m1,map,1.0,A/100/,4.0
ray renderer=2

isomesh m1,map,1.0,A/53/,carve=2.0
ray renderer=2

isomesh m1,map,1.0,A/53/,carve=100
ray renderer=2

isomesh m1,map,1.0,A/53/,carve=-0.5
ray renderer=2

isomesh m1,map,1.0,A/53/,carve=10
ray renderer=2

isomesh m1,map,1.0,A/53/,carve=0.0
ray renderer=2

/print "END-LOG"



