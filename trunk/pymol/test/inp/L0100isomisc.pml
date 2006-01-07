# -c

/print "BEGIN-LOG"

load lrg/map.ccp4

isosurf m1,map,1.0
ray renderer=2
isosurf m1,map,0.5
ray renderer=2
isosurf m1,map,500
ray renderer=2
isosurf m1,map,-3000
ray renderer=2

dele m1

/for a in range(1,10): cmd.isosurface("m2","map",1.0+a/10.0,state=a,quiet=0)

ray renderer=2

frame 3
ray renderer=2

frame 8
ray renderer=2

dele all

load lrg/map.xplor

isosurf m1,map,1.0
ray renderer=2
isosurf m1,map,0.5
ray renderer=2
isosurf m1,map,500
ray renderer=2
isosurf m1,map,-3000
ray renderer=2

dele m1

/for a in range(1,10): cmd.isosurface("m2","map",1.0+a/10.0,state=a,quiet=0)

ray renderer=2

frame 3
ray renderer=2

frame 8
ray renderer=2

dele all

load lrg/map.ccp4

load lrg/ref.pdb

isosurf m1,map,1.0,(ref)
ray renderer=2

isosurf m1,map,1.0,A/10:20/
ray renderer=2

isosurf m1,map,1.0,A/50:60/
ray renderer=2

isosurf m1,map,1.0,A/100/
ray renderer=2

isosurf m1,map,1.0,A/100/,2.0
ray renderer=2

isosurf m1,map,1.0,A/100/,4.0
ray renderer=2

isosurf m1,map,1.0,A/53/,carve=2.0
ray renderer=2

isosurf m1,map,1.0,A/53/,carve=100
ray renderer=2

isosurf m1,map,1.0,A/53/,carve=-0.5
ray renderer=2

isosurf m1,map,1.0,A/53/,carve=10
ray renderer=2

isosurf m1,map,1.0,A/53/,carve=0.0
ray renderer=2

dele all

load lrg/map.ccp4

isodot m1,map,1.0
ray renderer=2
isodot m1,map,0.5
ray renderer=2
isodot m1,map,500
ray renderer=2
isodot m1,map,-3000
ray renderer=2

dele m1

/for a in range(1,10): cmd.isodot("m2","map",1.0+a/10.0,state=a,quiet=0)

ray renderer=2

frame 3
ray renderer=2

frame 8
ray renderer=2

dele all

load lrg/map.xplor

isodot m1,map,1.0
ray renderer=2
isodot m1,map,0.5
ray renderer=2
isodot m1,map,500
ray renderer=2
isodot m1,map,-3000
ray renderer=2

dele m1

/for a in range(1,10): cmd.isodot("m2","map",1.0+a/10.0,state=a,quiet=0)

ray renderer=2

frame 3
ray renderer=2

frame 8
ray renderer=2

dele all

load lrg/map.ccp4

load lrg/ref.pdb

isodot m1,map,1.0,(ref)
ray renderer=2

isodot m1,map,1.0,A/10:20/
ray renderer=2

isodot m1,map,1.0,A/50:60/
ray renderer=2

isodot m1,map,1.0,A/100/
ray renderer=2

isodot m1,map,1.0,A/100/,2.0
ray renderer=2

isodot m1,map,1.0,A/100/,4.0
ray renderer=2

isodot m1,map,1.0,A/53/,carve=2.0
ray renderer=2

isodot m1,map,1.0,A/53/,carve=100
ray renderer=2

isodot m1,map,1.0,A/53/,carve=-0.5
ray renderer=2

isodot m1,map,1.0,A/53/,carve=10
ray renderer=2

isodot m1,map,1.0,A/53/,carve=0.0
ray renderer=2



/print "END-LOG"



