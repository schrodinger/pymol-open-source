# -c

feedback enable, repsurface, blather
feedback enable, isomesh, blather

/print "BEGIN-LOG"

set auto_zoom,0
fragment phe
map_new t1,coulomb,0.5,phe,5
fragment phe,phe
map_new t2,coulomb,0.5,phe,5
map_new t3,coulomb,0.5,phe,5,state=-3

isomesh m1,t1,5
refresh
isomesh m2,t2,5
refresh
isomesh m3,t3,5
refresh


/print "END-LOG"