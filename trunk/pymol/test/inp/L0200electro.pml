# -c

/print "BEGIN-LOG"

load lrg/molb.phi,e_map
load lrg/molb.pdb,molb
hide
show surf,molb
ramp_new e_pot,e_map,[-10,0,10]
refresh
isomesh pos,e_map,0.4
color blue,pos
isomesh neg,e_map,-0.4
refresh
color red,neg
color e_pot,molb

ray renderer=2

/print "END-LOG"