# -c

/print "BEGIN-LOG"

load lrg/molb.phi,e_map
load lrg/molb.pdb,molb
hide
show surf,molb
ramp_new e_pot,e_map,[-10,0,10]
isomesh pos,e_map,0.4
color blue,pos
isomesh neg,e_map,-0.4
color red,neg
color e_pot,molb

/print "END-LOG"