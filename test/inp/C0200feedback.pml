# -c

/print "BEGIN-LOG"

feedback 
feedback ?
feedback enable
feedback enable, selector
feedback enable, selector, everything

feedback enable, repcartoon extrude, everything
load dat/pept.pdb
show cartoon
refresh
feedback disable, repcartoon extrude, everything
rebuild
refresh

/print "END-LOG"
