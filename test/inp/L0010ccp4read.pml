# -c

# CCP4 map l

/print "BEGIN-LOG"

load lrg/map.ccp4
load lrg/map.ccp4,map2,format=ccp4

set normalize_ccp4_maps = 0
load lrg/map.ccp4,map3

/extent=cmd.get_extent("map")
cmd._dump_floats(extent[0])
cmd._dump_floats(extent[1])

del map
del map2

print cmd.get_names('objects')

enable map2

ray renderer=2

/print "END-LOG"



