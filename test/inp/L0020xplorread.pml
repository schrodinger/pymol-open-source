# -c

# XPLOR map l

/print "BEGIN-LOG"

load lrg/map.xplor
load lrg/map.xplor,map2,format=xplor

/extent=cmd.get_extent("map")
cmd._dump_floats(extent[0])
cmd._dump_floats(extent[1])

del map

print cmd.get_names('objects')

enable map2

ray renderer=2

/print "END-LOG"



