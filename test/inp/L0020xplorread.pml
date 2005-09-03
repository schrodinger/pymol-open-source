# -c

# XPLOR map l

/print "BEGIN-LOG"

load lrg/map.xplor
load lrg/map.xplor,map2,format=xplor

/extent=cmd.get_extent("map")
cmd._dump_floats(extent[0])
cmd._dump_floats(extent[1])

dele map

print cmd.get_names('objects')

enable map2

ray renderer=2
reinit
f=open("lrg/map.xplor")
content=f.read()
f.close()

cmd.load_raw(content, 'xplor', 'map',quiet=0)
dele all
cmd.load_raw(content, 'xplor', 'map',quiet=1)
isomesh m1, map
ray renderer=2

/print "END-LOG"



