# -c

feedback enable, repsurface, blather
feedback enable, isomesh, blather

/print "BEGIN-LOG"

load dat/embed06.p1m
set ray_default_renderer,2
ray

print cmd.get_names()

/print "END-LOG"
