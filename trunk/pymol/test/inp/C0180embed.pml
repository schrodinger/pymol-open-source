# -c

/print "BEGIN-LOG"

load dat/embed01.p1m
set ray_default_renderer,2
ray

load dat/embed02.p1m
ray

load dat/embed03.p1m
ray

load dat/embed04.p1m
ray

load dat/embed05.p1m
ray

print cmd.get_state()
print cmd.get_names()

/print "END-LOG"
