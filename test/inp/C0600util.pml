# -c

/print "BEGIN-LOG"

load dat/il2.pdb
load dat/pept.pdb
load dat/3al1.pdb

util.ss
util.cbac
util.cbay
util.cbaw
util.cbas
util.cbc

hide
show car
zoom
set ray_default_renderer=2
ray

util.rainbow("pept")

/print "END-LOG"
