# -c

/print "BEGIN-LOG"

load dat/pept.pkl
print cmd.get_names()

delete pept
print cmd.get_names()

load dat/pept.pkl
delete pept

load dat/pept.pkl
delete pept

load dat/il2.pdb
load dat/pept.pkl

print cmd.get_names()

delete il2
delete pept

print cmd.get_names()

load dat/il2.pdb
load dat/pept.pkl
delete il2
delete pept


load dat/il2.pdb
load dat/pept.pkl

delete il2
delete pept

print cmd.get_names()

/print "END-LOG"
