# -c

/print "BEGIN-LOG"

# null cases

order does_not_exist
print cmd.get_names()

order *
print cmd.get_names()

order *,location=top
print cmd.get_names()

order *,location=bottom
print cmd.get_names()

order *,location=current
print cmd.get_names()

# singleton cases

load dat/pept.pdb
order pept
print cmd.get_names()

order *
print cmd.get_names()

order *,location=top
print cmd.get_names()

order *,location=bottom
print cmd.get_names()

order *,location=current
print cmd.get_names()

# doubleton

load dat/il2.pdb
print cmd.get_names()

order *
print cmd.get_names()

order il2,location=current
print cmd.get_names()

order il2,location=top
print cmd.get_names()

order il2,location=bottom
print cmd.get_names()

order il2 pept
print cmd.get_names()

load dat/small01.mol
print cmd.get_names()

order small01 pept

print cmd.get_names()

order small01 pept il2

print cmd.get_names()

order *
print cmd.get_names()

create aaa,pept
order *,sort=1
print cmd.get_names()

create zzz,pept
order *,sort=1
print cmd.get_names()

create pept2, pept
create pept1, pept
order pept*
print cmd.get_names()

order pept*,sort=1
print cmd.get_names()

order pept*,location=top
print cmd.get_names()

order pept*,location=bottom
print cmd.get_names()

order aaa,location=bottom
order zzz,location=top
print cmd.get_names()

/print "END-LOG"
