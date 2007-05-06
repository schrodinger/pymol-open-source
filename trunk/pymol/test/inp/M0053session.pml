# -c

/print "BEGIN-LOG"

from pymol2 import PyMOL

lst = []

for i in range(0,3): lst.append( PyMOL() )

for pi in lst: pi.start()

cmd.load("dat/1tii.pdb",quiet=0)
lst[0].cmd.load("dat/pept.pdb",quiet=0)
lst[1].cmd.load("dat/3al1.pdb",quiet=0)

print cmd.count_atoms()
for pi in lst: print pi.cmd.count_atoms()

cmd._dump_floats(lst[0].cmd.get_view())
lst[0].cmd.ray(renderer=2)
lst[0].cmd.view("F1","store")
lst[0].cmd.turn("x",45)
cmd._dump_floats(lst[0].cmd.get_view())
lst[0].cmd.show_as("cartoon")
lst[0].cmd.ray(renderer=2)
lst[0].cmd.view("F2","store")

cmd._dump_floats(lst[1].cmd.get_view())
lst[1].cmd.ray(renderer=2)
lst[1].cmd.view("F1","store")
lst[1].cmd.turn("y",45)
cmd._dump_floats(lst[1].cmd.get_view())
lst[1].cmd.show_as("cartoon")
lst[1].cmd.ray(renderer=2)
lst[1].cmd.view("F2","store")

cmd._dump_floats(cmd.get_view())
for pi in lst: cmd._dump_floats(pi.cmd.get_view())

lst[0].cmd.hide("everything")
lst[1].cmd.hide("everything")

lst[0].cmd.save("tmp/M0053sess0.pse",quiet=0)
lst[1].cmd.save("tmp/M0053sess1.pse",quiet=0)

lst[0].cmd.load("tmp/M0053sess1.pse",quiet=0)
lst[1].cmd.load("tmp/M0053sess0.pse",quiet=0)

cmd.ray(renderer=2)
for pi in lst: pi.cmd.ray(renderer=2)

print cmd.count_atoms()
or pi in lst: print pi.cmd.count_atoms()

cmd._dump_floats(cmd.get_view())
for pi in lst: cmd._dump_floats(pi.cmd.get_view())

lst[0].cmd.show_as("lines")
lst[1].cmd.show_as("lines")

lst[0].cmd.view("F1",animate=0)
lst[1].cmd.view("F1",animate=0)
cmd._dump_floats(cmd.get_view())
for pi in lst: cmd._dump_floats(pi.cmd.get_view())
cmd.ray(renderer=2)
for pi in lst: pi.cmd.ray(renderer=2)

lst[0].cmd.show_as("cartoon")
lst[1].cmd.show_as("cartoon")

lst[0].cmd.view("F2",animate=0)
lst[1].cmd.view("F2",animate=0)
cmd._dump_floats(cmd.get_view())
for pi in lst: cmd._dump_floats(pi.cmd.get_view())
cmd.ray(renderer=2)
for pi in lst: pi.cmd.ray(renderer=2)


for pi in lst: pi.stop()

del lst

/print "END-LOG"
