# -c

/print "BEGIN-LOG"

from pymol2 import PyMOL

lst = []

for i in range(0,3): lst.append( PyMOL() )

for pi in lst: pi.start()

cmd._dump_floats(cmd.get_view())
for pi in lst: cmd._dump_floats(pi.cmd.get_view())

lst[0].cmd.load("dat/pept.pdb",quiet=0)
lst[1].cmd.load("dat/3al1.pdb",quiet=0)

print cmd.count_atoms()
for pi in lst: print pi.cmd.count_atoms()

cmd._dump_floats(cmd.get_view())
for pi in lst: cmd._dump_floats(pi.cmd.get_view())

cmd._dump_floats(lst[0].cmd.get_view())
lst[0].cmd.view("F1","store")
lst[0].cmd.turn("x",45)
cmd._dump_floats(lst[0].cmd.get_view())
lst[0].cmd.view("F2","store")
cmd._dump_floats(cmd.get_view())

lst[0].cmd.view("F1",animate=0)
cmd._dump_floats(lst[0].cmd.get_view())
cmd._dump_floats(cmd.get_view())
lst[0].cmd.view("F2",animate=0)
cmd._dump_floats(lst[0].cmd.get_view())
cmd._dump_floats(cmd.get_view())

cmd._dump_floats(lst[1].cmd.get_view())
lst[1].cmd.view("F1","store")
lst[1].cmd.turn("y",45)
cmd._dump_floats(lst[1].cmd.get_view())
lst[1].cmd.view("F2","store")

lst[1].cmd.view("F1",animate=0)
cmd._dump_floats(lst[1].cmd.get_view())
cmd._dump_floats(cmd.get_view())
lst[1].cmd.view("F2",animate=0)
cmd._dump_floats(lst[1].cmd.get_view())
cmd._dump_floats(cmd.get_view())

for pi in lst: cmd._dump_floats(pi.cmd.get_view())

for pi in lst: pi.stop()

del lst

/print "END-LOG"
