
from pymol import cmd
from pymol.cgo import *
from pymol.vfont import plain

cgo = []

axes = [[2.0,0.0,0.0],[0.0,2.0,0.0],[0.0,0.0,2.0]]

pos = [0.0,0.0,0.0]
wire_text(cgo,plain,pos,'Hello World',axes)

pos = [0.0,-3.0,0.0]
cyl_text(cgo,plain,pos,'Hello Universe',0.10,axes=axes)

cmd.set("cgo_line_radius",0.03)
cmd.load_cgo(cgo,'txt')
cmd.zoom("all",2.0)


         
