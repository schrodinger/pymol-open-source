
# Python script extend_demo01.py

# A recipe for creating a new PyMOL command which calls your own PyMOL code

from pymol import cmd

# step 1: define your Python function

def my_color(color="red",sele="name ca"):
    print "I am coloring "+color+ " all atoms in selection: "+sele
    cmd.color(color, sele)

# step 2: add you command to the PyMOL command language

cmd.extend("my_color", my_color)

# step 3: run this script from within PymOL

#PyMOL> run extend_demo01.py

# step 4: try it out!

#PyMOL> my_color
#PyMOL> my_color yellow
#PyMOL> my_color blue, elem c

