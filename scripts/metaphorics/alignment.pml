from pymol import m4x

# quash output

#feedback disable, all, everything

# assign secondary structure (if not already assigned)

dss preserve=1

# hide selection indicators

deselect

# customize settings

set sphere_scale, 0.25
set stick_radius, 0.15
set dot_width, 1.5

# hide default representations

hide lines
hide nonbonded

# show overall proteins as CA-trace

show ribbon, m4x_aligned
set ribbon_sampling, 1

# hide other stuff

hide dashes
hide labels

# disable color blending on cartoon 

#set cartoon_discrete_colors,1  

# don't smooth loops and sheets -- lets use the real coordinates

#set cartoon_smooth_loops, 0
#set cartoon_flat_sheets, 0

# show valencies

set valence,1
set stick_ball,1

# orient entire system along principal axes

orient m4x_aligned

# get context information

context_info=m4x.get_context_info()

# setup the contexts

m4x.setup_alignment_contexts(context_info)

# set Ctrl Q to Quit 

_ cmd.set_key('CTRL-Q',cmd.quit)
_ feedback enable,python,output
_ print ' '

