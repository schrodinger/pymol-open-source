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

# hide all representations

hide

# show overall protein cartoon

show cartoon

# set everything to cyan

util.cbac

# disable color blending on cartoon

set cartoon_discrete_colors,1  

# don't smooth loops and sheets -- lets use the real coordinates

set cartoon_smooth_loops, 0
set cartoon_flat_sheets, 0

# show valencies

set valence,1
set stick_ball,1

# orient entire system along principal axes

orient

# get context information

context_info=m4x.get_context_info()

# setup the contexts

m4x.setup_contexts(context_info)

# color protein secondary structure (C-alpha atoms only)

_ color white, name ca and ss ''+L and( rep cartoon and not rep sticks)
_ color orange, name ca and ss S and rep cartoon and ((not rep sticks) or (nbr. nbr. nbr. rep cartoon))
_ color magenta, name ca and ss H and rep cartoon and ((not rep sticks) or (nbr. nbr. nbr. rep cartoon))

# set Ctrl Q to Quit 

_ cmd.set_key('CTRL-Q',cmd.quit)
_ feedback enable,python,output
_ print ' '
