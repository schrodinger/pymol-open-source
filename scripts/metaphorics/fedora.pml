
# quash output

feedback disable, all, everything

# assign secondary structure (if not already assigned)

dss preserve=1

# define selections (may need to be improved)

select protein, not hetatm
select ligand, not ( protein or hoh+wat/ )
select water, ( hoh+wat/ within 5 of ligand )
select contacts, byres ( protein within 4 of ligand )

# hide selection indicators

deselect

# customize settings

set sphere_scale, 0.25
set stick_radius, 0.1
set dot_width, 1.5

# display desired representations

hide
show spheres, ligand
show sticks, ligand or contacts
show cartoon, protein
show nonbonded, water

# draw hydrogen bonds

dist hbond, not ligand, ligand, mode=2
hide labels, hbond

# color protein by secondary structure

color white, protein and ss ''+L
color orange, protein and ss S
color magenta, protein and ss H

# disable color blending on cartoon

set cartoon_discrete_colors,1  

# don't smooth loops -- let's see the real coordinates

set cartoon_smooth_loops, 0
set cartoon_flat_sheets, 0

# color visible atoms by type

util.cbag ligand
util.cbac contacts and not name ca

# now bind some zoom actions to F1-F4

cmd.set_key('F1',lambda :cmd.zoom())
cmd.set_key('F2',lambda :cmd.zoom("ligand",8))
cmd.set_key('F3',lambda :cmd.zoom("ligand"))
cmd.set_key('F4',lambda :cmd.zoom("ligand",-8))

# orient system along principal axes

orient



