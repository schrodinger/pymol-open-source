load 1FFK.pdb,1ffk

# create object with only nucleic acid

create nuc = (1ffk and not n;ca,cd)
del 1ffk

# color contiguous chains

edit (c;0 & n;P & i;713)
color blue,(pkchain)

edit (c;0 & n;P & i;2250)
color green,(pkchain)

edit (c;0 & n;P & i;1876)
color magenta,(pkchain)

edit (c;0 & n;P & i;2596)
color violet,(pkchain)

edit (c;0 & n;P & i;1993)
color cyan,(pkchain)

edit (c;0 & n;P & i;2764)
color white,(pkchain)

edit (c;0 & n;P & i;1552)
color salmon,(pkchain)

edit (c;0 & n;P & i;857)
color pink,(pkchain)

edit (c;0 & n;P & i;1149)
color orange,(pkchain)

edit (c;9 & n;P & i;21)
color yellow,(pkchain)

# orient the "hand"

reset
turn y,-60
turn z,180
move y,-15

# show spheres

show spheres

# ray trace

set antialias=1
ray

# write output

png ribosome.png
