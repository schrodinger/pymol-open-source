# load pdb file which has a CRYST1 record

load 1DN2.pdb,1dn2

# create neighbor with contact near a certain atom

symexp s,1dn2,(c; F and i;12 and n;N),4.0

# hide everything outside region of interest

hide (!(byres ((1dn2 and c; F and i;12 and n;N) x;12)))

# color-by-atom-cyan on 1dn2 to distinguish it from symmetry-related copy

util.cbac 1dn2

# restore view saved using from get_view in a previous session

set_view (\
    -0.042664494,   -0.003409438,    0.999083281,\
    -0.617746830,    0.786019623,   -0.023697896,\
    -0.785219014,   -0.618191123,   -0.035642438,\
    -1.056621075,   -1.729980350,  -39.743129730,\
    29.859001160,   38.402999878,   19.087999344,\
    27.443143845,   46.443149567,    1.000000000 )                              

# show the single hydrogen bond in the interface

dist ( s01000000 & segi '' & chain B & resn TYR & resi 296 & name O ),\
     ( 1dn2 & segi '' & chain F & resn CYS & resi 12 & name N )

# image

set antialias=1
ray
png contact.png
