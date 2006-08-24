# -c

/print "BEGIN-LOG"

set auto_zoom,0
fragment phe
map_new pot,coulomb,0.5,phe,5

slice_new slice, pot
ramp_new ramp, pot
color ramp,slice

set ray_default_renderer,2

refresh
ray

set slice_track_camera, on

turn x,45
turn z,45
refresh
ray

set slice_track_camera,off

refresh
reset

turn y,45

ray

set slice_height_map 
ray

set slice_grid, 1
ray

set slice_grid, 0.2
ray

ramp_new ramp, pot, color=sludge
ray

ramp_new ramp, pot, color=rainbow
ray

ramp_new ramp, pot, color=rainbow, selection=none
ray

/print "END-LOG"
