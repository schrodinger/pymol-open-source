
# example script for creation of an image with a slice region

load $PYMOL_PATH/test/dat/1tii.pdb

orient

# must disable depth cue and shadows

unset depth_cue
unset ray_shadows

# this controls the z depth of the slice plane
# (sets it halfway between the clipping planes)

fraction = 0.40
view = cmd.get_view()
near_dist = fraction*(view[16]-view[15])
far_dist = (view[16]-view[15]) - near_dist
cmd.clip("near", -near_dist)

# render opaque background image

as surface
set ray_interior_color, grey80
set opaque_background
set surface_color, white
ray
save image_back.png

cmd.clip("near", near_dist)
cmd.clip("far", far_dist)

# render the foreground image

as cartoon
util.cbc
unset opaque_background
ray
save image_front.png

# now use Photoshop, Gimp, or ImageMagick to combine the images

system composite image_front.png image_back.png image_merged.png
system display image_merged.png
