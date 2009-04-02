# Contributed functionality which is either
# untested or considered unsupported by
# DeLano Scientific

# contrib.stpng(filename) generates stereo pair of images
# also merges left and right images if Python Image library PIL is available

def stpng(filename="stereo.png"):
    import Image
    filename=str(filename)
    leftfile="%s_left.png"%filename
    rightfile="%s_right.png"%filename
# render left image and save as 'filename_left.png'
    print "Rendering %s"%leftfile
    cmd.turn("y","3")
    cmd.ray()
    cmd.png(leftfile)
# render right image and save as 'filename_right.png'
    print "Rendering %s"%rightfile
    cmd.turn("y","-6")
    cmd.ray()
    cmd.png(rightfile)
# reset original view
    cmd.turn("y","3")
# merge images - can be done with PIL or an external program
############################################################
# Requires Python PIL
# Open left and right images and determine the image size
    leftim=Image.open(leftfile)
    rightim=Image.open(rightfile)
    size=leftim.size
# Set size of new stereo figure and define location of left and right
# stereo images in the new figure
    sizex=size[0]
    sizey=size[1]
    stsizex=2*sizex
    stsizey=sizey
    stsize=(stsizex,stsizey)
    leftbox=(0,0,sizex,sizey)
    rightbox=(sizex,0,stsizex,stsizey)
# Create new 'empty' stereo image
    stereoim=Image.new("RGB",stsize)
# Paste left and right images into new stereo image
    stereoim.paste(leftim,leftbox)
    stereoim.paste(rightim,rightbox)
# Save stereo image - file format is automatically derived from the
# extension of the filename (i.e. png,jpg,tif,pdf)
    stereoim.save(filename)
    print "Saved stereo image:%s"%filename
########################################################
# Merge images using ImageMagick's montage command
#   cmd_montage="montage +frame +label -geometry %s+0+0! -scene 0 %s
# %s %s"%(size,leftfile,rightfile,filename)
#   print "command:%s"%cmd_montage
#   os.system(cmd_montage)
#   print "Saved stereo image:%s"%filename
########################################################
# Remove temporary left and right image files
    os.remove(leftfile)
    os.remove(rightfile)

