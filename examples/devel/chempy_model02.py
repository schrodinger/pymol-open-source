from pymol import cmd
from random import random
import time

# this shows how you can efficiently update the coordinates
# of an existing model for real-time viewing

# for asychronous execution,  you may want to run this using the "spawn"
# command from inside PyMOL or with the "-l" option from the unix shell

# first we need a model

cmd.load("$PYMOL_PATH/test/dat/pept.pdb","demo")

# let's dress it up a little bit 

cmd.show("sticks","demo")

cmd.show("spheres","resi 10")

cmd.color("yellow","resi 5 and element c")

# turn off some of the chatter about reloading the object...

cmd.feedback("disable","executive","actions")

# now loop, updating the coordinates and reloading the model into
# state 1 of the "demo" object

m = cmd.get_model()
while 1:
   time.sleep(0.05)
   try:
      cmd.set("suspend_updates","1") # only necessary if multithreading...
      for a in m.atom:
         a.coord[0]+=(random()-0.5)*0.1
         a.coord[1]+=(random()-0.5)*0.1
         a.coord[2]+=(random()-0.5)*0.1
      cmd.load_model(m,"demo",1)
   except:
      cmd.set("suspend_updates","0") # only necessary if multithreading...
      traceback.print_exc()
   cmd.refresh()

# Summary: this is portable, safe, but inefficient.  For real-time visualization
# of coordinate changes, there is a way to do this by passing in an opaque
# C data structure...

# Cheers, warren@delanoscientific.com



