# This is an example of firing up PyMOL inside of a subordinate
# process via an "import pymol"
# 
# NOTE: for this to work, PyMOL must be installed in a
# Python-dependent fashion (e.g. pymol-0_98-bin-win32-py23) etc.
#
# WARNING: stability issues have been known to occur with this
# approach, so anticipate problems...take-down is messy.
#
# WARNING: Right now, there is no way for the main process to know
# when PyMOL is actually initialized and ready to go, so we simply
# sleep a second after importing.


import __main__

# note that passing in a "-z" option would keep the window hidden
# until you called pymol.cmd.window("show").

__main__.pymol_argv= "pymol -qxiF  -X 300 -Y 100 -H 400 -W 400".split()
import pymol

# give PyMOL enough time to initialize (we need to find a safe and
# robust alternative to this stupid delay especially since the
# pymol.finish_launching() method now seems to be broken)

import time
time.sleep(1)


# put up some content

if 1:
   pymol.cmd.set("sweep_mode",3)
   pymol.cmd.rock()
   pymol.cmd.turn("x",180)
   pymol.cmd.load("$TUT/1hpv.pdb")
   pymol.preset.pretty("1hpv")
   pymol.cmd.orient()
   pymol.cmd.turn("y",85)
   pymol.cmd.zoom("all",20)
   pymol.cmd.orient("organic & e. N+O",animate=10)
   pymol.cmd.show("sticks","organic")

# play peek-a-boo with the window

if 1:
   time.sleep(5)
   pymol.cmd.window("hide")
   print("Peek-a-boo!")
   time.sleep(1)
   pymol.cmd.window("show")
   time.sleep(5)
   pymol.cmd.window("hide")
   print("Peek-a-boo!")
   time.sleep(1)
   pymol.cmd.window("show")
   time.sleep(5)
   pymol.cmd.window("hide")
   print("Peek-a-boo!")
   time.sleep(1)
   pymol.cmd.window("show")

# now quit 

   print("Quitting...")
   time.sleep(1)
   print("3...")
   time.sleep(1)
   print("2...")
   time.sleep(1)
   print("1...")
   time.sleep(1)
   print("Die!")

# note, we cannot let the main thread terminate without first calling
# pymol.cmd.quit() which will take-down PyMOL

   pymol.cmd.quit()

