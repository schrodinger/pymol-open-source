# Demonstration code for a new external gui appliction, which could be
# implemented using any cross-platform Python GUI (Tkinter, wxPython,
#
# To get started:
#
# (1) Copy this file into PyMOL/modules
#   (or create a subdirectory with prefix and rename this file to
#    __init__.py inside it:  "my_ext_gui/__init__.py"  )
#
# (2) Launch PyMOL with options that force it to use this GUI instead:
#
# unix: ./pymol.com -qiF -N my_ext_gui -X 100 -Y 100 -H 400 -W 400 
#
# win: pymolwin.exe +2 -qiF -N my_ext_gui -X 100 -Y 100 -H 400 -W 400 
#  (NOTE: +2 option keeps open the console window for debugging)
#
#
# PyMOL will first import this module, and then call an __init__
# method which should fire off a thread and return, as shown below

import threading

def run(pymol):
   print "\nNow start my custom gui here..."

   # next step might be to import wxPython, set the pymol property,
   # open windows, etc.

   # just for kicks, let's put up some up animated content
   # NOTE that the PyMOL API is thread-safe, so your gui thread
   # can message it asynchronously while the user is able to 
   # interact with the display window

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

def __init__(pymol,poll=0):
   t = threading.Thread(target=run,args=(pymol,))
   t.setDaemon(1)
   t.start()

# note that in order to maintain future compatibility with PyMOL, the
# PyMOL object passed into __init__ should be the exclusive means by
# which you access the API.  Do NOT import pymol directly, since in
# the future, pymol will be an object, and not available as a module.
#
# So store it in your GUI, and then use use as follows: 
#
# pymol.cmd.load(...)
# pymol.cmd.zoom(...)
#
# or perhaps:
# 
# self.cmd = pymol.cmd
# ...
# self.cmd.load(...)
# self.cmd.zoom(...)




