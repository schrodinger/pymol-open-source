# Add your pmg_tk startup scripts to this directory...

# here are two examples: (copy to new ".py" files and remove comment marks)

# === BEGIN EXAMPLE: myplugin.py ===
#
#from Tkinter import *
#from pymol import cmd
#
#def __init__(self):
#   self.menuBar.addcascademenu('Plugin', 'MyPlugin', 'Sample Plugin',
#                               label='Sample Plugin')
#   
#   self.menuBar.addmenuitem('MyPlugin', 'command',
#                      'Set White Background',
#                      label='Set White Background',
#                      command = lambda : cmd.set("bg_rgb","[1,1,1]"))
#
#   self.menuBar.addmenuitem('MyPlugin', 'command',
#                      'Set Black Background',
#                      label='Set Black Background',
#                      command = lambda : cmd.set("bg_rgb","[0,0,0]"))
# === END EXAMPLE




# === BEGIN EXAMPLE: newmenu.py ===
#
#from Tkinter import *
#from pymol import cmd
#
#def __init__(self):
#   self.menuBar.addmenu('NewMenu', 'Sample Menu')
#
#   self.menuBar.addmenuitem('NewMenu', 'command',
#                      'White background',
#                      label='White Background',
#                     command = lambda : cmd.set("bg_rgb","[1,1,1]"))
#
#   self.menuBar.addmenuitem('NewMenu', 'command',
#                      'Black background',
#                      label='Black Background',
#                     command = lambda : cmd.set("bg_rgb","[0,0,0]"))
#
# === END EXAMPLE ===
