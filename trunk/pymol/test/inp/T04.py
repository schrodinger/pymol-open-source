# 

# full blown threading stability test, higher enent rate...
#
      
import threading
import time
import whrandom
from pymol import cmd,util
#cmd.feedback("ena","thread","debug")

cmd.rock()
cmd.load("dat/1tii.pdb","obj1")
cmd.hide()

def turns():
   while 1:
      time.sleep(whrandom.random()*0.01)
      cmd.turn('x',whrandom.random()*10-5)
      time.sleep(whrandom.random()*0.01)
      cmd.turn('y',whrandom.random()*10-5)
      time.sleep(whrandom.random()*0.01)
      cmd.turn('z',whrandom.random()*10-5)

t = threading.Thread(target=turns)
t.setDaemon(1)
t.start()

def centers():
   while 1:
      try:
         resi = int(whrandom.random()*150)
         cmd.center("(resi %d)"%resi)
         time.sleep(whrandom.random()*0.10)
      except:
         pass

t = threading.Thread(target=centers)
t.setDaemon(1)
t.start()

def sets():
   while 1:
      time.sleep(whrandom.random()*0.15)
      if whrandom.random()>0.5:
         value=1
      else:
         value=0
      resi = int(whrandom.random()*150)
      cmd.center("(resi %d)"%resi)
      cmd.set('cartoon_fancy_helices',str(value))
      if whrandom.random()>0.5:
         value=1
      else:
         value=0
      cmd.set('cartoon_smooth_loop',str(value))
      if whrandom.random()>0.5:
         value=1
      else:
         value=0
      cmd.set('cartoon_round_helices',str(value))
      if whrandom.random()>0.5:
         value=1
      else:
         value=0
      cmd.set('cartoon_smooth_loops',str(value))
      if whrandom.random()>0.5:
         value=1
      else:
         value=0
      cmd.set('cartoon_flat_sheets',str(value))
      if whrandom.random()>0.5:
         value=1
      else:
         value=0
      cmd.set('cull_spheres',str(value))
      if whrandom.random()>0.5:
         value=1
      else:
         value=0
      cmd.set('stick_radius',whrandom.random()*0.2+0.1)
      if whrandom.random()>0.5:
         value=1
      else:
         value=0
      cmd.set('sphere_scale',whrandom.random()*0.5+0.75)


t = threading.Thread(target=sets)
t.setDaemon(1)
t.start()

def carts():
   while 1:
      resi = int(whrandom.random()*150)
      cmd.show('sticks',"(resi %d)"%resi)
      time.sleep(whrandom.random()*0.02)

      resi = int(whrandom.random()*150)
      cmd.hide('sticks',"(resi %d)"%resi)
      time.sleep(whrandom.random()*0.02)

      resi = int(whrandom.random()*150)
      cmd.show('spheres',"(resi %d)"%resi)
      time.sleep(whrandom.random()*0.02)

      resi = int(whrandom.random()*150)
      cmd.hide('spheres',"(resi %d)"%resi)
      time.sleep(whrandom.random()*0.02)

      resi = int(whrandom.random()*150)
      cmd.show('cartoon',"(resi %d)"%resi)
      time.sleep(whrandom.random()*0.02)

      resi = int(whrandom.random()*150)
      cmd.hide('cartoon',"(resi %d)"%resi)
      time.sleep(whrandom.random()*0.02)

      resi = int(whrandom.random()*150)
      cmd.show('lines',"(resi %d)"%resi)
      time.sleep(whrandom.random()*0.02)
      
      resi = int(whrandom.random()*150)
      cmd.hide('lines',"(resi %d)"%resi)
      time.sleep(whrandom.random()*0.02)

      resi = int(whrandom.random()*150)
      cmd.show('dots',"(resi %d)"%resi)
      time.sleep(whrandom.random()*0.02)
      
      resi = int(whrandom.random()*150)
      cmd.hide('dots',"(resi %d)"%resi)
      time.sleep(whrandom.random()*0.02)

      resi = int(whrandom.random()*150)
      cmd.show('ribbon',"(resi %d)"%resi)
      time.sleep(whrandom.random()*0.02)

      resi = int(whrandom.random()*150)
      cmd.hide('ribbon',"(resi %d)"%resi)
      time.sleep(whrandom.random()*0.02)

t = threading.Thread(target=carts)
t.setDaemon(1)
t.start()

