# 

# full blown threading stability test, higher enent rate...
#
      
import threading
import time
import random
from pymol import cmd,util
#cmd.feedback("ena","thread","debug")

cmd.load("dat/1tii.pdb","obj1")
cmd.hide()
cmd.set("cull_spheres",0)

def turns():
   while 1:
      time.sleep(random.random()*0.02)
      cmd.turn('x',random.random()*10-5)
      cmd.refresh()
      time.sleep(random.random()*0.02)
      cmd.turn('y',random.random()*10-5)
      time.sleep(random.random()*0.02)
      cmd.turn('z',random.random()*10-5)

t = threading.Thread(target=turns)
t.setDaemon(1)
t.start()

def centers():
   while 1:
      try:
         resi = int(random.random()*150)
         cmd.center("(resi %d)"%resi)
         time.sleep(random.random()*0.30)
      except: 
         print "exception"

t = threading.Thread(target=centers)
t.setDaemon(1)
t.start()

def sets():
   while 1:
      time.sleep(random.random()*0.02)
      if random.random()>0.5:
         value=1
      else:
         value=0
      resi = int(random.random()*150)
      cmd.center("(resi %d)"%resi)
      cmd.set('cartoon_fancy_helices',str(value))

      time.sleep(random.random()*0.02)
      if random.random()>0.5:
         value=1
      else:
         value=0
      cmd.set('cartoon_smooth_loop',str(value))

      time.sleep(random.random()*0.02)
      if random.random()>0.5:
         value=1
      else:
         value=0
      cmd.set('cartoon_round_helices',str(value))

      time.sleep(random.random()*0.02)      
      if random.random()>0.5:
         value=1
      else:
         value=0
      cmd.set('cartoon_smooth_loops',str(value))

      time.sleep(random.random()*0.02)
      if random.random()>0.5:
         value=1
      else:
         value=0
      cmd.set('cartoon_flat_sheets',str(value))

      time.sleep(random.random()*0.02)         
      if random.random()>0.5:
         value=1
      else:
         value=0
      cmd.set('stick_radius',random.random()*0.2+0.1)
      
      time.sleep(random.random()*0.02)      
      if random.random()>0.5:
         value=1
      else:
         value=0
      cmd.set('sphere_scale',random.random()*0.5+0.75)


t = threading.Thread(target=sets)
t.setDaemon(1)
t.start()

def carts():
   while 1:
      try:

         resi = int(random.random()*150)
         cmd.hide('everything',"(resi %d)"%resi)
         cmd.show('sticks',"(resi %d)"%resi)
         time.sleep(random.random()*0.03)

         resi = int(random.random()*150)
         cmd.hide('everything',"(resi %d)"%resi)
         cmd.show('spheres',"(resi %d)"%resi)
         time.sleep(random.random()*0.03)

         resi = int(random.random()*150)
         cmd.hide('everything',"(resi %d)"%resi)
         cmd.show('cartoon',"(resi %d)"%resi)
         time.sleep(random.random()*0.03)

         resi = int(random.random()*150)
         cmd.hide('everything',"(resi %d)"%resi)
         cmd.show('lines',"(resi %d)"%resi)
         time.sleep(random.random()*0.03)

   #      resi = int(random.random()*150)
   #      cmd.show('dots',"(resi %d)"%resi)
   #      time.sleep(random.random()*0.03)

   #      resi = int(random.random()*150)
   #      cmd.hide('dots',"(resi %d)"%resi)
   #      time.sleep(random.random()*0.03)

         resi = int(random.random()*150)
         cmd.hide('everything',"(resi %d)"%resi)
         cmd.show('ribbon',"(resi %d)"%resi)
         time.sleep(random.random()*0.03)

         resi = int(random.random()*150)
         cmd.hide('everything',"(resi %d)"%resi)
         cmd.show('nonbonded',"(resi %d)"%resi)
         time.sleep(random.random()*0.03)

         resi = int(random.random()*150)
         cmd.hide('everything',"(resi %d)"%resi)
         cmd.show('nb_spheres',"(resi %d)"%resi)
         time.sleep(random.random()*0.03)

      except:
         print "exception"

t = threading.Thread(target=carts)
t.setDaemon(1)
t.start()

def colors():
   while 1:
      color = 0
      for a in cmd.index("name ca"):
         cmd.color(str(color),"byres %s`%d"%a,quiet=1)
         color = color+1
         if color > 50:
            color = 0
         time.sleep(0.001)
         
t = threading.Thread(target=colors)
t.setDaemon(1)
t.start()

