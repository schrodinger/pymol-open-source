# 

# full blown threading stability test, higher enent rate...
#
      
import threading
import time
import random
from pymol import cmd
#cmd.feedback("ena","thread","debug")

cmd.rock()
cmd.load("dat/pept.pdb","obj1")

def bg_rgb():
   while 1:
      time.sleep(random.random()*0.05)
      cmd.set('bg_rgb','%8.3f %8.3f %8.3f'%(
              random.random()/3,
              random.random()/3,
              random.random()/3))

t = threading.Thread(target=bg_rgb)
t.setDaemon(1)
t.start()

def reps():
   while 1:
      time.sleep(random.random()*0.1)
      cmd.show('sticks')
      time.sleep(random.random()*0.2)
      cmd.hide('sticks')

t = threading.Thread(target=reps)
t.setDaemon(1)
t.start()
   
def selector():
   while 1:
      cmd.delete("sel1")
      cmd.select("sel1","(name c)")
      time.sleep(random.random()*0.1)
      cmd.delete("sel2")
      cmd.select("sel2","(i; 10 x; 5)")
      time.sleep(random.random()*0.1)
      cmd.delete("sel3")
      cmd.select("sel3","(obj1)")
      time.sleep(random.random()*0.1)

t = threading.Thread(target=selector)
t.setDaemon(1)
t.start()


def viewport():
   while 1:
      time.sleep(random.random()*0.5)
      cmd.viewport(640,480)
      time.sleep(random.random()*0.5)
      cmd.viewport(800,600)

t = threading.Thread(target=viewport)
t.setDaemon(1)
t.start()

def sets():
   while 1:
      time.sleep(random.random()*0.05)
      if random.random()>0.5:
         ortho=1
      else:
         ortho=0
      cmd.set('ortho',str(ortho))

t = threading.Thread(target=sets)
t.setDaemon(1)
t.start()


def turns():
   while 1:
      time.sleep(random.random()*0.05)
      cmd.turn('x',random.random()*10-5)
      time.sleep(random.random()*0.05)
      cmd.turn('y',random.random()*10-5)
      time.sleep(random.random()*0.05)
      cmd.turn('z',random.random()*10-5)

t = threading.Thread(target=turns)
t.setDaemon(1)
t.start()

