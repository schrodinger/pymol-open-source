# 

# full blown threading stability test, higher enent rate...
#
      
import threading
import time
import whrandom
from pymol import cmd
#cmd.feedback("ena","thread","debug")

cmd.rock()
cmd.load("dat/pept.pdb","obj1")

def bg_rgb():
   while 1:
      time.sleep(whrandom.random()*0.05)
      cmd.set('bg_rgb','%8.3f %8.3f %8.3f'%(
              whrandom.random()/3,
              whrandom.random()/3,
              whrandom.random()/3))

t = threading.Thread(target=bg_rgb)
t.setDaemon(1)
t.start()

def reps():
   while 1:
      time.sleep(whrandom.random()*0.1)
      cmd.show('sticks')
      time.sleep(whrandom.random()*0.2)
      cmd.hide('sticks')

t = threading.Thread(target=reps)
t.setDaemon(1)
t.start()
   
def selector():
   while 1:
      cmd.delete("sel1")
      cmd.select("sel1","(name c)")
      time.sleep(whrandom.random()*0.1)
      cmd.delete("sel2")
      cmd.select("sel2","(i; 10 x; 5)")
      time.sleep(whrandom.random()*0.1)
      cmd.delete("sel3")
      cmd.select("sel3","(obj1)")
      time.sleep(whrandom.random()*0.1)

t = threading.Thread(target=selector)
t.setDaemon(1)
t.start()


def viewport():
   while 1:
      time.sleep(whrandom.random()*0.5)
      cmd.viewport(640,480)
      time.sleep(whrandom.random()*0.5)
      cmd.viewport(800,600)

t = threading.Thread(target=viewport)
t.setDaemon(1)
t.start()

def sets():
   while 1:
      time.sleep(whrandom.random()*0.05)
      if whrandom.random()>0.5:
         ortho=1
      else:
         ortho=0
      cmd.set('ortho',str(ortho))

t = threading.Thread(target=sets)
t.setDaemon(1)
t.start()


def turns():
   while 1:
      time.sleep(whrandom.random()*0.05)
      cmd.turn('x',whrandom.random()*10-5)
      time.sleep(whrandom.random()*0.05)
      cmd.turn('y',whrandom.random()*10-5)
      time.sleep(whrandom.random()*0.05)
      cmd.turn('z',whrandom.random()*10-5)

t = threading.Thread(target=turns)
t.setDaemon(1)
t.start()

