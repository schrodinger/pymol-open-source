# 

# full blown threading stability test
#
# The user should make it a point of interacting with the GUI
# while this test is running, including using menus in both the
# external and internal guis, issuing commands, and of course
# terminating the script when it is done.
      
import threading
import time
import whrandom
from pymol import cmd

cmd.rock()

def load_save1():
   while 1:
      time.sleep(whrandom.random())
      cmd.delete("obj1")
      cmd.load("dat/pept.pdb","obj1")
      time.sleep(whrandom.random())
      cmd.save("tmp/T01a.pdb","obj1")
      
t = threading.Thread(target=load_save1)
t.setDaemon(1)
t.start()
   
def bg_rgb():
   while 1:
      time.sleep(whrandom.random()*5)
      cmd.set('bg_rgb','%8.3f %8.3f %8.3f'%(
              whrandom.random()/3,
              whrandom.random()/3,
              whrandom.random()/3))

t = threading.Thread(target=bg_rgb)
t.setDaemon(1)
t.start()

def reps():
   while 1:
      time.sleep(whrandom.random()*5)
      cmd.show('sticks')
      time.sleep(whrandom.random()*5)
      cmd.show('surface')
      time.sleep(whrandom.random()*5)
      cmd.show('mesh')

t = threading.Thread(target=reps)
t.setDaemon(1)
t.start()
   
def load_save2():
   while 1:
      time.sleep(whrandom.random())
      cmd.delete("obj2")
      cmd.load("dat/water.pdb","obj2")
      time.sleep(whrandom.random())
      cmd.save("tmp/T01b.pdb","obj2")

t = threading.Thread(target=load_save2)
t.setDaemon(1)
t.start()

def selector():
   while 1:
      cmd.delete("sel1")
      cmd.select("sel1","(name c)")
      time.sleep(whrandom.random())
      cmd.delete("sel2")
      cmd.select("sel2","(i; 10 x; 5)")
      time.sleep(whrandom.random())
      cmd.delete("sel3")
      cmd.select("sel3","(obj1)")
      time.sleep(whrandom.random())

t = threading.Thread(target=selector)
t.setDaemon(1)
t.start()


def viewport():
   while 1:
      time.sleep(whrandom.random()*5)
      cmd.viewport(640,480)
      time.sleep(whrandom.random()*5)
      cmd.viewport(800,600)

t = threading.Thread(target=viewport)
t.setDaemon(1)
t.start()

