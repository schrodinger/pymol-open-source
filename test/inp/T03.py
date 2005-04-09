# 

# full blown threading stability test, higher enent rate...
#
from pymol import util
import threading
import time
import random
from pymol import cmd
#cmd.feedback("ena","thread","debug")

cmd.rock()
cmd.load("dat/il2.pdb","obj1")
cmd.hide()
cmd.show("ribbon")
cmd.show("car")
util.ss()

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

def sets():
   while 1:
      time.sleep(random.random()*0.15)
      if random.random()>0.5:
         value=1
      else:
         value=0
      cmd.set('cartoon_fancy_helices',str(value))
      if random.random()>0.5:
         value=1
      else:
         value=0
      cmd.set('cartoon_smooth_loop',str(value))
      if random.random()>0.5:
         value=1
      else:
         value=0
      cmd.set('cartoon_round_helices',str(value))
      if random.random()>0.5:
         value=1
      else:
         value=0
      cmd.set('cartoon_smooth_loops',str(value))
      if random.random()>0.5:
         value=1
      else:
         value=0
      cmd.set('cartoon_flat_sheets',str(value))

t = threading.Thread(target=sets)
t.setDaemon(1)
t.start()

def carts():
   while 1:
      resi = int(random.random()*150)
      cmd.cartoon('loop',"(resi %d)"%resi)
      time.sleep(random.random()*0.05)
      resi = int(random.random()*150)
      cmd.cartoon('oval',"(resi %d)"%resi)
      cmd.cartoon('oval',"(resi %d)"%(resi+1))
      time.sleep(random.random()*0.05)
      resi = int(random.random()*150)
      cmd.cartoon('auto',"(resi %d)"%resi)
      cmd.cartoon('auto',"(resi %d)"%(resi+1))      
      time.sleep(random.random()*0.05)
      resi = int(random.random()*150)
      cmd.cartoon('tube',"(resi %d)"%resi)
      cmd.cartoon('tube',"(resi %d)"%(resi+1))
      time.sleep(random.random()*0.05)
      resi = int(random.random()*150)
      cmd.cartoon('rect',"(resi %d)"%resi)
      cmd.cartoon('rect',"(resi %d)"%(resi+1))      
      resi = int(random.random()*150)
      cmd.cartoon('oval',"(resi %d)"%resi)
      time.sleep(random.random()*0.05)
      resi = int(random.random()*150)
      cmd.cartoon('auto',"(resi %d)"%resi)
      time.sleep(random.random()*0.05)
      resi = int(random.random()*150)
      cmd.cartoon('tube',"(resi %d)"%resi)
      time.sleep(random.random()*0.05)
      resi = int(random.random()*150)
      cmd.cartoon('rect',"(resi %d)"%resi)
      time.sleep(random.random()*0.05)
      resi = int(random.random()*150)
      cmd.hide('car',"(resi %d)"%resi)
      time.sleep(random.random()*0.05)
      resi = int(random.random()*150)
      cmd.show('car',"(resi %d)"%resi)
      time.sleep(random.random()*0.05)
      resi = int(random.random()*150)
      cmd.show('car',"(resi %d)"%resi)
      time.sleep(random.random()*0.05)
      resi = int(random.random()*150)
      cmd.show('car',"(resi %d)"%resi)
      time.sleep(random.random()*0.05)
      resi = int(random.random()*150)
      cmd.show('car',"(resi %d)"%resi)
      time.sleep(random.random()*0.05)
      
t = threading.Thread(target=carts)
t.setDaemon(1)
t.start()

