
from pymol.wizard import Wizard
from pymol import cmd
import pymol
import copy

default_map = [ '', '', '']
default_level = [ 1.0, 3.0, 1.0]
default_radius = 8.0
default_track = 0

class Density(Wizard):

   def __init__(self):

      Wizard.__init__(self)
      
      # mode selection subsystem
      
      self.radius = default_radius
      self.map = copy.deepcopy(default_map)
      self.level = copy.deepcopy(default_level)
      self.track = copy.deepcopy(default_track)
      self.avail_maps = []
      
      self.menu['radius'] = [
                             [1, '4.0 A Radius','cmd.get_wizard().set_radius(4)'],
                             [1, '5.0 A Radius','cmd.get_wizard().set_radius(5)'],
                             [1, '6.0 A Radius','cmd.get_wizard().set_radius(6)'],
                             [1, '8.0 A Radius','cmd.get_wizard().set_radius(8)'],
                             [1, '10.0 A Radius','cmd.get_wizard().set_radius(10)'],
                             [1, '15.0 A Radius','cmd.get_wizard().set_radius(15)'],
                             [1, '20.0 A Radius','cmd.get_wizard().set_radius(20)'],
                             [1, '50.0 A Radius','cmd.get_wizard().set_radius(50)'],
                             ]
                              
      self.menu['map0'] = []
      self.menu['map1'] = []
      self.menu['map2'] = []

      level_menu = lambda x:[ [1, '1.0 sigma','cmd.get_wizard().set_level(%d,1.0)'%x],
                              [1, '1.5 sigma','cmd.get_wizard().set_level(%d,1.5)'%x],
                              [1, '2.0 sigma','cmd.get_wizard().set_level(%d,2.0)'%x],
                              [1, '3.0 sigma','cmd.get_wizard().set_level(%d,3.0)'%x],
                              [1, '5.0 sigma','cmd.get_wizard().set_level(%d,5.0)'%x]]
      
      self.menu['level0'] = level_menu(0)
      self.menu['level1'] = level_menu(1)
      self.menu['level2'] = level_menu(2)
      
      self.menu['track'] = [
      [ 1, "Track & Zoom", 'cmd.get_wizard().set_track(0)'],
      [ 1, "Track & Set Origin", 'cmd.get_wizard().set_track(1)'],
      [ 1, "Track Off", 'cmd.get_wizard().set_track(2)'],
      ]

      if (self.map[0] == '') and (self.map[1] == '') and (self.map[2]==''):
         for a in cmd.get_names(): # automatically load first map we find
            if cmd.get_type(a)=='object:map':
               self.map[0]=a
               break
      self.update_map_menus()
      
   def update_map_menus(self):

      self.avail_maps = []
      
      for a in cmd.get_names('objects'):
         if cmd.get_type(a)=='object:map':
            self.avail_maps.append(a)

      c = 0
      for a in self.map:
         map_kee = 'map'+str(c)
         level_kee = 'level'+str(c)
         self.menu[map_kee] = [[2,'Select Map','']]
         for a in self.avail_maps:
            self.menu[map_kee].append([ 1,a,'cmd.get_wizard().set_map(%d,"%s")'%(c,a) ])
         self.menu[map_kee].append([1,'(none)','cmd.get_wizard().set_map(%d,"")'%c])
         c = c + 1

   def set_track(self,track):
      self.track = track
      cmd.refresh_wizard()
      
   def set_level(self,map,level):
      self.level[map] = level
      self.update_maps()
      cmd.refresh_wizard()

   def set_map(self,map,map_name):
      self.map[map] = map_name
      cmd.refresh_wizard()

   def set_radius(self,radius):
      self.radius = radius
      self.update_maps()
      cmd.refresh_wizard()

   def update_maps(self):
      if 'pk1' in cmd.get_names('selections'):
         save = cmd.get_setting_text('auto_zoom')
         save = cmd.set('auto_zoom',0,quiet=1)                     
         c = 0
         for a in self.map:
            oname = 'w'+str(c+1)+'_'+a
            if oname not in cmd.get_names():
               color = 1
            else:
               color = 0
            if len(a):
               if cmd.get_type(a)=='object:map':
                  cmd.isomesh(oname,a,self.level[c],
                              "(pk1)",self.radius,state=1)
               if color:
                  if c == 0:
                     cmd.color('blue',oname)
                  elif c == 1:
                     cmd.color('white',oname)
                  else:
                     cmd.color('magenta',oname)
            c = c + 1
         save = cmd.set('auto_zoom',save,quiet=1)            
         if self.track<1:
            cmd.zoom("(pk1)",self.radius)
            cmd.clip("slab",self.radius-1)
         elif self.track<2:
            cmd.origin("(pk1)")
      cmd.refresh_wizard()      
# generic set routines

   def zoom(self):
      if 'pk1' in cmd.get_names('selections'):
         cmd.zoom("(pk1)",self.radius)
         cmd.clip("slab",self.radius-1)
      else:
         c = 0
         for a in self.map:
            oname = 'w'+str(c+1)+'_'+a
            if len(a):
               if a in cmd.get_names('objects'):
                  cmd.zoom(oname)
                  cmd.clip("slab",self.radius-1)
            c = c + 1
            
   def get_panel(self):
      self.update_map_menus()
      return [
         [ 1, 'Density Map Wizard',''],
         [ 2, 'Update Maps' , 'cmd.get_wizard().update_maps()'],
         [ 2, 'Zoom' , 'cmd.get_wizard().zoom()'],         
         [ 3, "Radius: %3.1f A"%self.radius,'radius'],
         [ 3, "Map 1: "+self.map[0],'map0'],
         [ 3, "       @ "+str(self.level[0])+" sigma",'level0'],
         [ 3, "Map 2: "+self.map[1],'map1'],
         [ 3, "       @ "+str(self.level[1])+" sigma",'level1'],
         [ 3, "Map 3: "+self.map[2],'map2'],
         [ 3, "       @ "+str(self.level[2])+" sigma",'level2'],
         [ 3, self.menu['track'][self.track][1], 'track' ],
         [ 2, 'Done','cmd.set_wizard()'],
         ]

   def cleanup(self):
      global default_radius,default_map,default_level
      default_radius = self.radius
      default_map = self.map
      default_level = self.level
      default_track = self.track
      self.clear()
      
   def clear(self):
      pass

   def do_pick(self,bondFlag):
      global dist_count
      if not bondFlag:
         if self.track<2:
            self.update_maps()
         
