#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information. 
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-* 
#-* 
#-*
#Z* -------------------------------------------------------------------

# pmm.py 
# This section contains the menu contents and the associated commands
#

def mol_show(s):
   return [[ 2, 'Show:'     ,''                             ],
           [ 1, 'lines'     ,'pm.show("lines"     ,"'+s+'")'],
           [ 1, 'sticks'    ,'pm.show("sticks"    ,"'+s+'")'],
           [ 1, 'ribbon'    ,'pm.show("ribbon"    ,"'+s+'")'],
           [ 0, ''          ,''                             ],
           [ 1, 'dots'      ,'pm.show("dots"      ,"'+s+'")'],
           [ 1, 'spheres'   ,'pm.show("spheres"   ,"'+s+'")'],
           [ 0, ''          ,''                             ],
           [ 1, 'mesh'      ,'pm.show("mesh"      ,"'+s+'")'],
           [ 1, 'surface'   ,'pm.show("surface"   ,"'+s+'")']]

def mol_hide(s):
   return [[ 2, 'Hide:'     ,''                             ],
           [ 1, 'lines'     ,'pm.hide("lines"     ,"'+s+'")'],
           [ 1, 'sticks'    ,'pm.hide("sticks"    ,"'+s+'")'],
           [ 1, 'ribbon'    ,'pm.hide("ribbon"    ,"'+s+'")'],
           [ 0, ''          ,''                             ],
           [ 1, 'dots'      ,'pm.hide("dots"      ,"'+s+'")'],
           [ 1, 'spheres'   ,'pm.hide("spheres"   ,"'+s+'")'],
           [ 0, ''          ,''                             ],
           [ 1, 'mesh'      ,'pm.hide("mesh"      ,"'+s+'")'],
           [ 1, 'surface'   ,'pm.hide("surface"   ,"'+s+'")'],
           [ 0, ''          ,''                             ],
           [ 1, 'everything','pm.hide("('+s+')")'           ]]

def dist_show(s):
   return [[ 2, 'Show:'     ,''                             ],
           [ 1, 'dashes'    ,'pm.show("dashes"    ,"'+s+'")'],
           [ 1, 'labels'    ,'pm.show("labels"    ,"'+s+'")']]

def dist_hide(s):
   return [[ 2, 'Hide:'     ,''                             ],
           [ 1, 'dashes'    ,'pm.hide("dashes"    ,"'+s+'")'],
           [ 1, 'labels'    ,'pm.hide("labels"    ,"'+s+'")']]

