7#A* -------------------------------------------------------------------
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
           [ 1, 'hydrogens' ,'pm.hide("('+s+' and hydro)")' ],
           [ 1, 'everything','pm.hide("('+s+')")'           ],
           [ 1, 'unselected','pm.hide("(not '+s+')")'      ],
           ]

def dist_show(s):
   return [[ 2, 'Show:'     ,''                             ],
           [ 1, 'dashes'    ,'pm.show("dashes"    ,"'+s+'")'],
           [ 1, 'labels'    ,'pm.show("labels"    ,"'+s+'")']]

def dist_hide(s):
   return [[ 2, 'Hide:'     ,''                             ],
           [ 1, 'dashes'    ,'pm.hide("dashes"    ,"'+s+'")'],
           [ 1, 'labels'    ,'pm.hide("labels"    ,"'+s+'")']]

def mol_color(s):
   return [[ 2, 'Color:'     ,''                        ],
           [ 1, '`292C`777H`229N`922O`950S`905*'     ,'pmu.cbag("'+s+'")'],
           [ 1, '`099C`777H`229N`922O`950S`905*'     ,'pmu.cbac("'+s+'")'],
           [ 1, '`990C`777H`229N`922O`950S`905*'     ,'pmu.cbay("'+s+'")'],
           [ 1, '`955C`777H`229N`922O`950S`905*'     ,'pmu.cbas("'+s+'")'],
           [ 1, '`777C`777H`229N`922O`950S`905*'     ,'pmu.cbaw("'+s+'")'],
           [ 0, ''          ,''                         ],
           [ 1, '`900red'         ,'pm.color("red","'+s+'")'  ],
           [ 1, '`090green'       ,'pm.color("green","'+s+'")'  ],
           [ 1, '`009blue'        ,'pm.color("blue","'+s+'")'  ],
           [ 1, '`990yellow'      ,'pm.color("yellow","'+s+'")'  ],
           [ 1, '`909violet'  ,'pm.color("violet","'+s+'")'  ],
           [ 1, '`099cyan'  ,'pm.color("cyan","'+s+'")'  ],           
           [ 1, '`955salmon'      ,'pm.color("salmon","'+s+'")'  ],
           [ 1, '`595lime' ,'pm.color("lime","'+s+'")'  ],
           [ 1, '`559slate'  ,'pm.color("slate","'+s+'")'  ],
           [ 1, '`905magenta' ,'pm.color("magenta","'+s+'")'  ],
           [ 1, '`950orange'      ,'pm.color("orange","'+s+'")'  ],
           [ 1, '`059marine'      ,'pm.color("marine","'+s+'")'  ],
           [ 1, '`551olive'   ,'pm.color("olive","'+s+'")'  ],
           [ 1, '`515purple'  ,'pm.color("purple","'+s+'")'  ],
           [ 1, '`155teal'  ,'pm.color("teal","'+s+'")'  ],
           [ 1, '`151forest'  ,'pm.color("forest","'+s+'")'  ],
           [ 1, '`999white'       ,'pm.color("white","'+s+'")'  ],
           [ 1, '`555grey'    ,'pm.color("grey","'+s+'")'  ]
           ]

def dist_color(s):
   return [[ 2, 'Color:'      ,''                            ],      
           [ 1, 'yellow'      ,'pm.color("yellow","'+s+'")'  ],
           [ 1, 'green'       ,'pm.color("green","'+s+'")'   ],
           [ 1, 'orange'      ,'pm.color("orange","'+s+'")'  ],
           [ 1, 'magenta'     ,'pm.color("magenta","'+s+'")' ]
           ]

def sele_action(s):
   return [[ 2, 'Actions:'     ,''                      ],     
           [ 1, 'Origin'   ,'pm.origin("'+s+'")'    ],
           [ 1, 'Zoom'         ,'pm.zoom("'+s+'")'      ],
           [ 1, 'Orient'       ,'pm.orient("'+s+'")'      ],
           [ 0, ''          ,''                         ],
           [ 1, 'Delete'       ,'pm.delete("'+s+'")'    ],
           [ 0, ''          ,''                         ],
           [ 1, 'Duplicate'    ,'pm.select("('+s+')")'  ],
           [ 0, ''          ,''                         ],
           [ 1, 'Full Residues' ,'pm.select("'+s+'","(byres '+s+')")' ],
           [ 1, 'Full Residues + 3 A' ,'pm.select("'+s+'","(byres ('+s+' expand 3))")' ],
           [ 1, 'Full Residues + 7 A' ,'pm.select("'+s+'","(byres ('+s+' expand 7))")' ],
           [ 1, 'Full Residues + 10 A' ,'pm.select("'+s+'","(byres ('+s+' expand 10))")' ],
           [ 0, ''          ,''                         ],
           [ 1, 'Expand By 3 A'  ,'pm.select("'+s+'","('+s+' expand 3)")' ],
           [ 1, 'Expand By 7 A'  ,'pm.select("'+s+'","('+s+' expand 7)")' ],
           [ 1, 'Expand By 10 A'  ,'pm.select("'+s+'","('+s+' expand 10)")' ],
           [ 0, ''          ,''                         ],
           [ 1, 'Invert'  ,'pm.select("'+s+'","(not '+s+')")' ],
           ]

def mol_action(s):
   return [[ 2, 'Actions:'     ,''                      ],     
           [ 1, 'Origin'   ,'pm.origin("'+s+'")'    ],
           [ 1, 'Zoom'         ,'pm.zoom("'+s+'")'      ],
           [ 1, 'Orient'       ,'pm.orient("'+s+'")'      ],
           [ 0, ''          ,''                         ],
           [ 1, 'Delete'       ,'pm.delete("'+s+'")'    ],
           ]

def all_action(s):
   return [[ 2, 'Actions:'     ,''                      ],     
           [ 1, 'Set Origin'   ,'pm.origin("'+s+'")'    ],
           [ 1, 'Zoom'         ,'pm.zoom("'+s+'")'      ],
           [ 0, ''          ,''                         ],
           [ 1, 'Delete'       ,'pm.delete("all")'    ]
           ]

