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

import cmd

def mol_show(s):
   return [[ 2, 'Show:'      , ''                               ],
           [ 1, 'lines'      , 'cmd.show("lines"     ,"'+s+'")' ],
           [ 1, 'nonbonded'  , 'cmd.show("nonbonded" ,"'+s+'")' ],
           [ 1, 'sticks'     , 'cmd.show("sticks"    ,"'+s+'")' ],
           [ 1, 'ribbon'     , 'cmd.show("ribbon"    ,"'+s+'")' ],
           [ 1, 'cartoon'    , 'cmd.show("cartoon"   ,"'+s+'")' ],
           [ 0, ''           , ''                               ],
           [ 1, 'labels'     , 'cmd.show("labels"    ,"'+s+'")' ],
           [ 1, 'cell'       , 'cmd.show("cell"      ,"'+s+'")' ],
           [ 0, ''           , ''                               ],
           [ 1, 'dots'       , 'cmd.show("dots"      ,"'+s+'")' ],
           [ 1, 'spheres'    , 'cmd.show("spheres"   ,"'+s+'")' ],
           [ 1, 'nb_spheres' , 'cmd.show("nb_spheres","'+s+'")' ],
           [ 0, ''           , ''                               ],
           [ 1, 'mesh'       , 'cmd.show("mesh"      ,"'+s+'")' ],
           [ 1, 'surface'    , 'cmd.show("surface"   ,"'+s+'")' ],
           [ 0, ''           , ''                               ],           
           [ 1, 'main chain' , 'cmd.show("lines","((byres ('+s+'))&n;ca,c,n,o,h)")' ],
           [ 1, 'side chain' , 'cmd.show("lines","((byres ('+s+'))&(!n;c,n,o,h))")' ],
           ]

def mol_hide(s):
   return [[ 2, 'Hide:'     , ''                                ],
           [ 1, 'everything', 'cmd.hide("everything","'+s+'")'  ],
           [ 0, ''          , ''                                ],
           [ 1, 'lines'     , 'cmd.hide("lines"     ,"'+s+'")'  ],
           [ 1, 'nonbonded' , 'cmd.hide("nonbonded" ,"'+s+'")'  ],           
           [ 1, 'sticks'    , 'cmd.hide("sticks"    ,"'+s+'")'  ],
           [ 1, 'ribbon'    , 'cmd.hide("ribbon"    ,"'+s+'")'  ],
           [ 1, 'cartoon'   , 'cmd.hide("cartoon"   ,"'+s+'")' ],           
           [ 0, ''          , ''                                ],
           [ 1, 'labels'    , 'cmd.hide("labels"    ,"'+s+'")'  ],
           [ 1, 'cell'      , 'cmd.hide("cell"      ,"'+s+'")' ],
           [ 0, ''          , ''                                ],
           [ 1, 'dots'      , 'cmd.hide("dots"      ,"'+s+'")'  ],
           [ 1, 'spheres'   , 'cmd.hide("spheres"   ,"'+s+'")'  ],
           [ 1, 'nb_spheres', 'cmd.hide("nb_spheres","'+s+'")'  ],           
           [ 0, ''          , ''                                ],
           [ 1, 'mesh'      , 'cmd.hide("mesh"      ,"'+s+'")'  ],
           [ 1, 'surface'   , 'cmd.hide("surface"   ,"'+s+'")'  ],
           [ 0, ''          , ''                                ],
           [ 1, 'main chain', 'cmd.hide("((byres ('+s+'))&n;c,n,o,h)")' ],
           [ 1, 'side chain', 'cmd.hide("((byres ('+s+'))&!n;ca,c,n,o,h)")' ],
           [ 1, 'waters'    , 'cmd.hide("(resn HOH+WAT and ('+s+'))")'     ],                      
           [ 0, ''          , ''                                ],
           [ 1, 'hydrogens' , 'cmd.hide("('+s+' and hydro)")'   ],
           [ 0, ''          , ''                                ],           
           [ 1, 'unselected', 'cmd.hide("(not '+s+')")'         ],
           ]

def dist_show(s):
   return [[ 2, 'Show:'     , ''                               ],
           [ 1, 'dashes'    , 'cmd.show("dashes"    ,"'+s+'")' ],
           [ 1, 'labels'    , 'cmd.show("labels"    ,"'+s+'")' ]
          ]   

def dist_hide(s):
   return [[ 2, 'Hide:'     , ''                                ],
           [ 1, 'dashes'    , 'cmd.hide("dashes"    ,"'+s+'")'  ],
           [ 1, 'labels'    , 'cmd.hide("labels"    ,"'+s+'")'  ]
          ]

def cgo_show(s):
   return [[ 2, 'Show:'     , ''                               ],
           [ 1, 'cgo'    , 'cmd.show("cgo"    ,"'+s+'")' ],
          ]   

def cgo_hide(s):
   return [[ 2, 'Hide:'     , ''                                ],
           [ 1, 'cgo'    , 'cmd.hide("cgo"    ,"'+s+'")'  ],
          ]

def simple_show(s):
   return [[ 2, 'Show:'       , ''                             ],
           [ 1, 'everything'  , 'cmd.show("everything","'+s+'")'          ]]

def simple_hide(s):
   return [[ 2, 'Hide:'     ,''                                ],
           [ 1, 'everything'    ,'cmd.hide("everything","'+s+'")'        ]]

def mesh_show(s):
   return [[ 2, 'Show:'       , ''                             ],
           [ 1, 'mesh'        , 'cmd.show("mesh","'+s+'")'     ],           
           [ 1, 'cell'        , 'cmd.show("cell","'+s+'")'     ],
           [ 1, 'everything'  , 'cmd.enable("'+s+'")'          ]]

def mesh_hide(s):
   return [[ 2, 'Hide:'       , ''                             ],
           [ 1, 'mesh'        , 'cmd.hide("mesh","'+s+'")'     ],                      
           [ 1, 'cell'        , 'cmd.hide("cell","'+s+'")'      ],           
           [ 1, 'everything'  , 'cmd.disable("'+s+'")'          ]]

def by_elem(s):
   return [
      [ 2, 'By Element:'     ,''                               ],
      [ 1, '\\292C\\777H\\229N\\922O\\950S\\905*'  , 'util.cbag("'+s+'")'],
   [ 1, '\\099C\\777H\\229N\\922O\\950S\\905*'  , 'util.cbac("'+s+'")'],
   [ 1, '\\909C\\777H\\229N\\922O\\950S\\905*'  , 'util.cbam("'+s+'")'],           
   [ 1, '\\990C\\777H\\229N\\922O\\950S\\905*'  , 'util.cbay("'+s+'")'],
   [ 1, '\\955C\\777H\\229N\\922O\\950S\\905*'  , 'util.cbas("'+s+'")'],
   [ 1, '\\777C\\777H\\229N\\922O\\950S\\905*'  , 'util.cbaw("'+s+'")'],
   [ 1, '\\559C\\777H\\229N\\922O\\950S\\905*'  , 'util.cbab("'+s+'")'],
   [ 1, '\\972C\\777H\\229N\\922O\\950S\\905*'  , 'util.cbao("'+s+'")'],
   [ 1, ' \\777H\\229N\\922O\\950S\\905*'  , 'util.cnc("'+s+'")']]                     

def by_ss(s):
   return [
            [ 2, 'By Secondary Structure:'     ,''                               ],
   [ 1, '\\900Helix \\990Sheet \\090Loop'  , 'util.cbss("'+s+'","red","yellow","green")'],
   [ 1, '\\099Helix \\909Sheet \\955Loop'  , 'util.cbss("'+s+'","cyan","magenta","salmon")'],
           ]

def spectrum(s):
   return [
            [ 2, 'Spectrum:'     ,''                               ],
           [ 1, '\\900r\\950a\\990i\\090n\\099b\\059o\\009w\\888(e. c)',
             'cmd.spectrum("count",selection="('+s+')&e. c")'],           
           [ 1, '\\900r\\950a\\990i\\090n\\099b\\059o\\009w\\888(*/ca)',
             'cmd.spectrum("count",selection="('+s+')&*/ca")'],
           [ 1, '\\900r\\950a\\990i\\090n\\099b\\059o\\009w',
             'cmd.spectrum("count",selection="'+s+'",byres=1)'],
           [ 0, ''                                , ''                 ],
           [ 1, 'b-factors'   , 'cmd.spectrum("b",selection=("'+s+'"))'         ],           
           ]

def by_chain(s):
   return [
      [ 2, 'By Chain:'     ,''                               ],
           [ 1, '\\900b\\950y \\090c\\099h\\059a\\009i\\705n\\888(e. c)',
             'util.color_chains("('+s+' and elem c)")'],
           [ 1, '\\900b\\950y \\090c\\099h\\059a\\009i\\705n\\888(*/ca)',
             'util.color_chains("('+s+' and name ca)")'],
           [ 1, '\\900b\\950y \\090c\\099h\\059a\\009i\\705n',
             'util.color_chains("('+s+')")'],
                 [ 0, ''                                , ''                 ],
           [ 1, '\\900c\\950h\\990a\\090i\\099n\\059b\\009o\\705w\\888s',
             'util.chainbow("('+s+')")'],                                 
      ]
def mol_color(s):
   return [[ 2, 'Color:'     ,''                               ],
           [ 1, 'by chain' , by_chain(s) ],
           [ 1, 'by element'  , by_elem(s) ],
           [ 1, 'by ss  '  , by_ss(s) ],
           [ 1, '\\900s\\950p\\990e\\090c\\099t\\059r\\009u\\555m', spectrum(s) ],
           [ 0, ''                                , ''                 ],           
           [ 1, '\\900red'         , 'cmd.color("red","'+s+'")'     ],
           [ 1, '\\090green'       , 'cmd.color("green","'+s+'")'   ],
           [ 1, '\\009blue'        , 'cmd.color("blue","'+s+'")'    ],
           [ 1, '\\990yellow'      , 'cmd.color("yellow","'+s+'")'  ],
           [ 1, '\\909magenta'     , 'cmd.color("magenta","'+s+'")' ],
           [ 1, '\\099cyan'        , 'cmd.color("cyan","'+s+'")'    ],           
           [ 1, '\\955salmon'      , 'cmd.color("salmon","'+s+'")'  ],
           [ 1, '\\595lime'        , 'cmd.color("lime","'+s+'")'    ],
           [ 1, '\\967pink'  ,'cmd.color("pink","'+s+'")'  ],
           [ 1, '\\559slate'       , 'cmd.color("slate","'+s+'")'   ],
           [ 1, '\\949violet'      , 'cmd.color("violet","'+s+'")'  ],
           [ 1, '\\950orange'      , 'cmd.color("orange","'+s+'")'  ],
           [ 1, '\\059marine'      , 'cmd.color("marine","'+s+'")'  ],
           [ 1, '\\905hotpink' ,'cmd.color("hotpink","'+s+'")'  ],
#           [ 1, '\\551olive'       , 'cmd.color("olive","'+s+'")'   ],
#           [ 1, '\\626purple'      , 'cmd.color("purple","'+s+'")'  ],
#           [ 1, '\\266teal'        , 'cmd.color("teal","'+s+'")'    ],
#           [ 1, '\\151forest'      , 'cmd.color("forest","'+s+'")'  ],
#           [ 1, '\\611firebrick'   , 'cmd.color("firebrick","'+s+'")'  ],
#           [ 1, '\\631chocolate'   , 'cmd.color("chocolate","'+s+'")'  ],           
           [ 1, '\\999white'       , 'cmd.color("white","'+s+'")'   ],
           [ 1, '\\987wheat'       , 'cmd.color("wheat","'+s+'")'   ],
           [ 1, '\\555grey'        , 'cmd.color("grey","'+s+'")'    ],
           [ 1, '\\222black'    ,'cmd.color("grey","'+s+'")'  ]
           ]

def general_color(s):
   return [[ 2, 'Color:'     ,''                        ],
           [ 1, '\\900red'         ,'cmd.color("red","'+s+'")'  ],
           [ 1, '\\090green'       ,'cmd.color("green","'+s+'")'  ],
           [ 1, '\\009blue'        ,'cmd.color("blue","'+s+'")'  ],
           [ 1, '\\990yellow'      ,'cmd.color("yellow","'+s+'")'  ],
           [ 1, '\\909magenta' ,'cmd.color("magenta","'+s+'")'  ],
           [ 1, '\\099cyan'  ,'cmd.color("cyan","'+s+'")'  ],           
           [ 1, '\\955salmon'      ,'cmd.color("salmon","'+s+'")'  ],
           [ 1, '\\595lime' ,'cmd.color("lime","'+s+'")'  ],
           [ 1, '\\967pink'  ,'cmd.color("pink","'+s+'")'  ],
           [ 1, '\\559slate'  ,'cmd.color("slate","'+s+'")'  ],
           [ 1, '\\949violet'  ,'cmd.color("violet","'+s+'")'  ],
           [ 1, '\\950orange'      ,'cmd.color("orange","'+s+'")'  ],
           [ 1, '\\059marine'      ,'cmd.color("marine","'+s+'")'  ],
           [ 1, '\\905hotpink' ,'cmd.color("hotpink","'+s+'")'  ],
#           [ 1, '\\551olive'   ,'cmd.color("olive","'+s+'")'  ],
#           [ 1, '\\626purple'  ,'cmd.color("purple","'+s+'")'  ],
#           [ 1, '\\266teal'  ,'cmd.color("teal","'+s+'")'  ],
#           [ 1, '\\151forest'  ,'cmd.color("forest","'+s+'")'  ],
#           [ 1, '\\611firebrick'   , 'cmd.color("firebrick","'+s+'")'  ],
#           [ 1, '\\631chocolate'   , 'cmd.color("chocolate","'+s+'")'  ],
           [ 1, '\\116density'     ,'cmd.color("density","'+s+'")'  ],
           [ 1, '\\999white'       ,'cmd.color("white","'+s+'")'  ],
           [ 1, '\\987wheat'       , 'cmd.color("wheat","'+s+'")'   ],
           [ 1, '\\555grey'    ,'cmd.color("grey","'+s+'")'  ],
           [ 1, '\\222black'    ,'cmd.color("grey","'+s+'")'  ]
           ]

def preset_ligand_sites(s):
   return [[ 2, 'Ligand Sites:', ''],
           [ 1, 'solid'   , 'preset.ligand_sites("'+s+'")'          ],
           [ 1, 'solid (better)'   , 'preset.ligand_sites_hq("'+s+'")'          ],
           [ 0, '', ''],
           [ 1, 'transparent'   , 'preset.ligand_sites_trans("'+s+'")'          ],
           [ 1, 'transparent (better)'   , 'preset.ligand_sites_trans_hq("'+s+'")'          ],
           [ 0, '', ''],
           [ 1, 'dot surface'   , 'preset.ligand_sites_dots("'+s+'")'          ],
           [ 0, '', ''],
           [ 1, 'mesh surface'   , 'preset.ligand_sites_mesh("'+s+'")'          ]]

def presets(s):
   return [[ 2, 'Preset:'       ,''                        ],     
           [ 1, 'simple'   ,'preset.simple("'+s+'")'          ],
           [ 1, 'simple (no solvent)'   ,'preset.simple_no_solv("'+s+'")'          ],           
           [ 1, 'technical'   , 'preset.technical("'+s+'")'          ],
           [ 1, 'ligands'   , 'preset.ligands("'+s+'")'          ],
           [ 1, 'ligand_sites'   , preset_ligand_sites(s)         ],
           [ 1, 'pretty ', 'preset.pretty("'+s+'")'          ],
           [ 1, 'pretty (with solvent)'     , 'preset.pretty_solv("'+s+'")'          ],
           [ 1, 'publication '   , 'preset.publication("'+s+'")'          ],
           [ 1, 'publication (with solvent)'   , 'preset.pub_solv("'+s+'")'          ],
           [ 0, ''               ,''                             ],                      
           [ 1, 'default'   ,'preset.default("'+s+'")'          ],           
           ]

#def hydrogens(s):
#   return [[ 2, 'Hydrogens:'       ,''                        ],     
#           [ 1, 'add'   ,'cmd.h_add("'+s+'")'          ],           
#           [ 1, 'remove'   ,'cmd.remove("('+s+') and hydro")'          ],
#           ]

def state(s):
   return [[ 2, 'State:'       ,''                        ],
           [ 1, 'freeze'  ,'cmd.set("state",cmd.get_state(),"'+s+'")'        ],
           [ 1, 'thaw'  ,'cmd.set("state",0,"'+s+'")'        ],           
           ]
   
def movement(s):
   return [[ 2, 'Movement:'       ,''                        ],     
           [ 1, 'protect'   ,'cmd.protect("'+s+'")'          ],
           [ 1, 'deprotect'   ,'cmd.deprotect("'+s+'")'          ],           
           ]
   
def selection(s):
   return [[ 2, 'Selection:'       ,''                        ],     
           [ 1, 'mask'   ,'cmd.mask("'+s+'")'          ],
           [ 1, 'unmask'   ,'cmd.unmask("'+s+'")'          ],           
           ]

def compute(s):
   return [[ 2, 'Compute:'       ,''                        ],     
           [ 1, 'atom count'   ,'cmd.count_atoms("'+s+'",quiet=0)'          ],
           [ 0, ''               ,''                             ],           
           [ 1, 'formal charge sum'   ,'util.sum_formal_charges("'+s+'",quiet=0)'          ],
           [ 1, 'partial charges sum'   ,'util.sum_partial_charges("'+s+'",quiet=0)'          ],                      
           ]

def vacuum(s):
   return [[ 2, 'Vacuum Electrostatics:'       ,''                        ],
           [ 2, '\\955WARNING:\\595 unvalidated and experimental!', '' ],
           [ 1, 'protein surface potential (absolute)', 'util.protein_vacuum_esp("'+s+'",quiet=0)'          ],
           [ 1, 'protein surface potential (relative)', 'util.protein_vacuum_esp("'+s+'",absolute=0,quiet=0)'          ],
           [ 2, '\\955NOTE:\\559 vacuum electrostatic potentials for', '' ],
           [ 2, '\\559macromolecules are not really appropriate', '' ],
           [ 2, '\\559or meaningful, so view with extreme caution!', '' ],
           ]
   
def mol_assign(s):
   return [[ 2, 'Assign:'       ,''                        ],     
           [ 1, 'Amber 99 atomic properties',  'util.assign_amber99("'+s+'")' ],
           ]
   
def mol_generate(s):
   return [[ 2, 'Generate:'       ,''                        ],
           [ 1, 'vacuum electrostatics', vacuum(s) ],           
#           [ 1, 'assign', mol_assign(s) ],
           ]
   
def invert(s):
   return [[ 2, 'Invert:'       ,''                        ],     
           [ 1, 'within object(s)'     ,'cmd.select("'+s+'","((byobj '+s+') and not '+s+')",show=1)'    ],
           [ 1, 'within segments(s)'     ,'cmd.select("'+s+'","((byseg '+s+') and not '+s+')",show=1)'    ],           
           [ 1, 'within chains(s)'     ,'cmd.select("'+s+'","((bychain '+s+') and not '+s+')",show=1)'    ],
           [ 1, 'within residues(s)'   ,'cmd.select("'+s+'","((byres '+s+') and not '+s+')",show=1)'    ],                      
           [ 0, ''               ,''                             ],
           [ 1, 'within molecule(s)'     ,'cmd.select("'+s+'","((bymol '+s+') and not '+s+')",show=1)'    ],
           [ 0, ''               ,''                             ],
           [ 1, 'within any'     ,'cmd.select("'+s+'","(not '+s+')",show=1)'    ],
           ]

def complete(s):
   return [[ 2, 'Complete:'       ,''                        ],     

           [ 1, 'resides'  ,'cmd.select("'+s+'","(byres '+s+')",show=1)'      ],
           [ 1, 'chains'  ,'cmd.select("'+s+'","(bychain '+s+')",show=1)'      ],
           [ 1, 'segments'  ,'cmd.select("'+s+'","(byseg '+s+')",show=1)'      ],
           [ 1, 'objects'  ,'cmd.select("'+s+'","(byobj '+s+')",show=1)'      ],
           [ 0, ''               ,''                             ],           
           [ 1, 'molecules'  ,'cmd.select("'+s+'","(bymol '+s+')",show=1)'      ],
           [ 0, ''               ,''                             ],
           [ 1, 'C-alphas'  ,'cmd.select("'+s+'","(bycalpha '+s+')",show=1)'      ],           
           ]

def expand(s):
   return [[ 2, 'Expand:'       ,''                        ],     
           [ 1, 'by 4 A'  ,'cmd.select("'+s+'","('+s+' expand 4)",show=1)' ],
           [ 1, 'by 6 A'  ,'cmd.select("'+s+'","('+s+' expand 6)",show=1)' ],           
           [ 1, 'by 8 A'  ,'cmd.select("'+s+'","('+s+' expand 8)",show=1)' ],
           [ 1, 'by 12 A'  ,'cmd.select("'+s+'","('+s+' expand 12)",show=1)' ],
           [ 1, 'by 16 A'  ,'cmd.select("'+s+'","('+s+' expand 16)",show=1)' ],
           [ 1, 'by 20 A'  ,'cmd.select("'+s+'","('+s+' expand 20)",show=1)' ],           
           [ 0, ''               ,''                             ],
           [ 1, 'by 4 A, residues'  ,'cmd.select("'+s+'","(byres ('+s+' expand 4))",show=1)' ],
           [ 1, 'by 6 A, residues'  ,'cmd.select("'+s+'","(byres ('+s+' expand 6))",show=1)' ],
           [ 1, 'by 8 A, residues'   ,'cmd.select("'+s+'","(byres ('+s+' expand 8))",show=1)' ],
           [ 1, 'by 12 A, residues'  ,'cmd.select("'+s+'","(byres ('+s+' expand 12))",show=1)' ],
           [ 1, 'by 16 A, residues'  ,'cmd.select("'+s+'","(byres ('+s+' expand 16))",show=1)' ],
           [ 1, 'by 20 A, residues'  ,'cmd.select("'+s+'","(byres ('+s+' expand 20))",show=1)' ],
           ]

def around(s):
   return [[ 2, 'Around:'       ,''                        ],     
           [ 1, 'atoms within 4 A'  ,'cmd.select("'+s+'","('+s+' around 4)",show=1)' ],
           [ 1, 'atoms within 6 A'  ,'cmd.select("'+s+'","('+s+' around 6)",show=1)' ],           
           [ 1, 'atoms within 8 A'  ,'cmd.select("'+s+'","('+s+' around 8)",show=1)' ],
           [ 1, 'atoms within 12 A'  ,'cmd.select("'+s+'","('+s+' around 12)",show=1)' ],
           [ 1, 'atoms within 16 A'  ,'cmd.select("'+s+'","('+s+' around 16)",show=1)' ],
           [ 1, 'atoms within 20 A'  ,'cmd.select("'+s+'","('+s+' around 20)",show=1)' ],                      
           [ 0, ''               ,''                             ],           
           [ 1, 'residues within 4 A'  ,'cmd.select("'+s+'","(byres ('+s+' around 4))",show=1)' ],
           [ 1, 'residues within 6 A'  ,'cmd.select("'+s+'","(byres ('+s+' around 6))",show=1)' ],
           [ 1, 'residues within 8 A'  ,'cmd.select("'+s+'","(byres ('+s+' around 8))",show=1)' ],
           [ 1, 'residues within 16 A'  ,'cmd.select("'+s+'","(byres ('+s+' around 12))",show=1)' ],
           [ 1, 'residues within 20 A'  ,'cmd.select("'+s+'","(byres ('+s+' around 20))",show=1)' ],                                 
           ]
   
def extend(s):
   return [[ 2, 'Extend:'       ,''                        ],     
           [ 1, 'by 1 bond'  ,'cmd.select("'+s+'","('+s+' extend 1)",show=1)' ],
           [ 1, 'by 2 bonds'  ,'cmd.select("'+s+'","('+s+' extend 2)",show=1)' ],           
           [ 1, 'by 3 bonds'  ,'cmd.select("'+s+'","('+s+' extend 3)",show=1)' ],
           [ 1, 'by 4 bonds'  ,'cmd.select("'+s+'","('+s+' extend 4)",show=1)' ],
           [ 1, 'by 5 bonds'  ,'cmd.select("'+s+'","('+s+' extend 5)",show=1)' ],
           [ 1, 'by 6 bonds'  ,'cmd.select("'+s+'","('+s+' extend 6)",show=1)' ],           
           [ 0, ''               ,''                             ],
           [ 1, 'by 1 bond, residues'  ,'cmd.select("'+s+'","(byres ('+s+' extend 1))",show=1)' ],
           [ 1, 'by 2 bonds, residues'  ,'cmd.select("'+s+'","(byres ('+s+' extend 2))",show=1)' ],
           [ 1, 'by 3 bonds, residues'   ,'cmd.select("'+s+'","(byres ('+s+' extend 3))",show=1)' ],
           [ 1, 'by 4 bonds, residues'  ,'cmd.select("'+s+'","(byres ('+s+' extend 4))",show=1)' ],
           [ 1, 'by 5 bonds, residues'  ,'cmd.select("'+s+'","(byres ('+s+' extend 5))",show=1)' ],
           [ 1, 'by 6 bonds, residues'  ,'cmd.select("'+s+'","(byres ('+s+' extend 6))",show=1)' ],
           ]

def sele_action(s):
   return [[ 2, 'Actions:'       ,''                        ],     
           [ 1, 'delete selection', 'cmd.delete("'+s+'")'          ],
           [ 1, 'rename selection', 'cmd.wizard("renaming","'+s+'")'          ],
           [ 0, ''               ,''                             ],
           [ 1, 'zoom'           ,'cmd.zoom("'+s+'")'            ],
           [ 1, 'center'         ,'cmd.center("'+s+'")'            ],           
           [ 1, 'origin'         ,'cmd.origin("'+s+'")'          ],
           [ 1, 'orient'         ,'cmd.orient("'+s+'")'          ],
           [ 0, ''               ,''                             ],
           [ 1, 'preset'         ,presets(s)         ],
           [ 0, ''               ,''                             ],
           [ 1, 'remove atoms'   ,'cmd.remove("'+s+'");cmd.delete("'+s+'")'          ],
           [ 0, ''               ,''                             ],
           [ 1, 'around'         , around(s)         ],           
           [ 1, 'expand'         , expand(s)         ],
           [ 1, 'extend'         , extend(s)         ],
           [ 1, 'invert'         , invert(s)         ],
           [ 1, 'complete'       , complete(s)         ],
           [ 0, ''          ,''                                              ],
           [ 1, 'duplicate selection'      ,'cmd.select("'+s+'")'          ],
           [ 1, 'create object'  ,'cmd.create(None,"'+s+'")'     ],           
           [ 1, 'find polar contacts'  ,
             'cmd.dist("'+s+'_polar_conts","'+s+'","'+s+'",quiet=1,mode=2,labels=0)'
             ],                      
           [ 0, ''          ,''                                  ],
           [ 1, 'selection'      , selection(s)         ],           
           [ 1, 'movement'       , movement(s)         ],
           [ 1, 'compute'        , compute(s)         ],           
           ]


def mol_action(s):
   return [[ 2, 'Actions:'     , ''                       ],     
           [ 1, 'zoom'         , 'cmd.zoom("'+s+'")'      ],
           [ 1, 'center'         ,'cmd.center("'+s+'")'            ],
           [ 1, 'origin'       , 'cmd.origin("'+s+'")'    ],
           [ 1, 'orient'       , 'cmd.orient("'+s+'")'    ],
           [ 0, ''          ,''                                              ],
           [ 1, 'preset'  , presets(s)       ],
           [ 0, ''          ,''                                              ],           
           [ 1, 'generate'  , mol_generate(s)       ],           
           [ 0, ''               ,''                             ],
           [ 1, 'assign sec. struc.'  ,'cmd.dss("'+s+'")'        ],
           [ 1, 'find polar contacts'  ,
             'cmd.dist("'+s+'_polar_conts","'+s+'","'+s+'",quiet=1,mode=2,labels=0)'
             ],                      
           [ 0, ''             , ''                       ],
           [ 1, 'rename object', 'cmd.wizard("renaming","'+s+'")'          ],
           [ 1, 'duplicate object'    ,'cmd.create(None,"'+s+'")'     ],           
           [ 1, 'delete object'       , 'cmd.delete("'+s+'")'    ],
           [ 0, ''          ,''                                              ],
           [ 1, 'add hydrogens'  ,'cmd.h_add("'+s+'")'     ],           
           [ 1, 'remove hydrogens'  ,'cmd.remove("(hydro and ('+s+'))")'     ],           
           [ 1, 'remove waters'  ,'cmd.remove("(resn HOH+WAT and ('+s+'))")'     ],
           [ 0, ''          ,''                                              ],
           [ 1, 'state'      , state(s)         ],                      
           [ 1, 'selection'      , selection(s)         ],           
           [ 1, 'movement'       , movement(s)         ],           
           [ 1, 'compute'        , compute(s)         ],
           ]

def simple_action(s):
   return [[ 2, 'Actions:'     , ''                       ],
           [ 1, 'zoom'         , 'cmd.zoom("'+s+'")'      ],
           [ 1, 'center'       , 'cmd.center("'+s+'")'    ],           
           [ 1, 'origin'       , 'cmd.origin("'+s+'")'    ],
           [ 0, ''             , ''                       ],
           [ 1, 'delete'       , 'cmd.delete("'+s+'")'    ],
           ]

def test1(s):
      return [[ 2, 'Test1:'     , ''                      ],     
           [ 1, 'zoom'         , 'cmd.zoom("all")'     ],
           [ 1, 'center'   , 'cmd.center("all")'   ],           
           [ 1, 'origin'   , 'cmd.origin("all")'   ],
           ]

def test2(s):
      return [[ 2, 'Test2:'     , ''                      ],     
           [ 1, 'zoom'         , 'cmd.zoom("all")'     ],
           [ 1, 'center'   , 'cmd.center("all")'   ],           
           [ 1, 'origin'   , 'cmd.origin("all")'   ],
           ]

def all_action(s):
   return [[ 2, 'Actions:'     , ''                      ],     
           [ 1, 'zoom'         , 'cmd.zoom("all")'     ],
           [ 1, 'center'   , 'cmd.center("all")'   ],           
           [ 1, 'origin'   , 'cmd.origin("all")'   ],
           [ 0, ''             , ''                      ],           
           [ 1, 'preset'  , presets(s)     ],           
           [ 0, ''             , ''                      ],           
           [ 1, 'delete everything'  , 'cmd.delete("all")'     ],           
           [ 0, ''          ,''                                              ],
           [ 1, 'add hydrogens' ,'cmd.h_add("'+s+'")'     ],           
           [ 1, 'remove hydrogens'  ,'cmd.remove("(elem h and ('+s+'))")'     ],
           [ 1, 'remove waters'  ,'cmd.remove("(resn HOH+WAT and ('+s+'))")'     ],                      
           [ 0, ''          ,''                                              ],
           [ 1, 'selection'      , selection(s)         ],                      
           [ 1, 'movement'       , movement(s)         ],
           [ 1, 'compute'        , compute(s)         ],                      
           ]


def mol_labels(s):
   return [[ 2, 'Labels:'        , ''                                  ],
           [ 1, 'clear'          , 'cmd.label("'+s+'","\'\'")'         ],
           [ 0, ''               , ''                                  ],
           [ 1, 'residues'       ,
  """cmd.label('''(name ca+C1*+C1' and (byres("""+s+""")))''','''"%s-%s"%(resn,resi)''')"""  ],
           [ 1, 'chains'       ,   'util.label_chains("'+s+'")'  ],
           [ 1, 'segments'       ,   'util.label_segments("'+s+'")'  ],           
           [ 0, ''               , ''                                  ],           
           [ 1, 'atom name'      , 'cmd.label("'+s+'","name")'         ],
           [ 1, 'element symbol' , 'cmd.label("'+s+'","elem")'         ],           
           [ 1, 'residue name'  , 'cmd.label("'+s+'","resn")'         ],
           [ 1, 'residue identifier'    , 'cmd.label("'+s+'","resi")'         ],
           [ 1, 'chain identifier' , 'cmd.label("'+s+'","chain")'         ],
           [ 1, 'segment identifier'       , 'cmd.label("'+s+'","segi")'         ],           
           [ 0, ''               , ''                                  ],
           [ 1, 'b-factor'       , 'cmd.label("'+s+'","\'%1.2f\'%b")'  ],
           [ 1, 'occupancy'       , 'cmd.label("'+s+'","\'%1.2f\'%q")'  ],
           [ 1, 'vdw radius'       , 'cmd.label("'+s+'","\'%1.2f\'%vdw")'  ],
           [ 0, ''               , ''                                  ],
           [ 1, 'partial charge(.2f)' ,            
  'cmd.label("'+s+'","\'%.2f\'%partial_charge")'                      ],
           [ 1, 'partial charge(.4f)' , 
  'cmd.label("'+s+'","\'%.4f\'%partial_charge")'                      ],
           [ 1, 'formal charge' , 
  'cmd.label("'+s+'","\'%d\'%formal_charge")'                      ],
           [ 1, 'bohr radius'       , 'cmd.label("'+s+'","\'%1.2f\'%bohr")'  ],                                 
           [ 0, ''               , ''                                  ],
           [ 1, 'text type'      , 'cmd.label("'+s+'","text_type")'    ],
           [ 1, 'numeric type'   , 'cmd.label("'+s+'","numeric_type")' ],
           [ 0, ''               , ''                                  ],      
           [ 1, 'atom i.d.'             , 'cmd.label("'+s+'","id")' ],
           [ 1, 'atom i.d.(+1)'           , 'cmd.label("'+s+'","id+1")' ],
           ]


def mol_view(s):
   return [
      [ 1, 'zoom'           ,'cmd.zoom("'+s+'")'            ],
      [ 1, 'center'           ,'cmd.center("'+s+'")'            ],
      [ 1, 'origin'           ,'cmd.origin("'+s+'")'            ],
      [ 1, 'orient'           ,'cmd.orient("'+s+'")'            ],
      ]

def all_option(s):
   return [
      [ 2, '(all)'      , '' ],
      [ 1, 'show'      , mol_show(s) ],
      [ 1, 'hide'      , mol_hide(s) ],
      [ 1, 'color'      , mol_color(s) ],
      [ 1, 'view'      , mol_view(s) ],
      [ 1, 'preset'      , presets(s) ],
      [ 0, ''             , ''                      ],
      [ 1, 'zoom'           ,'cmd.zoom("'+s+'")'            ],
      [ 1, 'center'           ,'cmd.center("'+s+'")'            ],
      [ 1, 'origin'           ,'cmd.origin("'+s+'")'            ],
      [ 1, 'orient'           ,'cmd.orient("'+s+'")'            ],
      [ 0, ''             , ''                      ],
      [ 1, 'label'      , mol_labels(s) ],
      [ 0, '', '' ],
      [ 1, 'enable'         ,'cmd.enable("'+s+'")'            ],
      [ 1, 'disable'        ,'cmd.disable("'+s+'")'            ],
      ]

def enable_disable(enable):
   if enable:
      result = [[ 2, 'Enable', '' ]]
      cmmd = 'cmd.enable("'
   else:
      result = [[ 2, 'Disable', '']]
      cmmd = 'cmd.disable("'
   return result + map(lambda ob,cm=cmmd:[1,ob,cm+ob+'")'],['all']+cmd.get_names('objects'))
   
def main_menu():
   return [
      [ 2, 'Main Pop-Up'  , '' ],
      [ 1, 'zoom (vis)'           ,'cmd.zoom("visible")'            ],
      [ 1, 'center (vis)'           ,'cmd.center("visible")'            ],      
      [ 1, 'orient (vis)'           ,'cmd.orient("visible")'            ],
      [ 1, 'reset'           ,'cmd.reset()'            ],
      [ 0, ''             , ''                      ],           
      [ 1, '(all)'      , all_option("all") ],
      [ 1, '(visible)'      , all_option("visible") ],
      [ 0, ''             , ''                      ],
      [ 1, 'ray'           ,'cmd.ray()' ],
      [ 0, ''             , ''                      ],
      [ 1, 'enable', enable_disable(1) ],
      [ 1, 'disable', enable_disable(0) ],   
      [ 0, ''             , ''                      ],
      [ 1, 'delete all'           ,'cmd.delete("all")' ],
      [ 1, 'reinitialize'           ,'cmd.reinitialize()' ],
      [ 1, 'quit'           ,'cmd.quit()' ],
      ]

def pick_option(title,s,object=0):
   result = [
      [ 2, title, '' ],
      [ 1, 'color'      , mol_color(s) ],
      [ 1, 'show'      , mol_show(s) ],
      [ 1, 'hide'      , mol_hide(s) ],
      [ 1, 'preset'  , presets(s)       ],      
      [ 0, ''             , ''                      ],
      [ 1, 'zoom'           ,'cmd.zoom("'+s+'")'            ],
      [ 1, 'center'           ,'cmd.center("'+s+'")'            ],
      [ 1, 'origin'           ,'cmd.origin("'+s+'")'            ],
      [ 1, 'orient'           ,'cmd.orient("'+s+'")'            ],
      [ 1, 'indicate'        ,'cmd.indicate("'+s+'")'            ],
      [ 0, ''             , ''                      ],
      [ 1, 'labels'      , mol_labels(s) ],
      ]
   
   if object:
      result.extend([
         [ 0, ''             , ''                      ],         
         [ 1, 'disable'        ,'cmd.disable("'+s+'")'            ],
         [ 0, ''             , ''                      ],
         [ 1, 'delete'        ,'cmd.delete("'+s+'")'            ]         
         ])
   else:
      result.extend([
      [ 0, ''             , ''                      ],      
      [ 1, 'create object','cmd.create(None,"'+s+'")'            ],      
      [ 0, ''             , ''                      ],
      [ 1, 'remove'      , 'cmd.remove("'+s+'")' ],     
                    ])
   return result

def pick_option_rev(title,s,object=0):
   result = pick_option(title,s,object)[1:]
   result.reverse()
   return result

def pick_menu(s1,s2):
   if s1[-1]=='`':
      title = s1[0:-1]
   else:
      title = s1
   return [[ 2, title     , '' ],
           [ 1, 'atom'    , pick_option("Atom",s2) ],
           [ 1, 'residue' , pick_option("Residue","(byres ("+s2+"))") ],
           [ 1, 'chain'   , pick_option("Chain","(bychain ("+s2+"))") ],
           [ 1, 'segment' , pick_option("Segment","(byseg ("+s2+"))") ],
           [ 1, 'object'  , pick_option("Object","(byobject ("+s2+"))",1) ],
           [ 0, ''             , ''                      ],
           [ 1, 'molecule', pick_option("Molecule","(bymol ("+s2+"))") ],
           [ 0, ''             , ''                      ],
           [ 1, 'fragment', pick_option("Fragment","(byfrag ("+s2+"))") ],
           [ 1, 'fragment+joint(s)', pick_option("Fragment","((byfrag ("+s2+")) extend 1)") ],
           ]
      
def seq_menu(s2,s3):
   
   return [[ 2, 'Sequence'    , '' ],
           [ 1, 'selection', pick_option('('+s3+')',s3) ],
           [ 0, ''             , ''                      ],
           [ 1, 'residue' , pick_option("Residue","(byres ("+s2+"))") ],
           [ 1, 'chain'   , pick_option("Chain","(bychain ("+s2+"))") ],
           [ 1, 'segment' , pick_option("Segment","(byseg ("+s2+"))") ],
           [ 1, 'object'  , pick_option("Object","(byobject ("+s2+"))",1) ],
           [ 0, ''             , ''                      ],
           [ 1, 'molecule', pick_option("Molecule","(bymol ("+s2+"))") ],
           [ 0, ''             , ''                      ],
           [ 1, 'C-alpha'    , pick_option("C-alpha",s2) ],
           ]
      

def seq_option(title,s,object=0):
   c=len(title)-1
   while title[c]!='/':
      c = c-1
   title = title[0:c+1]
   
   result = [
      [ 2, title, '' ],
      [ 1, 'color'      , mol_color(s) ],
      [ 1, 'show'      , mol_show(s) ],
      [ 1, 'hide'      , mol_hide(s) ],
      [ 1, 'preset'  , presets(s)       ],      
      [ 0, ''             , ''                      ],
      [ 1, 'zoom'           ,'cmd.zoom("'+s+'")'            ],
      [ 1, 'center'           ,'cmd.center("'+s+'")'            ],
      [ 1, 'origin'           ,'cmd.origin("'+s+'")'            ],
      [ 1, 'orient'           ,'cmd.orient("'+s+'")'            ],
      [ 1, 'indicate'        ,'cmd.indicate("'+s+'")'            ],
      [ 0, ''             , ''                      ],
      [ 1, 'labels'      , mol_labels(s) ],
      ]
   
   if object:
      result.extend([
         [ 0, ''             , ''                      ],         
         [ 1, 'disable'        ,'cmd.disable("'+s+'")'            ],
         [ 0, ''             , ''                      ],
         [ 1, 'delete'        ,'cmd.delete("'+s+'")'            ]         
         ])
   else:
      result.extend([
      [ 0, ''             , ''                      ],      
      [ 1, 'create object','cmd.create(None,"'+s+'")'            ],            
      [ 0, ''             , ''                      ],
      [ 1, 'remove'      , 'cmd.remove("'+s+'")' ],     
                    ])
   return result
   
