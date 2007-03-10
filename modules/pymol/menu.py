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

def extract(s):
    return [[ 2, 'Extract', '' ],
            [ 1, 'object', 'cmd.create(None,"'+s+'",extract="'+s+'",zoom=0)' ],
            [ 1, 'extend 1', 'cmd.create(None,"('+s+') extend 1",extract="'+s+'",zoom=0)' ],
            [ 1, 'byres extend 1', 'cmd.create(None,"byres (('+s+') extend 1)",extract="'+s+'",zoom=0)' ],            
            ]
            
def all_motion(s):
    return [[ 2, 'Camera Motions:'     , ''                       ],     
              [ 1, 'store'         , 'cmd.mview("store")'      ],
              [ 1, 'clear'       ,   'cmd.mview("clear")'      ],
              [ 0, ''               ,''                             ],
              [ 1, 'interpolate'   , 'cmd.mview("interpolate")'   ],
              [ 1, 'reinterpolate'   , 'cmd.mview("reinterpolate")'   ],            
              ]

def mol_motion(s):
    return [[ 2, 'Object Motions:'     , ''                       ],     
              [ 1, 'store'         , 'cmd.mview("store",object="'+s+'")'      ],
              [ 1, 'clear'       ,   'cmd.mview("clear",object="'+s+'")'    ],
              [ 0, ''               ,''                             ],
              [ 1, 'interpolate'   ,   'cmd.mview("interpolate",object="'+s+'")'    ],
              [ 1, 'reinterpolate'   ,   'cmd.mview("reinterpolate",object="'+s+'")'    ],
              ]

def rep_action(action,s) :
    return [
        [ 1, 'lines'      , 'cmd.'+action+'("lines"     ,"'+s+'")' ],
        [ 1, 'sticks'     , 'cmd.'+action+'("sticks"    ,"'+s+'")' ],
        [ 1, 'ribbon'     , 'cmd.'+action+'("ribbon"    ,"'+s+'")' ],
        [ 1, 'cartoon'    , 'cmd.'+action+'("cartoon"   ,"'+s+'")' ],
        [ 0, ''           , ''                               ],
        [ 1, 'label'      , 'cmd.'+action+'("labels"    ,"'+s+'")' ],
        [ 1, 'cell'       , 'cmd.'+action+'("cell"      ,"'+s+'")' ],
        [ 0, ''           , ''                               ],
        [ 1, 'nonbonded'  , 'cmd.'+action+'("nonbonded" ,"'+s+'")' ],
        [ 1, 'dots'       , 'cmd.'+action+'("dots"      ,"'+s+'")' ],
        [ 1, 'spheres'    , 'cmd.'+action+'("spheres"   ,"'+s+'")' ],
        [ 1, 'nb_spheres' , 'cmd.'+action+'("nb_spheres","'+s+'")' ],
        [ 0, ''           , ''                               ],
        [ 1, 'mesh'       , 'cmd.'+action+'("mesh"      ,"'+s+'")' ],
        [ 1, 'surface'    , 'cmd.'+action+'("surface"   ,"'+s+'")' ],
        ]

def mol_as(s):
    return (
        [[ 2, 'As:'   , '']]
        +rep_action('show_as',s)
        )

def mol_toggle(s):
    return (
        [[ 2, 'As:'   , '']]
        +rep_action('toggle',s)
        )

def show_misc(s):
    return [[ 2, 'Show:', '' ],
              [ 1, 'lines', 'cmd.show("lines","'+s+'")'],
              [ 1, 'sticks', 'cmd.show("sticks","'+s+'")'],
              [ 1, 'spheres', 'cmd.show("spheres","'+s+'")'],
              ]

def mol_show(s):
    return (
        [[ 2, 'Show:'      , ''                               ],
         [ 1, 'as'    , mol_as(s) ],
         [ 0, '', '']]
        + rep_action('show',s) +
        [[ 0, '', ''],
         [ 1, 'organic' , show_misc('(organic and ('+s+'))') ],
         [ 1, 'main chain' , show_misc("((byres ("+s+"))&n;ca,c,n,o,h)") ],
         [ 1, 'side chain' , show_misc("((byres ("+s+"))&(!(n;c,o,h|(n. n&!r. pro))))") ],
         [ 1, 'disulfides' , show_misc("(byres ((("+s+
            ") & r. CYS+CYX & n. SG) & bound_to (("+s+") & r. CYS+CYX & n. SG))) & n. CA+CB+SG") ]
         ])

    cmd.show("lines","(byres (("+s+" & r. CYS+CYX & n. SG) & bound_to ("+s+" & r. CYS+CYX & n. SG))) & n. CA+CB+SG")
    
def mol_hide(s):
    return (
        [[ 2, 'Hide:'     , ''                                ],
         [ 1, 'everything', 'cmd.hide("everything","'+s+'")'  ],
         [ 0, ''          , ''                                ]]
        + rep_action('hide',s) +
        [[ 0, ''          , ''                                ],
         [ 1, 'main chain', 'cmd.hide("((byres ('+s+'))&(n. c,o,h|(n. n&!r. pro)))")' ],
         [ 1, 'side chain', 'cmd.hide("((byres ('+s+'))&!(n. ca,c,o,h|(n. n&!r. pro)))")' ],
         [ 1, 'waters'    , 'cmd.hide("(solvent and ('+s+'))")'     ],                      
         [ 0, ''          , ''                                ],
         [ 1, 'hydrogens' , 'cmd.hide("('+s+' and hydro)")'   ],
         [ 0, ''          , ''                                ],           
         [ 1, 'unselected', 'cmd.hide("(not '+s+')")'         ],
         ]
        )
        
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

def map_show(s):
    return [[ 2, 'Show:'       , ''                             ],
              [ 1, 'dots'        , 'cmd.show("dots","'+s+'")'     ],           
              [ 1, 'extent'        , 'cmd.show("extent","'+s+'")'     ],
              [ 1, 'everything'  , 'cmd.show("everything","'+s+'")'          ]]

def map_hide(s):
    return [[ 2, 'Hide:'     ,''                                ],
              [ 1, 'dots'        , 'cmd.hide("dots","'+s+'")'     ],           
              [ 1, 'extent'      , 'cmd.hide("extent","'+s+'")'     ],
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

def surface_show(s):
    return [[ 2, 'Show:'       , ''                             ],
              [ 1, 'surface'        , 'cmd.show("surface","'+s+'")'     ],           
              [ 1, 'cell'        , 'cmd.show("cell","'+s+'")'     ],
              [ 1, 'everything'  , 'cmd.enable("'+s+'")'          ]]

def surface_hide(s):
    return [[ 2, 'Hide:'       , ''                             ],
              [ 1, 'surface'        , 'cmd.hide("surface","'+s+'")'     ],                      
              [ 1, 'cell'        , 'cmd.hide("cell","'+s+'")'      ],           
              [ 1, 'everything'  , 'cmd.disable("'+s+'")'          ]]

def slice_show(s):
    return [[ 2, 'Show:'       , ''                             ],
              [ 1, 'slice'       , 'cmd.show("slice","'+s+'")'     ],
              ]

def slice_hide(s):
    return [[ 2, 'Hide:'       , ''                             ],
              [ 1, 'slice'        , 'cmd.hide("slice","'+s+'")'     ],
              ]

def by_elem2(s):
    return [
        [ 2, 'Atoms'     ,''                               ],
        [1,'\\494C\\777H\\229N\\922O\\950S...','util.cba(10,"'+s+'")'],# lime
        [1,'\\155C\\777H\\229N\\922O\\950S...','util.cba(5262,"'+s+'")'],# deepteal
        [1,'\\904C\\777H\\229N\\922O\\950S...','util.cba(12,"'+s+'")'],# hotpink
        [1,'\\983C\\777H\\229N\\922O\\950S...','util.cba(36,"'+s+'")'],# yelloworange
        [1,'\\525C\\777H\\229N\\922O\\950S...','util.cba(5271,"'+s+'")'],# violetpurple
        [1,'\\666C\\777H\\229N\\922O\\950S...','util.cba(124,"'+s+'")'],# grey70
        [1,'\\049C\\777H\\229N\\922O\\950S...','util.cba(17,"'+s+'")'],# marine
        [1,'\\760C\\777H\\229N\\922O\\950S...','util.cba(18,"'+s+'")'],# olive
        ]
        
def by_elem3(s):
    return [
        [ 2, 'Atoms'     ,''                               ],
        [1,'\\564C\\777H\\229N\\922O\\950S...','util.cba(5270,"'+s+'")'],# smudge
        [1,'\\077C\\777H\\229N\\922O\\950S...','util.cba(20,"'+s+'")'],# teal
        [1,'\\644C\\777H\\229N\\922O\\950S...','util.cba(5272,"'+s+'")'],# dirtyviolet
        [1,'\\976C\\777H\\229N\\922O\\950S...','util.cba(52,"'+s+'")'],# wheat
        [1,'\\944C\\777H\\229N\\922O\\950S...','util.cba(5258,"'+s+'")'],# deepsalmon
        [1,'\\978C\\777H\\229N\\922O\\950S...','util.cba(5274,"'+s+'")'],# lightpink
        [1,'\\499C\\777H\\229N\\922O\\950S...','util.cba(5257,"'+s+'")'],# aquamarine
        [1,'\\994C\\777H\\229N\\922O\\950S...','util.cba(5256,"'+s+'")'],# paleyellow
        ]
    
def by_elem4(s):
    return [
        [ 2, 'Atoms'     ,''                               ],
        [1,'\\094C\\777H\\229N\\922O\\950S...','util.cba(15,"'+s+'")'],# limegreen
        [1,'\\247C\\777H\\229N\\922O\\950S...','util.cba(5277,"'+s+'")'],# skyblue
        [1,'\\824C\\777H\\229N\\922O\\950S...','util.cba(5279,"'+s+'")'],# warmpink
        [1,'\\792C\\777H\\229N\\922O\\950S...','util.cba(5276,"'+s+'")'],# limon
        [1,'\\949C\\777H\\229N\\922O\\950S...','util.cba(53,"'+s+'")'],# violet
        [1,'\\889C\\777H\\229N\\922O\\950S...','util.cba(5278,"'+s+'")'],# bluewhite
        [1,'\\297C\\777H\\229N\\922O\\950S...','util.cba(5275,"'+s+'")'],# greencyan
        [1,'\\653C\\777H\\229N\\922O\\950S...','util.cba(5269,"'+s+'")'],# sand
        ]

def by_elem5(s):
    return [
        [ 2, 'Atoms'     ,''                               ],
[1,'\\252C\\777H\\229N\\922O\\950S...','util.cba(22,"'+s+'")'],# forest
[1,'\\466C\\777H\\229N\\922O\\950S...','util.cba(5266,"'+s+'")'],# lightteal
[1,'\\755C\\777H\\229N\\922O\\950S...','util.cba(5280,"'+s+'")'],# darksalmon
[1,'\\570C\\777H\\229N\\922O\\950S...','util.cba(5267,"'+s+'")'],# splitpea
[1,'\\634C\\777H\\229N\\922O\\950S...','util.cba(5268,"'+s+'")'],# raspberry
[1,'\\555C\\777H\\229N\\922O\\950S...','util.cba(104,"'+s+'")'],# grey50
[1,'\\226C\\777H\\229N\\922O\\950S...','util.cba(23,"'+s+'")'],# deepblue
[1,'\\632C\\777H\\229N\\922O\\950S...','util.cba(51,"'+s+'")'],# brown
              ]
    
def by_elem(s):
    return [
        [ 2, 'Atoms'     ,''                               ],
        [1,' \\777H\\229N\\922O\\950S...','util.cnc("'+s+'")'],

[1,'\\292C\\777H\\229N\\922O\\950S...','util.cba(26,"'+s+'")'],# carbon
[1,'\\099C\\777H\\229N\\922O\\950S...','util.cba(5,"'+s+'")'],# cyan
[1,'\\927C\\777H\\229N\\922O\\950S...','util.cba(154,"'+s+'")'],# lightmagenta
[1,'\\990C\\777H\\229N\\922O\\950S...','util.cba(6,"'+s+'")'],# yellow
[1,'\\955C\\777H\\229N\\922O\\950S...','util.cba(9,"'+s+'")'],# salmon
[1,'\\888C\\777H\\229N\\922O\\950S...','util.cba(29,"'+s+'")'],# hydrogen
[1,'\\449C\\777H\\229N\\922O\\950S...','util.cba(11,"'+s+'")'],# slate
[1,'\\962C\\777H\\229N\\922O\\950S...','util.cba(13,"'+s+'")'],# orange
        [ 1, 'set 2'     ,by_elem2(s)                    ],
        [ 1, 'set 3'     ,by_elem3(s)                    ],
        [ 1, 'set 4'     ,by_elem4(s)                    ],
        [ 1, 'set 5'     ,by_elem5(s)                    ],      
              ]

def by_ss(s):
    return [
                [ 2, 'By Secondary Structure:'     ,''                               ],
    [ 1, '\\900Helix \\990Sheet \\090Loop'  , 'util.cbss("'+s+'","red","yellow","green")'],
    [ 1, '\\099Helix \\909Sheet \\955Loop'  , 'util.cbss("'+s+'","cyan","magenta","salmon")'],
    [ 1, '\\099Helix \\900Sheet \\909Loop'  , 'util.cbss("'+s+'","cyan","red","magenta")',],
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
              [ 1, 'b-factors'   , 'cmd.spectrum("b",selection=("'+s+'"),quiet=0)'         ],
              [ 1, 'b-factors(*/ca)'   , 'cmd.spectrum("b",selection="(('+s+')&*/ca)",quiet=0)'         ],                       
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


def reds(s):
    return [
        [ 2, 'Reds'     ,''                               ],
        [1,'\\900red','cmd.color(4,"'+s+'")'],
        [1,'\\922tv_red','cmd.color(32,"'+s+'")'],
        [1,'\\634raspberry','cmd.color(5268,"'+s+'")'],
        [1,'\\755darksalmon','cmd.color(5280,"'+s+'")'],
        [1,'\\955salmon','cmd.color(9,"'+s+'")'],
        [1,'\\944deepsalmon','cmd.color(5258,"'+s+'")'],
        [1,'\\824warmpink','cmd.color(5279,"'+s+'")'],
        [1,'\\611firebrick','cmd.color(49,"'+s+'")'],
        [1,'\\522ruby','cmd.color(21,"'+s+'")'],
        [1,'\\521chocolate','cmd.color(50,"'+s+'")'],
        [1,'\\632brown','cmd.color(51,"'+s+'")'],
        ]

def greens(s):
    return [
        [ 2, 'Greens'     ,''                               ],
        [1,'\\090green','cmd.color(3,"'+s+'")'],
        [1,'\\292tv_green','cmd.color(33,"'+s+'")'],
        [1,'\\490chartreuse','cmd.color(14,"'+s+'")'],
        [1,'\\570splitpea','cmd.color(5267,"'+s+'")'],
        [1,'\\564smudge','cmd.color(5270,"'+s+'")'],
        [1,'\\686palegreen','cmd.color(5259,"'+s+'")'],
        [1,'\\094limegreen','cmd.color(15,"'+s+'")'],
        [1,'\\494lime','cmd.color(10,"'+s+'")'],
        [1,'\\792limon','cmd.color(5276,"'+s+'")'],      
        [1,'\\252forest','cmd.color(22,"'+s+'")'],
        ]

def blues(s):
    return [
        [ 2, 'Blues'     ,''                               ],
        [1,'\\009blue','cmd.color(2,"'+s+'")'],
        [1,'\\339tv_blue','cmd.color(34,"'+s+'")'],
        [1,'\\049marine','cmd.color(17,"'+s+'")'],
        [1,'\\449slate','cmd.color(11,"'+s+'")'],
        [1,'\\779lightblue','cmd.color(5263,"'+s+'")'],
        [1,'\\247skyblue','cmd.color(5277,"'+s+'")'],
        [1,'\\409purpleblue','cmd.color(16,"'+s+'")'],
        [1,'\\226deepblue','cmd.color(23,"'+s+'")'],
        [1,'\\115density','cmd.color(4155,"'+s+'")'],
        ]

def yellows(s):
    return [
        [ 2, 'Yellows'     ,''                               ],
        [1,'\\990yellow','cmd.color(6,"'+s+'")'],
        [1,'\\992tv_yellow','cmd.color(35,"'+s+'")'],
        [1,'\\994paleyellow','cmd.color(5256,"'+s+'")'],
        [1,'\\983yelloworange','cmd.color(36,"'+s+'")'],            
        [1,'\\792limon','cmd.color(5276,"'+s+'")'],
        [1,'\\976wheat','cmd.color(52,"'+s+'")'],
        [1,'\\653sand','cmd.color(5269,"'+s+'")'],
        ]

def magentas(s):
    return [
        [ 2, 'Magentas'     ,''                               ],
        [1,'\\909magenta','cmd.color(8,"'+s+'")'],
        [1,'\\927lightmagenta','cmd.color(154,"'+s+'")'],
        [1,'\\904hotpink','cmd.color(12,"'+s+'")'],
        [1,'\\968pink','cmd.color(48,"'+s+'")'],
        [1,'\\978lightpink','cmd.color(5274,"'+s+'")'],
        [1,'\\644dirtyviolet','cmd.color(5272,"'+s+'")'],
        [1,'\\949violet','cmd.color(53,"'+s+'")'],
        [1,'\\525violetpurple','cmd.color(5271,"'+s+'")'],
        [1,'\\707purple','cmd.color(19,"'+s+'")'],
        [1,'\\515deeppurple','cmd.color(5261,"'+s+'")'],
        ]

def cyans(s):
    return [
        [ 2, 'Cyans'     ,''                               ],
        [1,'\\099cyan','cmd.color(5,"'+s+'")'],
        [1,'\\799palecyan','cmd.color(5265,"'+s+'")'],
        [1,'\\499aquamarine','cmd.color(5257,"'+s+'")'],
        [1,'\\297greencyan','cmd.color(5275,"'+s+'")'],
        [1,'\\077teal','cmd.color(20,"'+s+'")'],
        [1,'\\155deepteal','cmd.color(5262,"'+s+'")'],
        [1,'\\466lightteal','cmd.color(5266,"'+s+'")'],
        ]

def oranges(s):
    return [
        [ 2, 'Oranges'     ,''                               ],
        [1,'\\950orange','cmd.color(13,"'+s+'")'],
        [1,'\\951tv_orange','cmd.color(37,"'+s+'")'],
        [1,'\\962brightorange','cmd.color(30,"'+s+'")'],
        [1,'\\985lightorange','cmd.color(5264,"'+s+'")'],      
        [1,'\\983yelloworange','cmd.color(36,"'+s+'")'],      
        [1,'\\760olive','cmd.color(18,"'+s+'")'],
        [1,'\\551deepolive','cmd.color(5260,"'+s+'")'],
        ]

def tints(s):
    return [
        [ 2, 'Tints'     ,''                               ],
        [1,'\\976wheat','cmd.color(52,"'+s+'")'],
        [1,'\\686palegreen','cmd.color(5259,"'+s+'")'],
        [1,'\\779lightblue','cmd.color(5263,"'+s+'")'],      
        [1,'\\994paleyellow','cmd.color(5256,"'+s+'")'],
        [1,'\\978lightpink','cmd.color(5274,"'+s+'")'],
        [1,'\\799palecyan','cmd.color(5265,"'+s+'")'],
        [1,'\\985lightorange','cmd.color(5264,"'+s+'")'],            
        [1,'\\889bluewhite','cmd.color(5278,"'+s+'")'],
        ]
    
def grays(s):
    return [
        [ 2, 'Grays'     ,''                               ],
        [ 1, '\\999white ', 'cmd.color("white","'+s+'")'  ],
        [ 1, '\\999gray90 ', 'cmd.color("grey90","'+s+'")'  ],
        [ 1, '\\888gray80 ', 'cmd.color("grey80","'+s+'")'  ],
        [ 1, '\\777gray70 ', 'cmd.color("grey70","'+s+'")'  ],
        [ 1, '\\666gray60 ', 'cmd.color("grey60","'+s+'")'  ],
        [ 1, '\\555gray50 ', 'cmd.color("grey50","'+s+'")'  ],
        [ 1, '\\444gray40 ', 'cmd.color("grey40","'+s+'")'  ],
        [ 1, '\\333gray30 ', 'cmd.color("grey30","'+s+'")'  ],
        [ 1, '\\222gray20 ', 'cmd.color("grey20","'+s+'")'  ],
        [ 1, '\\222gray10 ', 'cmd.color("grey10","'+s+'")'  ],
        [ 1, '\\222black ', 'cmd.color("black","'+s+'")'  ],
        ]

def all_colors(s):
    return [
    [ 1, '\\900reds'        ,reds(s) ],
    [ 1, '\\090greens'      ,greens(s) ],
    [ 1, '\\009blues'       ,blues(s) ],
    [ 1, '\\990yellows'      ,yellows(s) ],
    [ 1, '\\909magentas'    , magentas(s) ],
    [ 1, '\\099cyans'        , cyans(s) ],
    [ 1, '\\950oranges'        , oranges(s) ],   
    [ 1, '\\978tints'        ,tints(s) ],
    [ 1, '\\666grays'        ,grays(s) ],
#   [ 0, '', ''],
#   [ 1, '\\900red'         ,'cmd.color("red","'+s+'")'  ],
#   [ 1, '\\090green'       ,'cmd.color("green","'+s+'")'  ],
#   [ 1, '\\009blue'        ,'cmd.color("blue","'+s+'")'  ],
#   [ 1, '\\990yellow'      ,'cmd.color("yellow","'+s+'")'  ],
#   [ 1, '\\909magenta' ,'cmd.color("magenta","'+s+'")'  ],
#   [ 1, '\\099cyan'  ,'cmd.color("cyan","'+s+'")'  ],           
#   [ 1, '\\955salmon'      ,'cmd.color("salmon","'+s+'")'  ],
#   [1,  '\\940orange','cmd.color(13,"'+s+'")'],
#   
#   [ 1, '\\555gray'    ,'cmd.color("gray","'+s+'")'  ],
#   [ 1, '\\999white'       ,'cmd.color("white","'+s+'")'  ],
    
        ]

def color_auto(s):
    return [
        [ 2, 'Auto'     ,''                               ],
        [ 1, 'elem c', 'cmd.color("auto","('+s+') and elem c")' ],
        [ 0, ''                                , ''                 ],
        [ 1, 'all','cmd.color("auto","'+s+'")' ],                  
        [ 0, ''                                , ''                 ],
        [ 1, '\\900b\\950y \\090o\\099b\\059j\\999(e. c)',
          'util.color_objs("('+s+' and elem c)")'],
        [ 1, '\\900b\\950y \\090o\\099b\\059j',
          'util.color_objs("('+s+')")'],
        ]
   
def mol_color(s):
    return (
        [[ 2, 'Color:'     ,''                               ],
         [ 1, 'by element'  , by_elem(s) ],
         [ 1, 'by chain' , by_chain(s) ],
         [ 1, 'by ss  '  , by_ss(s) ],
         [ 1, '\\900s\\950p\\990e\\090c\\099t\\059r\\009u\\555m', spectrum(s) ],
         [ 0, ''                                , ''                 ],
         [ 1, 'auto', color_auto(s) ],
         [ 0, ''                                , ''                 ],         
         ] +
        all_colors(s))

def general_color(s):
    return [[ 2, 'Color:'     ,''                        ]] + all_colors(s)

def preset_ligand_sites(s):
    return [[ 2, 'Ligand Sites:', ''],
              [ 1, 'cartoon'   , 'preset.ligand_cartoon("'+s+'")'          ],
              [ 0, '', ''],
              [ 1, 'solid surface'   , 'preset.ligand_sites("'+s+'")'          ],
              [ 1, 'solid (better)'   , 'preset.ligand_sites_hq("'+s+'")'          ],
              [ 0, '', ''],
              [ 1, 'transparent surface'   , 'preset.ligand_sites_trans("'+s+'")'          ],
              [ 1, 'transparent (better)'   , 'preset.ligand_sites_trans_hq("'+s+'")'          ],
              [ 0, '', ''],
              [ 1, 'dot surface'   , 'preset.ligand_sites_dots("'+s+'")'          ],
              [ 0, '', ''],
              [ 1, 'mesh surface'   , 'preset.ligand_sites_mesh("'+s+'")'          ]]

def presets(s):
    return [[ 2, 'Preset:'       ,''                        ],     
              [ 1, 'simple'   ,'preset.simple("'+s+'")'          ],
              [ 1, 'simple (no solvent)'   ,'preset.simple_no_solv("'+s+'")'          ],           
              [ 1, 'ball and stick' , 'preset.ball_and_stick("'+s+'")' ],
              [ 1, 'b factor putty' , 'preset.b_factor_putty("'+s+'")' ],
              [ 1, 'technical'   , 'preset.technical("'+s+'")'          ],
              [ 1, 'ligands'   , 'preset.ligands("'+s+'")'          ],
              [ 1, 'ligand sites'   , preset_ligand_sites(s)         ],
              [ 1, 'pretty ', 'preset.pretty("'+s+'")'          ],
              [ 1, 'pretty (with solvent)'     , 'preset.pretty_solv("'+s+'")'          ],
              [ 1, 'publication '   , 'preset.publication("'+s+'")'          ],
              [ 1, 'publication (with solvent)'   , 'preset.pub_solv("'+s+'")'          ],
              [ 0, ''               ,''                             ],                      
              [ 1, 'default'   ,'preset.default("'+s+'")'          ],           
              ]

def hydrogens(s):
   return [[ 2, 'Hydrogens:'       ,''                        ],     
           [ 1, 'add'   ,'cmd.h_add("'+s+'")'          ],           
           [ 1, 'remove'   ,'cmd.remove("('+s+') and hydro")'          ],
           ]

def state(s):
    return [[ 2, 'State:'       ,''                        ],
              [ 1, 'freeze'  ,'cmd.set("state",cmd.get_state(),"'+s+'")'        ],
              [ 1, 'thaw'  ,'cmd.set("state",cmd.get("state","'+s+'"));cmd.unset("state","'+s+'")'        ],           
              ]
    
def movement(s):
    return [[ 2, 'Movement:'       ,''                        ],     
              [ 1, 'protect'   ,'cmd.protect("'+s+'")'          ],
              [ 1, 'deprotect'   ,'cmd.deprotect("'+s+'")'          ],           
              ]

def sequence(s):
    return [[ 2, 'Sequence:'       ,''                        ],     
              [ 1, 'include'   ,'cmd.set("seq_view","on","'+s+'")'          ],
              [ 1, 'exclude'   ,'cmd.set("seq_view","off","'+s+'")'          ],
              [ 0, ''               ,''                             ],                      
              [ 1, 'default'   ,'cmd.unset("seq_view","'+s+'")'          ],                      
              ]

def masking(s):
    return [[ 2, 'Masking:'       ,''                        ],     
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
#              [ 2, '\\955WARNING:\\595 Unvalidated and experimental code!', '' ],
              [ 1, 'protein contact potential (local)', 'util.protein_vacuum_esp("'+s+'",mode=2,quiet=0)'          ],
#           [ 1, 'protein surface potential (absolute)', 'util.protein_vacuum_esp("'+s+'",mode=0,quiet=0)'          ],
#           [ 1, 'protein surface potential (relative)', 'util.protein_vacuum_esp("'+s+'",mode=1,quiet=0)'          ],
              [ 2, '\\955NOTE:\\559 Due to short cutoffs, truncations, and', ''],
              [ 2, '\\559lack of solvent "screening", these computed ', ''],
              [ 2, '\\559potentials are only qualitatively useful.', ''],
              [ 2, '\\559Please view with skepticism!', '' ],
              ]

def symmetry(s):
    return [[ 2, 'Symmetry Mates:'       ,''                        ],
              [ 2, '\\955 +/- one unit cell and...', '' ],
              [ 1, 'within 4 A', 'cmd.symexp("'+s+'_","'+s+'","'+s+'",cutoff=4,segi=1)'          ],
              [ 1, 'within 5 A', 'cmd.symexp("'+s+'_","'+s+'","'+s+'",cutoff=5,segi=1)'          ],
              [ 1, 'within 6 A', 'cmd.symexp("'+s+'_","'+s+'","'+s+'",cutoff=6,segi=1)'          ],
              [ 1, 'within 8 A', 'cmd.symexp("'+s+'_","'+s+'","'+s+'",cutoff=8,segi=1)'          ],
              [ 1, 'within 12 A', 'cmd.symexp("'+s+'_","'+s+'","'+s+'",cutoff=12,segi=1)'          ],
              [ 1, 'within 20 A', 'cmd.symexp("'+s+'_","'+s+'","'+s+'",cutoff=20,segi=1)'          ],
              [ 1, 'within 50 A', 'cmd.symexp("'+s+'_","'+s+'","'+s+'",cutoff=50,segi=1)'          ],
              [ 1, 'within 100 A', 'cmd.symexp("'+s+'_","'+s+'","'+s+'",cutoff=100,segi=1)'          ],
              [ 1, 'within 250 A', 'cmd.symexp("'+s+'_","'+s+'","'+s+'",cutoff=250,segi=1)'          ],
              [ 1, 'within 1000 A', 'cmd.symexp("'+s+'_","'+s+'","'+s+'",cutoff=1000,segi=1)'          ]]

def mol_assign(s):
    return [[ 2, 'Assign:'       ,''                        ],     
              [ 1, 'Amber 99 atomic properties',  'util.assign_amber99("'+s+'")' ],
              ]

def selection(s):
    return [[ 2, 'Selections:', '' ],
              [ 1, 'all', 'cmd.select("'+s+'_all","'+s+'")'],
              [ 1, 'polymer', 'cmd.select("'+s+'_polymer","('+s+') and polymer")'],
              [ 1, 'organic', 'cmd.select("'+s+'_organic","('+s+') and organic")'],
              [ 1, 'solvent', 'cmd.select("'+s+'_solvent","('+s+') and solvent")'],           
              ]

def mol_generate(s):
    return [[ 2, 'Generate:'       ,''                        ],
              [ 1, 'selection', selection(s) ],
              [ 1, 'symmetry mates', symmetry(s) ],
              [ 1, 'vacuum electrostatics', vacuum(s) ],
#           [ 1, 'assign', mol_assign(s) ],
              ]
    
def invert(s):
    return [[ 2, 'Invert:'       ,''                        ],     
              [ 1, 'within object(s)'     ,'cmd.select("'+s+'","((byobj '+s+') and not '+s+')",enable=1)'    ],
              [ 1, 'within segments(s)'     ,'cmd.select("'+s+'","((byseg '+s+') and not '+s+')",enable=1)'    ],           
              [ 1, 'within chains(s)'     ,'cmd.select("'+s+'","((bychain '+s+') and not '+s+')",enable=1)'    ],
              [ 1, 'within residues(s)'   ,'cmd.select("'+s+'","((byres '+s+') and not '+s+')",enable=1)'    ],                      
              [ 0, ''               ,''                             ],
              [ 1, 'within molecule(s)'     ,'cmd.select("'+s+'","((bymol '+s+') and not '+s+')",enable=1)'    ],
              [ 0, ''               ,''                             ],
              [ 1, 'within any'     ,'cmd.select("'+s+'","(not '+s+')",enable=1)'    ],
              ]

def complete(s):
    return [[ 2, 'Complete:'       ,''                        ],     

              [ 1, 'resides'  ,'cmd.select("'+s+'","(byres '+s+')",enable=1)'      ],
              [ 1, 'chains'  ,'cmd.select("'+s+'","(bychain '+s+')",enable=1)'      ],
              [ 1, 'segments'  ,'cmd.select("'+s+'","(byseg '+s+')",enable=1)'      ],
              [ 1, 'objects'  ,'cmd.select("'+s+'","(byobj '+s+')",enable=1)'      ],
              [ 0, ''               ,''                             ],           
              [ 1, 'molecules'  ,'cmd.select("'+s+'","(bymol '+s+')",enable=1)'      ],
              [ 0, ''               ,''                             ],
              [ 1, 'C-alphas'  ,'cmd.select("'+s+'","(bycalpha '+s+')",enable=1)'      ],           
              ]

def modify_by_object(s,op):
    list = cmd.get_names("public_objects",1)[0:25] # keep this practical
    list = filter(lambda x:cmd.get_type(x)=="object:molecule",list)
    result = [[ 2, 'Object:', '']]
    for a in list:
        if a!=s:
            result.append([1,a,
                                'cmd.select("'+s+'","('+s+') '+op+' ('+a+')",enable=1)'])
    return result

def modify_by_sele(s,op):
    list = cmd.get_names("public_selections",0)[0:25] # keep this practical
    result = [[ 2, 'Selection:', '']]
    for a in list:
        if a!=s:
            result.append([1,a,
                                'cmd.select("'+s+'","('+s+') '+op+' ('+a+')",enable=1)'])
    return result

def restrict(s):
    return [[ 2, 'Restrict:'       ,''                        ],     
            [ 1, 'to object'   , modify_by_object(s,'and') ],
            [ 1, 'to selection' , modify_by_sele(s,'and') ],
            [ 0, ''               ,''                             ],           
            [ 1, 'to visible'   , 'cmd.select("'+s+'","('+s+') and vis",enable=1)'],
            [ 0, ''               ,''                             ],
            [ 1, 'to polymer'   , 'cmd.select("'+s+'","('+s+') and polymer",enable=1)'],
            [ 1, 'to solvent'   , 'cmd.select("'+s+'","('+s+') and solvent",enable=1)'],
            [ 1, 'to organic'   , 'cmd.select("'+s+'","('+s+') and organic",enable=1)'],
            [ 1, 'to inorganic'   , 'cmd.select("'+s+'","('+s+') and inorganic",enable=1)'],
            ]

def include(s):
    return [[ 2, 'Include:'       ,''                        ],     
              [ 1, 'object'   , modify_by_object(s,'or') ],
              [ 1, 'selection' , modify_by_sele(s,'or') ],
              [ 0, ''               ,''                             ],           
              [ 1, 'visible'   , 'cmd.select("'+s+'","('+s+') or vis",enable=1)'],
              ]

def exclude(s):
    return [[ 2, 'Exclude:'       ,''                        ],     
            [ 1, 'object'   , modify_by_object(s,'and not') ],
            [ 1, 'selection' , modify_by_sele(s,'and not') ],
            [ 0, ''               ,''                             ],           
            [ 1, 'polymer'   , 'cmd.select("'+s+'","('+s+') and not organic",enable=1)'],
            [ 1, 'solvent'   , 'cmd.select("'+s+'","('+s+') and not solvent",enable=1)'],
            [ 1, 'organic'   , 'cmd.select("'+s+'","('+s+') and not organic",enable=1)'],
            [ 1, 'inorganic' , 'cmd.select("'+s+'","('+s+') and not organic",enable=1)'],
            ]

def expand(s):
    return [[ 2, 'Expand:'       ,''                        ],     
              [ 1, 'by 4 A'  ,'cmd.select("'+s+'","('+s+' expand 4)",enable=1)' ],
              [ 1, 'by 5 A'  ,'cmd.select("'+s+'","('+s+' expand 5)",enable=1)' ],
              [ 1, 'by 6 A'  ,'cmd.select("'+s+'","('+s+' expand 6)",enable=1)' ],           
              [ 1, 'by 8 A'  ,'cmd.select("'+s+'","('+s+' expand 8)",enable=1)' ],
              [ 1, 'by 12 A'  ,'cmd.select("'+s+'","('+s+' expand 12)",enable=1)' ],
              [ 1, 'by 20 A'  ,'cmd.select("'+s+'","('+s+' expand 20)",enable=1)' ],           
              [ 0, ''               ,''                             ],
              [ 1, 'by 4 A, residues'  ,'cmd.select("'+s+'","(byres ('+s+' expand 4))",enable=1)' ],
              [ 1, 'by 5 A, residues'  ,'cmd.select("'+s+'","(byres ('+s+' expand 5))",enable=1)' ],
              [ 1, 'by 6 A, residues'  ,'cmd.select("'+s+'","(byres ('+s+' expand 6))",enable=1)' ],
              [ 1, 'by 8 A, residues'   ,'cmd.select("'+s+'","(byres ('+s+' expand 8))",enable=1)' ],
              [ 1, 'by 12 A, residues'  ,'cmd.select("'+s+'","(byres ('+s+' expand 12))",enable=1)' ],
              [ 1, 'by 20 A, residues'  ,'cmd.select("'+s+'","(byres ('+s+' expand 20))",enable=1)' ],
              ]

def around(s):
    return [[ 2, 'Around:'       ,''                        ],     
              [ 1, 'atoms within 4 A'  ,'cmd.select("'+s+'","('+s+' around 4)",enable=1)' ],
              [ 1, 'atoms within 5 A'  ,'cmd.select("'+s+'","('+s+' around 5)",enable=1)' ],           
              [ 1, 'atoms within 6 A'  ,'cmd.select("'+s+'","('+s+' around 6)",enable=1)' ],           
              [ 1, 'atoms within 8 A'  ,'cmd.select("'+s+'","('+s+' around 8)",enable=1)' ],
              [ 1, 'atoms within 12 A'  ,'cmd.select("'+s+'","('+s+' around 12)",enable=1)' ],
              [ 1, 'atoms within 20 A'  ,'cmd.select("'+s+'","('+s+' around 20)",enable=1)' ],
              [ 0, ''               ,''                             ],           
              [ 1, 'residues within 4 A'  ,'cmd.select("'+s+'","(byres ('+s+' around 4))",enable=1)' ],
              [ 1, 'residues within 5 A'  ,'cmd.select("'+s+'","(byres ('+s+' around 5))",enable=1)' ],           
              [ 1, 'residues within 6 A'  ,'cmd.select("'+s+'","(byres ('+s+' around 6))",enable=1)' ],
              [ 1, 'residues within 8 A'  ,'cmd.select("'+s+'","(byres ('+s+' around 8))",enable=1)' ],
              [ 1, 'residues within 12 A'  ,'cmd.select("'+s+'","(byres ('+s+' around 12))",enable=1)' ],
              [ 1, 'residues within 20 A'  ,'cmd.select("'+s+'","(byres ('+s+' around 20))",enable=1)' ],                                 
              ]
    
def extend(s):
    return [[ 2, 'Extend:'       ,''                        ],     
              [ 1, 'by 1 bond'  ,'cmd.select("'+s+'","('+s+' extend 1)",enable=1)' ],
              [ 1, 'by 2 bonds'  ,'cmd.select("'+s+'","('+s+' extend 2)",enable=1)' ],           
              [ 1, 'by 3 bonds'  ,'cmd.select("'+s+'","('+s+' extend 3)",enable=1)' ],
              [ 1, 'by 4 bonds'  ,'cmd.select("'+s+'","('+s+' extend 4)",enable=1)' ],
              [ 1, 'by 5 bonds'  ,'cmd.select("'+s+'","('+s+' extend 5)",enable=1)' ],
              [ 1, 'by 6 bonds'  ,'cmd.select("'+s+'","('+s+' extend 6)",enable=1)' ],           
              [ 0, ''               ,''                             ],
              [ 1, 'by 1 bond, residues'  ,'cmd.select("'+s+'","(byres ('+s+' extend 1))",enable=1)' ],
              [ 1, 'by 2 bonds, residues'  ,'cmd.select("'+s+'","(byres ('+s+' extend 2))",enable=1)' ],
              [ 1, 'by 3 bonds, residues'   ,'cmd.select("'+s+'","(byres ('+s+' extend 3))",enable=1)' ],
              [ 1, 'by 4 bonds, residues'  ,'cmd.select("'+s+'","(byres ('+s+' extend 4))",enable=1)' ],
              [ 1, 'by 5 bonds, residues'  ,'cmd.select("'+s+'","(byres ('+s+' extend 5))",enable=1)' ],
              [ 1, 'by 6 bonds, residues'  ,'cmd.select("'+s+'","(byres ('+s+' extend 6))",enable=1)' ],
              ]

def polar(s):
    return [[ 2, 'Polar Contacts:', ''],
              [ 1, 'within selection'  ,
                 'cmd.dist("'+s+'_polar_conts","'+s+'","'+s+'",quiet=1,mode=2,label=0,reset=1);cmd.enable("'+s+'_polar_conts")'],
              [ 1, 'involving side chains'  ,
                 'cmd.dist("'+s+'_polar_conts","('+s+')","('+s+
                 ') and polymer and not (name n,o,h)",quiet=1,mode=2,label=0,reset=1);cmd.enable("'+s+'_polar_conts")'],
              [ 1, 'involving solvent'  ,
                 'cmd.dist("'+s+'_polar_conts","('+s+') and solvent","('+s+
                 ') and not (solvent)",quiet=1,mode=2,label=0,reset=1);cmd.enable("'+s+'_polar_conts")'],
              [ 1, 'excluding solvent'  ,
                 'cmd.dist("'+s+'_polar_conts","('+s+') and not (solvent)","('+s+
                 ') and not (solvent)",quiet=1,mode=2,label=0,reset=1);cmd.enable("'+s+'_polar_conts")'],
              [ 1, 'excluding main chain'  ,
                 'cmd.dist("'+s+'_polar_conts","('+s+') and not (polymer and name n,o,h)","('+s+
                 ') and not (polymer and name n,o,h)",quiet=1,mode=2,label=0,reset=1);cmd.enable("'+s+'_polar_conts")'],
              [ 1, 'excluding intra-main chain'  ,
                 'cmd.dist("'+s+'_polar_conts","('+s+')","('+s+
                 ') and not (polymer and name n,o,h)",quiet=1,mode=2,label=0,reset=1);cmd.enable("'+s+'_polar_conts")'],
              [ 1, 'just intra-side chain'  ,
                 'cmd.dist("'+s+'_polar_conts","('+s+') and not (solvent or (polymer and name n,o,h))","('+s+
                 ') and not (solvent or (polymer and name n,o,h))",quiet=1,mode=2,label=0,reset=1);cmd.enable("'+s+
                 '_polar_conts")'],
              [ 1, 'just intra-main chain'  ,
                 'cmd.dist("'+s+'_polar_conts","('+s+') and not (solvent or (polymer and not name n,o,h))","('+s+
                 ') and not (solvent or (polymer and not name n,o,h))",quiet=1,mode=2,label=0,reset=1);cmd.enable("'+s+
                 '_polar_conts")'],
              [ 0, '', '' ],
              [ 1, 'to other atoms in object',
                 'cmd.dist("'+s+'_polar_conts","('+s+')","(byobj ('+s+')) and (not ('+s+
                 '))",quiet=1,mode=2,label=0,reset=1);cmd.enable("'+s+'_polar_conts")'],
              [ 1, 'to others excluding solvent',
                 'cmd.dist("'+s+'_polar_conts","('+s+') and not solvent","(byobj ('+s+')) and (not ('+s+
                 ')) and (not solvent)",quiet=1,mode=2,label=0,reset=1);cmd.enable("'+s+'_polar_conts")'],
              [ 1, 'to any atoms',
                 'cmd.dist("'+s+'_polar_conts","('+s+')","(not '+s+
                 ')",quiet=1,mode=2,label=0,reset=1);cmd.enable("'+s+'_polar_conts")'],
              [ 1, 'to any excluding solvent',
                 'cmd.dist("'+s+'_polar_conts","('+s+') and not solvent","(not ('+s+
                 ')) and not solvent",quiet=1,mode=2,label=0,reset=1);cmd.enable("'+s+'_polar_conts")'],
              ]

def polar_inter(s):
    return [[ 2, 'Polar Contacts:', ''],
              ]

def find(s):
    return [[ 2, 'Find:', ''],
              [ 1, 'polar contacts', polar(s) ],
              ]

def align_to_object(s):
    list = cmd.get_names("public_objects",1)[0:25] # keep this practical
    list = filter(lambda x:cmd.get_type(x)=="object:molecule",list)
    result = [[ 2, 'Object:', '']]
    for a in list:
        if a!=s:
            result.append([1,a,
                           'cmd.align("polymer and name ca and ('+s+')",'+
                           '"polymer and name ca and ('+a+')",max_gap=50,quiet=0,'+
                           'object="aln_%s_to_%s",reset=1)'%(s,a)])
    return result

def align_to_sele(s):
    list = cmd.get_names("public_selections",0)[0:25] # keep this practical
    result = [[ 2, 'Selection:', '']]
    for a in list:
        if a!=s:
            result.append([1,a,
                           'cmd.align("polymer and name ca and ('+s+')",'+
                           '"polymer and name ca and ('+a+')",max_gap=50,quiet=0)'+
                           'object="aln_%s_to_%s",reset=1)'%(s,a)])
    return result

def mat_tran(s,direction=0):
    list = cmd.get_names("public_objects",1)[0:25] # keep this practical
    list = filter(lambda x:cmd.get_type(x)=="object:molecule",list)
    result = [[ 2, 'Object:', '']]
    for a in list:
        if a!=s:
            if direction:
                result.append([1,a,
                                    'cmd.matrix_transfer("'+a+'","'+s+'");'])
            else:
                result.append([1,a,
                                    'cmd.matrix_transfer("'+s+'","'+a+'");'])
    return result


def sele_align(s):
    return [[ 2, 'Align:', ''],
              [ 1, 'to molecule', align_to_object(s) ],
              [ 1, 'to selection', align_to_sele(s) ],
              [ 0, '', None ],
              [ 1, 'enabled to this', 'util.mass_align("'+s+'",1)' ],                                 
              [ 1, 'all to this', 'util.mass_align("'+s+'",0)' ],
              [ 0, '', None ],
              [ 1, 'states (*/ca)', 'cmd.intra_fit("('+s+') and name ca")' ],                        
              [ 1, 'states', 'cmd.intra_fit("'+s+'")' ],
              ]

def mol_align(s):
    return [[ 2, 'Align:', ''],
              [ 1, 'to molecule', align_to_object(s) ],
              [ 1, 'to selection', align_to_sele(s) ],
              [ 0, '', None ],
              [ 1, 'enabled to this', 'util.mass_align("'+s+'",1)' ],                                 
              [ 1, 'all to this', 'util.mass_align("'+s+'",0)' ],
              [ 0, '', None ],
              [ 1, 'states (*/ca)', 'cmd.intra_fit("('+s+') and name ca")' ],                        
              [ 1, 'states', 'cmd.intra_fit("'+s+'")' ],
              [ 0, '', None ],
              [ 1, 'matrix from', mat_tran(s,1) ],
              [ 1, 'matrix to', mat_tran(s,0) ],
              [ 1, 'matrix reset', 'cmd.matrix_reset("'+s+'")'],
              ]

def modify_sele(s):
    return [[ 2, 'Modify:', ''],
              [ 1, 'around'         , around(s)         ],           
              [ 1, 'expand'         , expand(s)         ],
              [ 1, 'extend'         , extend(s)         ],
              [ 1, 'invert'         , invert(s)         ],
              [ 1, 'complete'       , complete(s)         ],
              [ 1, 'restrict'       , restrict(s)       ],
              [ 1, 'include'        , include(s)       ],
              [ 1, 'exclude'        , exclude(s)       ]]

              
def sele_action(s):
    return [[ 2, 'Actions:'       ,''                        ],     
              [ 1, 'delete selection', 'cmd.delete("'+s+'")'          ],
              [ 1, 'rename selection', 'cmd.wizard("renaming","'+s+'")'          ],
              [ 0, ''               ,''                             ],
              [ 1, 'zoom'           ,'cmd.zoom("'+s+'",animate=-1)'            ],
              [ 1, 'orient'         ,'cmd.orient("'+s+'",animate=-1)'          ],
              [ 1, 'center'         ,'cmd.center("'+s+'",animate=-1)'            ],           
              [ 1, 'origin'         ,'cmd.origin("'+s+'")'          ],
              [ 0, ''               ,''                             ],
              [ 1, 'drag'       , 'cmd.drag("'+s+'")'    ],                        
              [ 0, ''               ,''                             ],
              [ 1, 'modify', modify_sele(s) ],
              [ 1, 'preset'         ,presets(s)         ],
              [ 1, 'find', find(s) ],
              [ 1, 'align', sele_align(s) ],
              [ 0, ''               ,''                             ],
              [ 1, 'remove atoms'   ,'cmd.remove("'+s+'");cmd.delete("'+s+'")'          ],
              [ 0, ''          ,''                                              ],
              [ 1, 'duplicate'      ,'cmd.select(None,"'+s+'")'          ], # broken...
              [ 1, 'copy to object' ,'cmd.create(None,"'+s+'",zoom=0)'     ],
              [ 1, 'extract object' ,'cmd.extract(None,"'+s+'",zoom=0)' ],
              [ 0, ''          ,''                                  ],
              [ 1, 'masking'        , masking(s)         ],
              [ 1, 'movement'       , movement(s)         ],
              [ 1, 'compute'        , compute(s) ],
              ]


def sele_action2(s):
    return [[ 2, 'Actions:'       ,''                        ],     
              [ 1, 'delete selection', 'cmd.delete("'+s+'")'          ],
              [ 1, 'rename selection', 'cmd.wizard("renaming","'+s+'")'          ],
              [ 0, ''               ,''                             ],
              [ 1, 'preset'         ,presets(s)         ],
              [ 1, 'find', find(s) ],
              [ 0, ''               ,''                             ],
              [ 1, 'remove atoms'   ,'cmd.remove("'+s+'");cmd.delete("'+s+'")'          ],
              [ 0, ''               ,''                             ],
              [ 1, 'around'         , around(s)         ],           
              [ 1, 'expand'         , expand(s)         ],
              [ 1, 'extend'         , extend(s)         ],
              [ 1, 'invert'         , invert(s)         ],
              [ 1, 'complete'       , complete(s)         ],
              [ 0, ''          ,''                                              ],
              [ 1, 'duplicate selection'      ,'cmd.select(None,"'+s+'")'          ],
              [ 1, 'copy to object'  ,'cmd.create(None,"'+s+'",zoom=0)'     ],           
              [ 1, 'extract object' ,'cmd.extract(None,"'+s+'",zoom=0)' ],
            [ 0, ''          ,''                                  ],
              [ 1, 'masking'      , masking(s)         ],
              [ 1, 'movement'       , movement(s)         ],
              [ 1, 'compute'        , compute(s)         ],           
              ]


    
def mol_action(s):
    return [[ 2, 'Actions:'     , ''                       ],     
              [ 1, 'zoom'         , 'cmd.zoom("'+s+'",animate=-1)'      ],
              [ 1, 'orient'       , 'cmd.orient("'+s+'",animate=-1)'    ],
              [ 1, 'center'         ,'cmd.center("'+s+'",animate=-1)'            ],
              [ 1, 'origin'       , 'cmd.origin("'+s+'")'    ],
              [ 0, ''               ,''                             ],
              [ 1, 'drag'       , 'cmd.drag("'+s+'")'    ],            
              [ 0, ''          ,''                                              ],
              [ 1, 'preset'  ,   presets(s)       ],
              [ 1, 'find',     find(s) ],
              [ 1, 'align',     mol_align(s) ],                      
              [ 1, 'generate'  ,   mol_generate(s)       ],           
              [ 0, ''               ,''                             ],
              [ 1, 'assign sec. struc.'  ,'cmd.dss("'+s+'")'        ],
              [ 0, ''             , ''                       ],
              [ 1, 'rename object', 'cmd.wizard("renaming","'+s+'")'          ],
              [ 1, 'duplicate object'    ,'cmd.create(None,"'+s+'")'     ],           
              [ 1, 'delete object'       , 'cmd.delete("'+s+'")'    ],
              [ 0, ''          ,''                                              ],
              [ 1, 'hydrogens' , hydrogens(s)    ],           
              [ 1, 'remove waters'  ,'cmd.remove("(solvent and ('+s+'))")'     ],
              [ 0, ''          ,''                                              ],
              [ 1, 'state'          , state(s)         ],                      
              [ 1, 'masking'        , masking(s)         ],
              [ 1, 'sequence'       , sequence(s)         ],                      
              [ 1, 'movement'       , movement(s)         ],           
              [ 1, 'compute'        , compute(s)         ],
              ]

def slice_action(s):
    return [[ 2, 'Actions:'     , ''                       ],
              [ 1, 'zoom'         , 'cmd.zoom("'+s+'",animate=-1)'      ],
              [ 1, 'center'       , 'cmd.center("'+s+'",animate=-1)'    ],           
              [ 1, 'origin'       , 'cmd.origin("'+s+'")'    ],         
              [ 0, ''             , ''                       ],
              [ 1, 'tracking on' , 'cmd.set("slice_track_camera",1,"'+s+'")'      ],
              [ 1, 'tracking off' , 'cmd.set("slice_track_camera",0,"'+s+'")'      ],           
              [ 0, ''             , ''                       ],
              [ 1, 'height map on' , 'cmd.set("slice_height_map",1,"'+s+'")'    ],
              [ 1, 'height map off', 'cmd.set("slice_height_map",0,"'+s+'")'    ],                    
              [ 0, ''             , ''                       ],
              [ 1, 'dynamic grid on' , 'cmd.set("slice_dynamic_grid",1,"'+s+'")'    ],
              [ 1, 'dynamic grid off', 'cmd.set("slice_dynamic_grid",0,"'+s+'")'    ],                    
              [ 0, ''             , ''                       ],
              [ 1, 'rename'       , 'cmd.wizard("renaming","'+s+'")'          ],           
              [ 0, ''             , ''                       ],
              [ 1, 'delete'       , 'cmd.delete("'+s+'")'    ],
              ]

def simple_action(s):
    return [[ 2, 'Actions:'     , ''                       ],
              [ 1, 'zoom'         , 'cmd.zoom("'+s+'",animate=-1)'      ],
              [ 1, 'center'       , 'cmd.center("'+s+'",animate=-1)'    ],           
              [ 1, 'origin'       , 'cmd.origin("'+s+'")'    ],
              [ 0, ''             , ''                       ],
              [ 1, 'rename'       , 'cmd.wizard("renaming","'+s+'")'          ],           
              [ 0, ''             , ''                       ],
              [ 1, 'delete'       , 'cmd.delete("'+s+'")'    ],
              ]

def map_mesh(s):
    return [[ 2, 'Mesh:',  '' ],
            [ 1, '@ level 1.0'         , 'cmd.isomesh("'+s+'_mesh","'+s+'",1.0)'      ],
            [ 0, ''             , ''                       ],            
            [ 1, '@ level 2.0'         , 'cmd.isomesh("'+s+'_mesh","'+s+'",2.0)'      ],
            [ 1, '@ level 3.0'         , 'cmd.isomesh("'+s+'_mesh","'+s+'",3.0)'      ],            
            [ 0, ''             , ''                       ],            
            [ 1, '@ level 0.0'         , 'cmd.isomesh("'+s+'_mesh","'+s+'",0.0)'      ],
            [ 1, '@ level -1.0'         , 'cmd.isomesh("'+s+'_mesh","'+s+'",1.0)'      ],
            [ 1, '@ level -2.0'         , 'cmd.isomesh("'+s+'_mesh","'+s+'",2.0)'      ],
            [ 1, '@ level -3.0'         , 'cmd.isomesh("'+s+'_mesh","'+s+'",-3.0)'      ],
            ]

def map_surface(s):
    return [[ 2, 'Surface:',  '' ],
            [ 1, '@ level 1.0'         , 'cmd.isosurface("'+s+'_surf","'+s+'",1.0)'      ],
            [ 0, ''             , ''                       ],            
            [ 1, '@ level 2.0'         , 'cmd.isosurface("'+s+'_surf","'+s+'",2.0)'      ],
            [ 1, '@ level 3.0'         , 'cmd.isosurface("'+s+'_surf","'+s+'",3.0)'      ],            
            [ 0, ''             , ''                       ],
            [ 1, '@ level 0.0'         , 'cmd.isosurface("'+s+'_surf","'+s+'",0.0)'      ],            
            [ 1, '@ level -1.0'         , 'cmd.isosurface("'+s+'_surf","'+s+'",-1.0)'      ],
            [ 1, '@ level -2.0'         , 'cmd.isosurface("'+s+'_surf","'+s+'",-2.0)'      ],
            [ 1, '@ level -3.0'         , 'cmd.isosurface("'+s+'_surf","'+s+'",-3.0)'      ],
            ]

def map_gradient(s):
    return [[ 2, 'Gradient:',  '' ],
            [ 1, 'default'         , 'cmd.gradient("'+s+'_grad","'+s+'");cmd.ramp_new("'+s+
              '_grad_ramp","'+s+'");cmd.color("'+s+'_grad_ramp","'+s+'_grad");' ]
            ]
    
def map_slice(s):
    return [[ 2, 'Slice:',  '' ],
            [ 1, 'default'         , 'cmd.slice_new("'+s+'_slice","'+s+'");cmd.ramp_new("'+s+
              '_slice_ramp","'+s+'");cmd.color("'+s+'_slice_ramp","'+s+'_slice");'+
              'cmd.set("slice_track_camera",1,"'+s+'_slice");'+
              'cmd.set("slice_dynamic_grid",1,"'+s+'_slice")'],
            ]

def map_action(s):
    return [[ 2, 'Actions:'     , ''                       ],
              [ 1, 'mesh'         , map_mesh(s)  ],
              [ 1, 'surface'      , map_surface(s)  ],
              [ 1, 'slice'        , map_slice(s)  ],
              [ 1, 'gradient'     , map_gradient(s)  ],                                    
              [ 0, ''             , ''                       ],
              [ 1, 'zoom'         , 'cmd.zoom("'+s+'",animate=-1)'      ],
              [ 1, 'center'       , 'cmd.center("'+s+'",animate=-1)'    ],           
              [ 1, 'origin'       , 'cmd.origin("'+s+'")'    ],
              [ 0, ''             , ''                       ],
              [ 1, 'rename'       , 'cmd.wizard("renaming","'+s+'")'          ],           
              [ 0, ''             , ''                       ],
              [ 1, 'delete'       , 'cmd.delete("'+s+'")'    ],
              ]

def level(s):
    return [[ 2, 'Level',  '' ],
            [ 1, 'level 5.0'         , 'cmd.isolevel("'+s+'",5.0)'      ],            
            [ 1, 'level 4.0'         , 'cmd.isolevel("'+s+'",4.0)'      ],                        
            [ 1, 'level 3.0'         , 'cmd.isolevel("'+s+'",3.0)'      ],
            [ 1, 'level 2.0'         , 'cmd.isolevel("'+s+'",2.0)'      ],
            [ 1, 'level 1.5'         , 'cmd.isolevel("'+s+'",1.5)'      ],            
            [ 1, 'level 1.0'         , 'cmd.isolevel("'+s+'",1.0)'      ],
            [ 1, 'level 0.5'         , 'cmd.isolevel("'+s+'",0.5)'      ],
            [ 1, 'level 0.0'         , 'cmd.isolevel("'+s+'",0.0)'      ],
            [ 1, 'level -0.5'         , 'cmd.isolevel("'+s+'",-0.5)'      ],
            [ 1, 'level -1.0'         , 'cmd.isolevel("'+s+'",-1.0)'      ],
            [ 1, 'level -1.0'         , 'cmd.isolevel("'+s+'",-1.5)'      ],            
            [ 1, 'level -2.0'         , 'cmd.isolevel("'+s+'",-2.0)'      ],
            [ 1, 'level -3.0'         , 'cmd.isolevel("'+s+'",-3.0)'      ],
            [ 1, 'level -4.0'         , 'cmd.isolevel("'+s+'",-4.0)'      ],
            [ 1, 'level -5.0'         , 'cmd.isolevel("'+s+'",-5.0)'      ],            
            ]

def surface_action(s):
    return [[ 2, 'Actions:'     , ''                       ],
              [ 1, 'level'         , level(s)  ],
              [ 0, ''             , ''                       ],
              [ 1, 'zoom'         , 'cmd.zoom("'+s+'",animate=-1)'      ],
              [ 1, 'center'       , 'cmd.center("'+s+'",animate=-1)'    ],           
              [ 1, 'origin'       , 'cmd.origin("'+s+'")'    ],
              [ 0, ''             , ''                       ],
              [ 1, 'rename'       , 'cmd.wizard("renaming","'+s+'")'          ],           
              [ 0, ''             , ''                       ],
              [ 1, 'delete'       , 'cmd.delete("'+s+'")'    ],
              ]

def mesh_action(s):
    return [[ 2, 'Actions:'     , ''                       ],
              [ 1, 'level'         , level(s)  ],
              [ 0, ''             , ''                       ],
              [ 1, 'zoom'         , 'cmd.zoom("'+s+'",animate=-1)'      ],
              [ 1, 'center'       , 'cmd.center("'+s+'",animate=-1)'    ],           
              [ 1, 'origin'       , 'cmd.origin("'+s+'")'    ],
              [ 0, ''             , ''                       ],
              [ 1, 'rename'       , 'cmd.wizard("renaming","'+s+'")'          ],           
              [ 0, ''             , ''                       ],
              [ 1, 'delete'       , 'cmd.delete("'+s+'")'    ],
              ]

def ramp_action(s):
    return [[ 2, 'Actions:'     , ''                       ],
              [ 1, 'delete'       , 'cmd.delete("'+s+'")'    ],
              ]

def test1(s):
        return [[ 2, 'Test1:'     , ''                      ],     
              [ 1, 'zoom'         , 'cmd.zoom("all",animate=-1)'     ],
              [ 1, 'center'   , 'cmd.center("all",animate=-1)'   ],           
              [ 1, 'origin'   , 'cmd.origin("all")'   ],
              ]

def test2(s):
        return [[ 2, 'Test2:'     , ''                      ],     
              [ 1, 'zoom'         , 'cmd.zoom("all",animate=-1)'     ],
              [ 1, 'center'   , 'cmd.center("all",animate=-1)'   ],           
              [ 1, 'origin'   , 'cmd.origin("all")'   ],
              ]

def all_action(s):
    return [[ 2, 'Actions:'     , ''                      ],     
              [ 1, 'zoom'         , 'cmd.zoom("all",animate=-1)'     ],
              [ 1, 'center'   , 'cmd.center("all",animate=-1)'   ],           
              [ 1, 'origin'   , 'cmd.origin("all")'   ],
              [ 0, ''             , ''                      ],           
              [ 1, 'preset'  , presets("all")     ],
              [ 1, 'find', find("all") ],           
              [ 0, ''          ,''                                              ],
              [ 1, 'hydrogens' ,hydrogens(s)     ],
#              [ 1, 'add hydrogens' ,'cmd.h_add("'+s+'")'     ],           
#              [ 1, 'remove hydrogens'  ,'cmd.remove("(hydro and ('+s+'))")'     ],
              [ 1, 'remove waters'  ,'cmd.remove("(solvent and ('+s+'))")'     ],                      
              [ 0, ''             , ''                      ],
              [ 1, 'delete selections'  , 'map(cmd.delete,cmd.get_names("selections"))'     ],           
              [ 0, ''          ,''                                              ],
              [ 1, 'delete everything'  , 'cmd.delete("all")'     ],           
              [ 0, ''          ,''                                              ],
              [ 1, 'masking'      , masking(s)         ],                      
              [ 1, 'movement'       , movement(s)         ],
              [ 1, 'compute'        , compute(s)         ],                      
              ]

def label_props(s):
    return [[ 2, 'Other Properties:'       ,''                        ],     
                             
              [ 1, 'formal charge' , 
  'cmd.label("'+s+'","\'%d\'%formal_charge")'                      ],
              [ 0, ''               , ''                                  ],
              [ 1, 'partial charge (0.00)' ,            
  'cmd.label("'+s+'","\'%.2f\'%partial_charge")'                      ],
              [ 1, 'partial charge (0.0000)' , 
  'cmd.label("'+s+'","\'%.4f\'%partial_charge")'                      ],
              [ 0, ''               , ''                                  ],
              [ 1, 'elec. radius'       , 'cmd.label("'+s+'","\'%1.2f\'%elec_radius")'  ],                                 
              [ 0, ''               , ''                                  ],
              [ 1, 'text type'      , 'cmd.label("'+s+'","text_type")'    ],
              [ 1, 'numeric type'   , 'cmd.label("'+s+'","numeric_type")' ],
              ]

def label_ids(s):
    return [[ 2, 'Atom Identifiers:'       ,''                        ],     
              [ 1, 'rank'           , 'cmd.label("'+s+'","rank")' ],
              [ 1, 'ID'             , 'cmd.label("'+s+'","ID")' ],
              [ 1, 'index'          , 'cmd.label("'+s+'","index")' ],           
              ]
              
def mol_labels(s):
    return [[ 2, 'Label:'        , ''                                  ],
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
              [ 1, 'other properties' , label_props(s) ],
              [ 0, ''               , ''                                  ],
              [ 1, 'atom identifiers' , label_ids(s) ],
              ]


def mol_view(s):
    return [
        [ 1, 'zoom'           ,'cmd.zoom("'+s+'",animate=-1)'            ],
        [ 1, 'orient'           ,'cmd.orient("'+s+'",animate=-1)'            ],
        [ 1, 'center'           ,'cmd.center("'+s+'",animate=-1)'            ],
        [ 1, 'origin'           ,'cmd.origin("'+s+'")'            ],
        ]

def all_option(s):
    return [
        [ 2, '(all)'      , '' ],
        [ 1, 'show'      , mol_show(s) ],
        [ 1, 'hide'      , mol_hide(s) ],
        [ 1, 'color'      , mol_color(s) ],
#      [ 1, 'view'      , mol_view(s) ],
        [ 1, 'preset'      , presets(s) ],
        [ 0, ''             , ''                      ],
        [ 1, 'zoom'           ,'cmd.zoom("'+s+'",animate=-1)'            ],
        [ 1, 'orient'           ,'cmd.orient("'+s+'",animate=-1)'            ],
        [ 1, 'center'           ,'cmd.center("'+s+'",animate=-1)'            ],
        [ 1, 'origin'           ,'cmd.origin("'+s+'")'            ],
        [ 1, 'select'        ,'cmd.select("'+s+'",enable=1,merge=2)'            ],        
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
    result = result + map(lambda ob,cm=cmmd:[1,ob,cm+ob+'")'],['all']+cmd.get_names('objects'))
    if not enable:
        result.insert(2,[1, 'selections', "util.hide_sele()"])
    else:
        result2 = [[ 2, 'Selections', '']]
        
    return result

def main_menu():
    return [
        [ 2, 'Main Pop-Up'  , '' ],
        [ 1, 'zoom (vis)'           ,'cmd.zoom("visible",animate=-1)'            ],
        [ 1, 'orient (vis)'           ,'cmd.orient("visible",animate=-1)'            ],
        [ 1, 'center (vis)'           ,'cmd.center("visible",animate=-1)'            ],      
        [ 1, 'reset'           ,'cmd.reset()'            ],
        [ 0, ''             , ''                      ],
        [ 1, 'enable', enable_disable(1) ],
        [ 1, 'disable', enable_disable(0) ],   
        [ 0, ''             , ''                      ],           
        [ 1, '(all)'      , all_option("all") ],
        [ 1, '(visible)'      , all_option("visible") ],
        [ 0, ''             , ''                      ],
        [ 1, 'ray'           ,'cmd.ray()' ],
        [ 0, ''             , ''                      ],
        [ 1, 'delete all'           ,'cmd.delete("all")' ],
        [ 1, 'reinitialize'           ,'cmd.reinitialize()' ],
        [ 1, 'quit'           ,'cmd.quit()' ],
        ]

def pick_sele_sub(s):
    result = [
        [ 2, 'Actions'  , '' ],      
        [ 1, 'rename', 'cmd.wizard("renaming","'+s+'")'          ],
        [ 1, 'clear'    , 'cmd.select("'+s+'","none")' ],
        [ 1, 'delete selection', 'cmd.delete("'+s+'")' ],
        [ 1, 'copy to object','cmd.create(None,"'+s+'",zoom=0)'            ],
        [ 1, 'extract object' ,'cmd.extract(None,"'+s+'",zoom=0)' ],
        [ 1, 'remove atoms'  , 'cmd.remove("'+s+'")' ],     
        ]
    return result

def pick_sele(title,s):
    result = [
        [ 2, title, '' ],
        [ 1, 'disable'    , 'cmd.disable("'+s+'")' ],
        [ 0, ''             , ''                      ],
        [ 1, 'actions', sele_action2(s) ],  
        [ 0, ''             , ''                      ],      
        [ 1, 'color'      , mol_color(s) ],
        [ 1, 'show'      , mol_show(s) ],
        [ 1, 'hide'      , mol_hide(s) ],
        [ 1, 'preset'  , presets(s)       ],      
        [ 1, 'label'       , mol_labels(s) ],
        [ 0, ''             , ''                      ],
        [ 1, 'zoom'           ,'cmd.zoom("'+s+'",animate=-1)'            ],
        [ 1, 'orient'           ,'cmd.orient("'+s+'",animate=-1)'            ],
        [ 1, 'center'           ,'cmd.center("'+s+'",animate=-1)'            ],
        [ 1, 'origin'           ,'cmd.origin("'+s+'")'            ],
        [ 0, ''               ,''                             ],        
        [ 1, 'drag'             ,'cmd.drag("'+s+'")'            ],
        [ 0, ''               ,''                             ],        
        [ 1, 'remove'             ,'cmd.remove("'+s+'")'            ],
        ]
    return result
    
def pick_option(title,s,object=0):
    result = [
        [ 2, title, '' ],
        [ 1, 'color'      , mol_color(s) ],
        [ 1, 'show'      , mol_show(s) ],
        [ 1, 'hide'      , mol_hide(s) ],
        [ 1, 'preset'  , presets(s)       ],      
        [ 1, 'label'          , mol_labels(s) ],
        [ 0, ''             , ''                      ],
        [ 1, 'zoom'           ,'cmd.zoom("'+s+'",animate=-1)'            ],
        [ 1, 'orient'           ,'cmd.orient("'+s+'",animate=-1)'            ],
        [ 1, 'center'           ,'cmd.center("'+s+'",animate=-1)'            ],
        [ 1, 'origin'           ,'cmd.origin("'+s+'")'            ],
        [ 1, 'select'        ,'cmd.select("'+s+'",enable=1,merge=2)'            ],
        [ 0, ''               ,''                             ],        
        [ 1, 'drag'             ,'cmd.drag("'+s+'")'            ],
        [ 1, 'masking'        , masking(s)         ],
        [ 1, 'movement'       , movement(s)         ],
        ]
  
    if object:
        result.extend([
            [ 1, 'delete'        ,'cmd.delete("'+s+'")'            ],
            [ 0, ''             , ''                      ],         
            [ 1, 'disable'        ,'cmd.disable("'+s+'")'            ],
            ])
    else:
        result.extend([
            [ 1, 'remove atoms' , 'cmd.remove("'+s+'")' ],     
            [ 0, ''             , ''                      ],      
            [ 1, 'copy to object','cmd.create(None,"'+s+'",zoom=0)'            ],
            [ 1, 'extract object' ,'cmd.extract(None,"'+s+'",zoom=0)' ],
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
        [ 1, 'label'      , mol_labels(s) ],
        [ 0, ''             , ''                      ],
        [ 1, 'zoom'           ,'cmd.zoom("'+s+'",animate=-1)'            ],
        [ 1, 'orient'           ,'cmd.orient("'+s+'",animate=-1)'            ],
        [ 1, 'center'           ,'cmd.center("'+s+'",animate=-1)'            ],
        [ 1, 'origin'           ,'cmd.origin("'+s+'")'            ],
        [ 1, 'select'        ,'cmd.select("'+s+'",enable=1,merge=2)'            ],
        [ 0, ''               ,''                             ],        
        [ 1, 'drag'             ,'cmd.drag("'+s+'")'            ],
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
        [ 1, 'create object','cmd.create(None,"'+s+'",zoom=0)'            ],
        [ 1, 'extract object' ,'cmd.extract(None,"'+s+'",zoom=0)' ],
        [ 0, ''             , ''                      ],
        [ 1, 'remove atoms' , 'cmd.remove("'+s+'")' ],     
                          ])
    return result
    
