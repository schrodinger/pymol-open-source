#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* Copyright (c) Schrodinger, LLC.
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

# This module defines the menus and their built-in commands

# shared main chain/side chain atoms
sele_ms_shared = '(n. CA|n. N&r. PRO)'

# text colors for delete and remove items
del_col = r'\933'
rem_col = r'\933'

class menucontext(object):
    '''
    Singleton context manager class to store expensive metadata
    which is used several times in submenus.

    @ivar ramps: list with ramp object names
    @ivar props: list with atom properties keys
    '''
    count = 0
    def __call__(self, cmd, sele=''):
        if self.count == 0:
            names = cmd.get_names()
            self.ramps = [n for n in names if cmd.get_type(n) == 'object:ramp']
            self.props = []
        self.count += 1
        return self
    def __enter__(self):
        return self
    def __exit__(self, type, value, traceback):
        self.count -= 1

menucontext = menucontext()

def extract(self_cmd, sele):
    return [[ 2, 'Extract', '' ],
            [ 1, 'object', 'cmd.create(None,"'+sele+'",extract="'+sele+'",zoom=0)' ],
            [ 1, 'extend 1', 'cmd.create(None,"('+sele+') extend 1",extract="'+sele+'",zoom=0)' ],
            [ 1, 'byres extend 1', 'cmd.create(None,"byres (('+sele+') extend 1)",extract="'+sele+'",zoom=0)' ],
            ]

def camera_store_with_scene(self_cmd,frame):
    list = self_cmd.get_scene_list()[0:40] # keep this practical
    result = [[ 2, 'Scene:', '']]
    for a in list:
        result.append([1,a,'cmd.mview("store",scene="'+a+'",first=%s)'%frame])
    return result


def store_with_state(self_cmd,obj='',frame=0):
    n_state = self_cmd.count_states('%' + obj if obj else '')
    result = [[ 2, 'State:', ''],
              [ 1, 'current','cmd.mview("store",object="'+obj+'",state=-1,first=%s)'%(frame)],
              [ 0, ''               ,''                   ],
              [ 1, '1', 'cmd.mview("store",object="'+obj+'",state=1,first=%s)'%(frame)   ],
              ]
    if (n_state>1):
        result.extend([
              [ 1, str(n_state),'cmd.mview("store",object="'+obj+'",state=%d,first=%s)'%(n_state,frame)],
              [ 0, ''               ,''                   ],
              ])
        n_show = min(8, n_state)
        result.extend([
              [ 1, str(state), 'cmd.mview("store",object="'+obj+'",state=%d,first=%s)'%(state,frame) ]
              for state in [(n_state * i) // n_show for i in range(2, n_show)]
              ])

    return result

def mouse_config(self_cmd):
    result = [[ 1, '3-Button Motions',
                'cmd.config_mouse("three_button_motions")' ],
              [ 1, '3-Button Editing',
                'cmd.config_mouse("three_button_editing")' ],
              [ 1, '3-Button Viewing',
                'cmd.mouse("three_button_viewing")' ],
              [ 1, '3-Button Lights',
                'cmd.mouse("three_button_lights")' ],
              [ 1, '3-Button All Modes',
                'cmd.config_mouse("three_button_all_modes")' ],
              [ 0, '', ''],
              [ 1, '2-Button Editing',
                'cmd.config_mouse("two_button_editing")' ],
              [ 1, '2-Button Viewing',
                'cmd.config_mouse("two_button_viewing")' ],
              [ 1, '2-Button Lights',
                'cmd.mouse("two_button_lights")' ],
              ]
    return result

def smooth(self_cmd,extra=''):
    return [[ 1, 'a little'  ,   'cmd.mview("smooth"%s)'%extra            ],
            [ 1, 'more'     ,   'cmd.mview("smooth",window=15%s)'%extra ],
            [ 1, 'a lot'   ,   'cmd.mview("smooth",window=30%s)'%extra ]]

def camera_motion(self_cmd, frame="0"):
    return [[ 2, 'Camera Motion:'     , ''                       ],
            [ 1, 'store'         , 'cmd.mview("store",first='+frame+')'      ],
            [ 1, 'store with scene' , camera_store_with_scene(self_cmd,frame) ],
            [ 1, 'store with state' , store_with_state(self_cmd,'',frame) ],
            [ 1, 'clear'       ,   'cmd.mview("clear",first='+frame+')'      ],
            [ 0, ''               ,''                             ],
            [ 1, 'reset camera motions'   , 'cmd.mview("reset")'   ],
            [ 0, ''               ,''                             ],
            [ 1, 'purge entire movie'   , 'cmd.mset()'   ],
            [ 0, ''               ,''                             ],
            [ 1, 'smooth key frames'  ,   smooth(self_cmd)     ],
            [ 0, ''               ,''                             ],
            [ 1, 'interpolate'   , 'cmd.mview("interpolate")'   ],
            [ 1, 'reinterpolate'   , 'cmd.mview("reinterpolate")'   ],
            [ 1, 'uninterpolate'   , 'cmd.mview("uninterpolate")'   ],
            ]

def obj_motion(self_cmd, obj, frame="0"):
    return [[ 2, 'Object "'+obj+'" Motion:'     , ''                       ],
            [ 1, 'drag'       ,   'cmd.drag("'+obj+'")'    ],
            [ 0, ''               ,''                             ],
            [ 1, 'store'         , 'cmd.mview("store",object="'+obj+'",first='+frame+')'      ],
            [ 1, 'store with state' , store_with_state(self_cmd,obj,frame) ],
            [ 1, 'reset'       ,   ';cmd.reset(object="'+obj+'");'    ],
            [ 1, 'clear'       ,   'cmd.mview("clear",object="'+obj+'",first='+frame+')'    ],
            [ 0, ''               ,''                             ],
            [ 1, 'reset object motions'       ,   'cmd.mview("reset",object="'+obj+'")'    ],
            [ 1, 'purge object motions'       ,   'cmd.mview("purge",object="'+obj+'")'    ],
            [ 0, ''               ,''                             ],
            [ 1, 'smooth key frames' , smooth(self_cmd,',object="'+obj+'"') ],
            [ 0, ''               ,''                             ],
            [ 1, 'interpolate'   ,   'cmd.mview("interpolate",object="'+obj+'")'    ],
            [ 1, 'reinterpolate'   ,   'cmd.mview("reinterpolate",object="'+obj+'")'    ],
            [ 1, 'uninterpolate'   ,   'cmd.mview("uninterpolate",object="'+obj+'")'    ],
            ]

def rep_action(self_cmd, sele, action) :
    flag_ignore_action = "set" if action == "hide" else "clear"
    return [
        [ 1, 'wire'       , 'cmd.'+action+'("wire"      ,"'+sele+'")' ],
        [ 1, '  lines'    , 'cmd.'+action+'("lines"     ,"'+sele+'")' ],
        [ 1, '  nonbonded', 'cmd.'+action+'("nonbonded" ,"'+sele+'")' ],
        [ 0, ''           , ''                               ],
        [ 1, 'licorice'   , 'cmd.'+action+'("licorice"  ,"'+sele+'")' ],
        [ 1, '  sticks'   , 'cmd.'+action+'("sticks"    ,"'+sele+'")' ],
        [ 1, '  nb_spheres','cmd.'+action+'("nb_spheres","'+sele+'")' ],
        [ 0, ''           , ''                               ],
        [ 1, 'ribbon'     , 'cmd.'+action+'("ribbon"    ,"'+sele+'")' ],
        [ 1, 'cartoon'    , 'cmd.'+action+'("cartoon"   ,"'+sele+'")' ],
        [ 0, ''           , ''                               ],
        [ 1, 'label'      , 'cmd.'+action+'("labels"    ,"'+sele+'")' ],
        [ 1, 'cell'       , 'cmd.'+action+'("cell"      ,"'+sele+'")' ],
        [ 0, ''           , ''                               ],
        [ 1, 'dots'       , 'cmd.'+action+'("dots"      ,"'+sele+'")' ],
        [ 1, 'spheres'    , 'cmd.'+action+'("spheres"   ,"'+sele+'")' ],
        [ 0, ''           , ''                               ],
        [ 1, 'mesh'       , 'cmd.'+action+'("mesh"      ,"'+sele+'")' ],
        [ 1, 'surface'    , 'cmd.'+action+'("surface"   ,"'+sele+'")' ],
        [ 1, 'flag ignore', [
            [ 2, 'flag ignore', ''],
            [ 1, flag_ignore_action, f'cmd.flag("ignore",{sele!r},{flag_ignore_action!r});cmd.rebuild({sele!r})' ],
            [ 0, '', '' ],
            [ 2, '\\272Note:\\559 Atoms with the "ignore"', ''],
            [ 2, '\\559flag are ignored during surface', ''],
            [ 2, '\\559calculation. By default, all', ''],
            [ 2, '\\559non-polymer atoms are ignored.', ''],
        ]],
        ]

def mol_as(self_cmd, sele):
    return (
        [[ 2, 'As:'   , '']]
        +rep_action(self_cmd, sele,'show_as')
        )

def mol_toggle(self_cmd, sele):
    return (
        [[ 2, 'As:'   , '']]
        +rep_action(self_cmd, sele,'toggle')
        )

def show_misc(self_cmd, sele):
    return [[ 2, 'Show:', '' ],
            [ 1, 'lines', 'cmd.show("lines","'+sele+'")'],
            [ 1, 'sticks', 'cmd.show("sticks","'+sele+'")'],
            [ 1, 'spheres', 'cmd.show("spheres","'+sele+'")'],
            ]

def mol_show(self_cmd, sele):
    return (
        [[ 2, 'Show:'      , ''                               ],
         [ 1, 'as'    , mol_as(self_cmd, sele) ],
         [ 0, '', '']]
        + rep_action(self_cmd, sele,'show') +
        [[ 0, '', ''],
         [ 1, 'organic' , show_misc(self_cmd, '(organic and ('+sele+'))') ],
         [ 1, 'main chain' , show_misc(self_cmd, "((byres ("+sele+"))&(bb.))") ],
         [ 1, 'side chain' , show_misc(self_cmd, "((byres ("+sele+"))&(sc.|"+sele_ms_shared+"))") ],
         [ 1, 'disulfides' , show_misc(self_cmd, "(byres ((("+sele+
            ") & r. CYS+CYX & n. SG) & bound_to (("+sele+") & r. CYS+CYX & n. SG))) & n. CA+CB+SG") ]
         ] +
         [[ 0, '', ''],
           [ 1, 'valence', 'cmd.set_bond("valence", "1", "'+sele+'",quiet=1)'],
         ] )

    self_cmd.show("lines","(byres ((" + sele + " & r. CYS+CYX & n. SG) & bound_to ("
                  + sele + " & r. CYS+CYX & n. SG))) & n. CA+CB+SG")

def hide_hydro(self_cmd, sele):
    return ( [[ 2, 'Hide:'     , ''                                ],
              [ 1, 'all' , 'cmd.hide("('+sele+' and hydro)")'   ],
              [ 1, 'nonpolar' , 'cmd.hide("('+sele+' and hydro and (elem C extend 1))")' ],
              ] )

def mol_hide(self_cmd, sele):
    return (
        [[ 2, 'Hide:'     , ''                                ],
         [ 1, 'everything', 'cmd.hide("everything","'+sele+'")'  ],
         [ 0, ''          , ''                                ]]
        + rep_action(self_cmd, sele,'hide') +
        [[ 0, ''          , ''                                ],
         [ 1, 'main chain', 'cmd.hide("((byres ('+sele+'))&(bb.&!'+sele_ms_shared+'))")' ],
         [ 1, 'side chain', 'cmd.hide("((byres ('+sele+'))&(sc.&!'+sele_ms_shared+'))")' ],
         [ 1, 'waters'    , 'cmd.hide("(solvent and ('+sele+'))")'     ],
         [ 0, ''          , ''                                ],
         [ 1, 'hydrogens' , hide_hydro(self_cmd, sele) ],
         [ 0, ''          , ''                                ],
         [ 1, 'unselected', 'cmd.hide("(not '+sele+')")'         ],
         ]
         + [[ 0, '', ''],
           [ 1, 'valence', 'cmd.set_bond("valence", "0", "'+sele+'",quiet=1)'],
         ] )


def measurement_show(self_cmd, sele):
    return [[ 2, 'Show:'     , ''                               ],
              [ 1, 'dashes'    , 'cmd.show("dashes"    ,"'+sele+'")' ],
              [ 1, 'angles'    , 'cmd.show("angles"    ,"'+sele+'")' ],
              [ 1, 'dihedrals' , 'cmd.show("dihedrals" ,"'+sele+'")' ],
              [ 1, 'labels'    , 'cmd.show("labels"    ,"'+sele+'")' ]
             ]

def measurement_hide(self_cmd, sele):
    return [[ 2, 'Hide:'     , ''                                ],
              [ 1, 'dashes'    , 'cmd.hide("dashes"    ,"'+sele+'")'  ],
              [ 1, 'angles'    , 'cmd.hide("angles"    ,"'+sele+'")'  ],
              [ 1, 'dihedrals' ,  'cmd.hide("dihedrals"  ,"'+sele+'")'  ],
              [ 1, 'labels'    , 'cmd.hide("labels"    ,"'+sele+'")'  ]
             ]

def cgo_show(self_cmd, sele):
    return [[ 2, 'Show:'     , ''                               ],
              [ 1, 'cgo'    , 'cmd.show("cgo"    ,"'+sele+'")' ],
             ]

def cgo_hide(self_cmd, sele):
    return [[ 2, 'Hide:'     , ''                                ],
              [ 1, 'cgo'    , 'cmd.hide("cgo"    ,"'+sele+'")'  ],
             ]

def simple_show(self_cmd, sele):
    return [[ 2, 'Show:'       , ''                             ],
              [ 1, 'everything'  , 'cmd.show("everything","'+sele+'")'          ]]

def simple_hide(self_cmd, sele):
    return [[ 2, 'Hide:'     ,''                                ],
              [ 1, 'everything'    ,'cmd.hide("everything","'+sele+'")'        ]]

def map_show(self_cmd, sele):
    return [[ 2, 'Show:'       , ''                             ],
              [ 1, 'dots'        , 'cmd.show("dots","'+sele+'")'     ],
              [ 1, 'extent'        , 'cmd.show("extent","'+sele+'")'     ],
              [ 1, 'everything'  , 'cmd.show("everything","'+sele+'")'          ]]

def map_hide(self_cmd, sele):
    return [[ 2, 'Hide:'     ,''                                ],
              [ 1, 'dots'        , 'cmd.hide("dots","'+sele+'")'     ],
              [ 1, 'extent'      , 'cmd.hide("extent","'+sele+'")'     ],
              [ 1, 'everything'    ,'cmd.hide("everything","'+sele+'")'        ]]

def mesh_show(self_cmd, sele):
    return [[ 2, 'Show:'       , ''                             ],
              [ 1, 'mesh'        , 'cmd.show("mesh","'+sele+'")'     ],
              [ 1, 'cell'        , 'cmd.show("cell","'+sele+'")'     ],
              [ 1, 'everything'  , 'cmd.show("everything","'+sele+'")'          ]]

def mesh_hide(self_cmd, sele):
    return [[ 2, 'Hide:'       , ''                             ],
              [ 1, 'mesh'        , 'cmd.hide("mesh","'+sele+'")'     ],
              [ 1, 'cell'        , 'cmd.hide("cell","'+sele+'")'      ],
              [ 1, 'everything'  , 'cmd.hide("everything","'+sele+'")'          ]]

def surface_show(self_cmd, sele):
    return [[ 2, 'Show:'       , ''                             ],
              [ 1, 'surface'        , 'cmd.show("surface","'+sele+'")'     ],
              [ 1, 'cell'        , 'cmd.show("cell","'+sele+'")'     ],
              [ 1, 'everything'  , 'cmd.show("everything","'+sele+'")'          ]]

def surface_hide(self_cmd, sele):
    return [[ 2, 'Hide:'       , ''                             ],
              [ 1, 'surface'        , 'cmd.hide("surface","'+sele+'")'     ],
              [ 1, 'cell'        , 'cmd.hide("cell","'+sele+'")'      ],
              [ 1, 'everything'  , 'cmd.hide("everything","'+sele+'")'          ]]

def slice_show(self_cmd, sele):
    return [[ 2, 'Show:'       , ''                             ],
              [ 1, 'slice'       , 'cmd.show("slice","'+sele+'")'     ],
              ]

def slice_hide(self_cmd, sele):
    return [[ 2, 'Hide:'       , ''                             ],
              [ 1, 'slice'        , 'cmd.hide("slice","'+sele+'")'     ],
              ]

def volume_show(self_cmd, sele):
    return [[ 2, 'Show:'       , ''                             ],
              [ 1, 'volume'       , 'cmd.show("volume","'+sele+'")'     ],
              [ 1, 'extent'       , 'cmd.show("extent","'+sele+'")'     ],
              ]

def volume_hide(self_cmd, sele):
    return [[ 2, 'Hide:'       , ''                             ],
              [ 1, 'volume'        , 'cmd.hide("volume","'+sele+'")'     ],
              [ 1, 'extent'        , 'cmd.hide("extent","'+sele+'")'     ],
              ]

def by_elem2(self_cmd, sele):
    return [
        [ 2, 'Atoms'     ,''                               ],
        [1,'\\494C\\777H\\229N\\922O\\950S...','util.cba(10,"'+sele+'",_self=cmd)'],# lime
        [1,'\\155C\\777H\\229N\\922O\\950S...','util.cba(5262,"'+sele+'",_self=cmd)'],# deepteal
        [1,'\\904C\\777H\\229N\\922O\\950S...','util.cba(12,"'+sele+'",_self=cmd)'],# hotpink
        [1,'\\983C\\777H\\229N\\922O\\950S...','util.cba(36,"'+sele+'",_self=cmd)'],# yelloworange
        [1,'\\525C\\777H\\229N\\922O\\950S...','util.cba(5271,"'+sele+'",_self=cmd)'],# violetpurple
        [1,'\\666C\\777H\\229N\\922O\\950S...','util.cba(124,"'+sele+'",_self=cmd)'],# grey70
        [1,'\\049C\\777H\\229N\\922O\\950S...','util.cba(17,"'+sele+'",_self=cmd)'],# marine
        [1,'\\760C\\777H\\229N\\922O\\950S...','util.cba(18,"'+sele+'",_self=cmd)'],# olive
        ]

def by_elem3(self_cmd, sele):
    return [
        [ 2, 'Atoms'     ,''                               ],
        [1,'\\564C\\777H\\229N\\922O\\950S...','util.cba(5270,"'+sele+'",_self=cmd)'],# smudge
        [1,'\\077C\\777H\\229N\\922O\\950S...','util.cba(20,"'+sele+'",_self=cmd)'],# teal
        [1,'\\644C\\777H\\229N\\922O\\950S...','util.cba(5272,"'+sele+'",_self=cmd)'],# dirtyviolet
        [1,'\\976C\\777H\\229N\\922O\\950S...','util.cba(52,"'+sele+'",_self=cmd)'],# wheat
        [1,'\\944C\\777H\\229N\\922O\\950S...','util.cba(5258,"'+sele+'",_self=cmd)'],# deepsalmon
        [1,'\\978C\\777H\\229N\\922O\\950S...','util.cba(5274,"'+sele+'",_self=cmd)'],# lightpink
        [1,'\\499C\\777H\\229N\\922O\\950S...','util.cba(5257,"'+sele+'",_self=cmd)'],# aquamarine
        [1,'\\994C\\777H\\229N\\922O\\950S...','util.cba(5256,"'+sele+'",_self=cmd)'],# paleyellow
        ]

def by_elem4(self_cmd, sele):
    return [
        [ 2, 'Atoms'     ,''                               ],
        [1,'\\094C\\777H\\229N\\922O\\950S...','util.cba(15,"'+sele+'",_self=cmd)'],# limegreen
        [1,'\\247C\\777H\\229N\\922O\\950S...','util.cba(5277,"'+sele+'",_self=cmd)'],# skyblue
        [1,'\\824C\\777H\\229N\\922O\\950S...','util.cba(5279,"'+sele+'",_self=cmd)'],# warmpink
        [1,'\\792C\\777H\\229N\\922O\\950S...','util.cba(5276,"'+sele+'",_self=cmd)'],# limon
        [1,'\\949C\\777H\\229N\\922O\\950S...','util.cba(53,"'+sele+'",_self=cmd)'],# violet
        [1,'\\889C\\777H\\229N\\922O\\950S...','util.cba(5278,"'+sele+'",_self=cmd)'],# bluewhite
        [1,'\\297C\\777H\\229N\\922O\\950S...','util.cba(5275,"'+sele+'",_self=cmd)'],# greencyan
        [1,'\\653C\\777H\\229N\\922O\\950S...','util.cba(5269,"'+sele+'",_self=cmd)'],# sand
        ]

def by_elem5(self_cmd, sele):
    return [
        [ 2, 'Atoms'     ,''                               ],
[1,'\\252C\\777H\\229N\\922O\\950S...','util.cba(22,"'+sele+'",_self=cmd)'],# forest
[1,'\\466C\\777H\\229N\\922O\\950S...','util.cba(5266,"'+sele+'",_self=cmd)'],# lightteal
[1,'\\755C\\777H\\229N\\922O\\950S...','util.cba(5280,"'+sele+'",_self=cmd)'],# darksalmon
[1,'\\570C\\777H\\229N\\922O\\950S...','util.cba(5267,"'+sele+'",_self=cmd)'],# splitpea
[1,'\\634C\\777H\\229N\\922O\\950S...','util.cba(5268,"'+sele+'",_self=cmd)'],# raspberry
[1,'\\555C\\777H\\229N\\922O\\950S...','util.cba(104,"'+sele+'",_self=cmd)'],# grey50
[1,'\\226C\\777H\\229N\\922O\\950S...','util.cba(23,"'+sele+'",_self=cmd)'],# deepblue
[1,'\\632C\\777H\\229N\\922O\\950S...','util.cba(51,"'+sele+'",_self=cmd)'],# brown
              ]

def by_elem6(self_cmd, sele):
    return [
        [ 2, 'Atoms'     ,''                               ],
[1,'\\191C\\911H\\229N\\922O\\950S...','util.cbh("tv_red","'+sele+'",_self=cmd)'],# tv_red
[1,'\\191C\\917H\\229N\\922O\\950S...','util.cbh("lightmagenta","'+sele+'",_self=cmd)'],# lightmagenta
[1,'\\191C\\119H\\229N\\922O\\950S...','util.cbh("tv_blue","'+sele+'",_self=cmd)'],# tv_blue
[1,'\\191C\\940H\\229N\\922O\\950S...','util.cbh("orange","'+sele+'",_self=cmd)'],# orange
[1,'\\191C\\870H\\229N\\922O\\950S...','util.cbh("olive","'+sele+'",_self=cmd)'],# olive
[1,'\\191C\\088H\\229N\\922O\\950S...','util.cbh("teal","'+sele+'",_self=cmd)'],# teal
[1,'\\191C\\521H\\229N\\922O\\950S...','util.cbh("chocolate","'+sele+'",_self=cmd)'],# chocolate
[1,'\\191C\\000H\\229N\\922O\\950S...','util.cbh("black","'+sele+'",_self=cmd)'],# black
              ]

def by_elem(self_cmd, sele):
    return [
        [ 2, 'Atoms'     ,''                               ],
        [1,' \\777H\\229N\\922O\\950S...','util.cnc("'+sele+'",_self=cmd)'],

[1,'\\292C\\777H\\229N\\922O\\950S...','util.cba(33,"'+sele+'",_self=cmd)'],# tv_green
[1,'\\099C\\777H\\229N\\922O\\950S...','util.cba(5,"'+sele+'",_self=cmd)'],# cyan
[1,'\\927C\\777H\\229N\\922O\\950S...','util.cba(154,"'+sele+'",_self=cmd)'],# lightmagenta
[1,'\\990C\\777H\\229N\\922O\\950S...','util.cba(6,"'+sele+'",_self=cmd)'],# yellow
[1,'\\955C\\777H\\229N\\922O\\950S...','util.cba(9,"'+sele+'",_self=cmd)'],# salmon
[1,'\\888C\\777H\\229N\\922O\\950S...','util.cba(144,"'+sele+'",_self=cmd)'],# grey90
[1,'\\449C\\777H\\229N\\922O\\950S...','util.cba(11,"'+sele+'",_self=cmd)'],# slate
[1,'\\962C\\777H\\229N\\922O\\950S...','util.cba(13,"'+sele+'",_self=cmd)'],# orange
        [ 1, 'set 2'     ,by_elem2(self_cmd, sele)                    ],
        [ 1, 'set 3'     ,by_elem3(self_cmd, sele)                    ],
        [ 1, 'set 4'     ,by_elem4(self_cmd, sele)                    ],
        [ 1, 'set 5'     ,by_elem5(self_cmd, sele)                    ],
        [ 1, 'set 6/H'   ,by_elem6(self_cmd, sele)                    ],
              ]

def by_ss(self_cmd, sele):
    return [
                [ 2, 'By Secondary Structure:'     ,''                               ],
    [ 1, '\\900Helix \\990Sheet \\090Loop'  , 'util.cbss("'+sele+'","red","yellow","green",_self=cmd)'],
    [ 1, '\\099Helix \\909Sheet \\955Loop'  , 'util.cbss("'+sele+'","cyan","magenta","salmon",_self=cmd)'],
    [ 1, '\\099Helix \\900Sheet \\909Loop'  , 'util.cbss("'+sele+'","cyan","red","magenta",_self=cmd)'],
              ]

def by_rep(self_cmd, sele, mode=0):
    r = [[ 2, 'By Representation:', '' ]]
    r += [
        [1, rep, by_rep_sub(self_cmd, rep, setting, sele)] if rep else [0, '', '']
        for (rep, setting) in rep_setting_lists[mode]
    ]
    return r

def by_rep_sub(self_cmd, rep, setting, sele):
    expr = 'cmd.set("%s", "{0}", %s, quiet=0)' % (setting, repr(sele))
    r = [[ 2, rep.capitalize(), '' ]]
    r += all_colors_generic(self_cmd, expr)
    r += [
        [ 0 , '', '' ],
        [ 1 , 'unset', 'cmd.unset("%s", %s, quiet=0)' % (setting, repr(sele)) ],
    ]
    return r

def spectrum(self_cmd, sele):
    r = [
        [ 2, 'Spectrum:'     ,''                               ],
        [ 1, '\\900r\\950a\\990i\\090n\\099b\\059o\\009w\\888(elem C)',
          'cmd.spectrum("count",selection="('+sele+')&elem C")'],
        [ 1, '\\900r\\950a\\990i\\090n\\099b\\059o\\009w\\888(*/CA)',
          'cmd.spectrum("count",selection="('+sele+')&*/CA")'],
        [ 1, '\\900r\\950a\\990i\\090n\\099b\\059o\\009w',
          'cmd.spectrum("count",selection="'+sele+'",byres=1)'],
        [ 0, ''                                , ''                 ],
        [ 1, 'b-factors'   ,   'cmd.spectrum("b",selection=("'+sele+'"),quiet=0)'         ],
        [ 1, 'b-factors(*/CA)' , 'cmd.spectrum("b",selection="(('+sele+')&*/CA)",quiet=0)'         ],
        [ 0, ''                                , ''                 ],
        [ 1, 'area (molecular)', 'util.color_by_area(("'+sele+'"),"molecular")'         ],
        [ 1, 'area (solvent)'  , 'util.color_by_area(("'+sele+'"),"solvent")'         ],
        ]
    return r

def by_chain(self_cmd, sele):
    by_segi = r'\900b\950y \090s\099e\059g\705i\888 '
    return [
        [ 2, 'By Chain:'     ,''                               ],
              [ 1, '\\900b\\950y \\090c\\099h\\059a\\009i\\705n\\888(elem C)',
                 'util.color_chains("('+sele+' and elem C)",_self=cmd)'],
              [ 1, '\\900b\\950y \\090c\\099h\\059a\\009i\\705n\\888(*/CA)',
                 'util.color_chains("('+sele+' and name CA)",_self=cmd)'],
              [ 1, '\\900b\\950y \\090c\\099h\\059a\\009i\\705n',
                 'util.color_chains("('+sele+')",_self=cmd)'],
                      [ 0, ''                                , ''                 ],
              [ 1, '\\900c\\950h\\990a\\090i\\099n\\059b\\009o\\705w\\888s',
                 'util.chainbow("('+sele+')",_self=cmd)'],
              [ 0, '', '' ],
              [ 1, by_segi + '(elem C)', 'cmd.spectrum("segi","rainbow","('+sele+') & elem C")'],
              [ 1, by_segi, 'cmd.spectrum("segi","rainbow","' + sele + '")' ],
        ]

rep_setting_lists = [
    [
        # molecule
        ('lines', 'line_color'),
        ('sticks', 'stick_color'),
        ('ribbon', 'ribbon_color'),
        ('cartoon', 'cartoon_color'),
        ('', ''),
        ('labels', 'label_color'),
        ('', ''),
        ('dots', 'dot_color'),
        ('spheres', 'sphere_color'),
        ('', ''),
        ('mesh', 'mesh_color'),
        ('surface', 'surface_color'),
    ],
    [
        # measurement
        ('dashes', 'dash_color'),
        ('angles', 'angle_color'),
        ('dihedrals', 'dihedral_color'),
        ('labels', 'label_color'),
    ],
    [
        # extra
        ('', 'cartoon_highlight_color'),
        ('', 'cartoon_ladder_color'),
        ('', 'cartoon_nucleic_acid_color'),
        ('', 'cartoon_ring_color'),
        ('', 'ellipsoid_color'),
        ('', 'label_outline_color'),
        ('', 'ray_interior_color'),
        ('', 'ray_trace_color'),
        ('', 'stick_ball_color'),
    ],
]

all_colors_list = [
    ('reds', [
        ('900', 'red'),
        ('922', 'tv_red'),
        ('634', 'raspberry'),
        ('755', 'darksalmon'),
        ('955', 'salmon'),
        ('944', 'deepsalmon'),
        ('824', 'warmpink'),
        ('611', 'firebrick'),
        ('522', 'ruby'),
        ('521', 'chocolate'),
        ('632', 'brown'),
    ]),
    ('greens', [
        ('090', 'green'),
        ('292', 'tv_green'),
        ('490', 'chartreuse'),
        ('570', 'splitpea'),
        ('564', 'smudge'),
        ('686', 'palegreen'),
        ('094', 'limegreen'),
        ('494', 'lime'),
        ('792', 'limon'),
        ('252', 'forest'),
    ]),
    ('blues', [
        ('009', 'blue'),
        ('339', 'tv_blue'),
        ('049', 'marine'),
        ('449', 'slate'),
        ('779', 'lightblue'),
        ('247', 'skyblue'),
        ('409', 'purpleblue'),
        ('226', 'deepblue'),
        ('115', 'density'),
    ]),
    ('yellows', [
        ('990', 'yellow'),
        ('992', 'tv_yellow'),
        ('994', 'paleyellow'),
        ('983', 'yelloworange'),
        ('792', 'limon'),
        ('976', 'wheat'),
        ('653', 'sand'),
    ]),
    ('magentas', [
        ('909', 'magenta'),
        ('927', 'lightmagenta'),
        ('904', 'hotpink'),
        ('968', 'pink'),
        ('978', 'lightpink'),
        ('644', 'dirtyviolet'),
        ('949', 'violet'),
        ('525', 'violetpurple'),
        ('707', 'purple'),
        ('515', 'deeppurple'),
    ]),
    ('cyans', [
        ('099', 'cyan'),
        ('799', 'palecyan'),
        ('499', 'aquamarine'),
        ('297', 'greencyan'),
        ('077', 'teal'),
        ('155', 'deepteal'),
        ('466', 'lightteal'),
    ]),
    ('oranges', [
        ('950', 'orange'),
        ('951', 'tv_orange'),
        ('962', 'brightorange'),
        ('985', 'lightorange'),
        ('983', 'yelloworange'),
        ('760', 'olive'),
        ('551', 'deepolive'),
    ]),
    ('tints', [
        ('976', 'wheat'),
        ('686', 'palegreen'),
        ('779', 'lightblue'),
        ('994', 'paleyellow'),
        ('978', 'lightpink'),
        ('799', 'palecyan'),
        ('985', 'lightorange'),
        ('889', 'bluewhite'),
    ]),
    ('grays', [
        ('999', 'white'),
        ('999', 'gray90'),
        ('888', 'gray80'),
        ('777', 'gray70'),
        ('666', 'gray60'),
        ('555', 'gray50'),
        ('444', 'gray40'),
        ('333', 'gray30'),
        ('222', 'gray20'),
        ('222', 'gray10'),
        ('222', 'black'),
    ]),
]

def colorramps(self_cmd, expr):
    with menucontext(self_cmd) as mc:
        return [[ 1, name, expr.format(name) ]
                for name in mc.ramps]

def all_colors_generic(self_cmd, expr):
    r = [
        [ 1, '\\' + c_list[0][0] + gn, [[2, gn.capitalize(), '']] + [
            [1, '\\' + c3 + cn, expr.format(cn)]
            for (c3, cn) in c_list
        ]]
        for (gn, c_list) in all_colors_list
    ]
    ramps = colorramps(self_cmd, expr)
    if ramps:
        r += [[ 0, '', '' ], [ 1, 'ramps', ramps ]]
    return r

def all_colors(self_cmd, sele):
    expr = 'cmd.color_deep("{0}", ' + repr(sele) + ', 0)'
    with menucontext(self_cmd, sele):
        return all_colors_generic(self_cmd, expr)

def vol_color(self_cmd, sele):
    from pymol.colorramping import namedramps
    rsele = repr(sele)
    return [
        [2, 'Coloring:', '' ],
        [1, 'panel', 'cmd.volume_panel(%s)' % (rsele)],
        [0, '', ''],
    ] + [
        [1, p, 'cmd.volume_color(%s, "%s")' % (rsele, p) ]
        for p in sorted(namedramps)
    ]

def slice_color(self_cmd, sele):
    expr = 'cmd.color("{0}", %s)' % (repr(sele))
    return colorramps(self_cmd, expr)

def color_auto(self_cmd, sele):
    return [
        [ 2, 'Auto'     ,''                               ],
        [ 1, 'elem C', 'cmd.color("auto","('+sele+') and elem C")' ],
        [ 0, ''                                , ''                 ],
        [ 1, 'all','cmd.color("auto","'+sele+'")' ],
        [ 0, ''                                , ''                 ],
        [ 1, '\\900b\\950y \\090o\\099b\\059j\\999(elem C)',
          'util.color_objs("('+sele+' and elem C)",_self=cmd)'],
        [ 1, '\\900b\\950y \\090o\\099b\\059j',
          'util.color_objs("('+sele+')",_self=cmd)'],
        ]

def mol_color(self_cmd, sele):
    with menucontext(self_cmd, sele):
      return (
        [[ 2, 'Color:'     ,''                               ],
         [ 1, 'by element'  , by_elem(self_cmd, sele) ],
         [ 1, 'by chain' , by_chain(self_cmd, sele) ],
         [ 1, 'by ss  '  , by_ss(self_cmd, sele) ],
         [ 1, 'by rep'  , by_rep(self_cmd, sele) ],
         [ 1, '\\900s\\950p\\990e\\090c\\099t\\059r\\009u\\555m', spectrum(self_cmd, sele) ],
         [ 0, ''                                , ''                 ],
         [ 1, 'auto', color_auto(self_cmd, sele) ],
         [ 0, ''                                , ''                 ],
         ] +
        all_colors(self_cmd, sele))

def measurement_color(self_cmd, sele):
    r = [
        [ 2, 'Color:', '' ],
        [ 1, 'by rep', by_rep(self_cmd, sele, 1) ],
        [ 0, '', '' ],
    ]
    r += all_colors(self_cmd, sele)
    return r

def mesh_color(self_cmd, name, rep='mesh'):
    expr = ('cmd.set("%s_negative_visible",1,"%s",quiet=0);'
            'cmd.set("%s_negative_color","{0}","%s",quiet=0)' % (rep, name, rep, name))
    with menucontext(self_cmd, name):
        negative = all_colors_generic(self_cmd, expr)
        return [
            [ 2, 'Color:', '' ],
            [ 1, 'negative'  , [
                [ 2, 'Negative Color:', '' ],
                [ 1, 'off', 'cmd.set("%s_negative_visible",0,"%s",quiet=0);' % (rep, name) ],
                [ 0, '', '' ],
            ] + negative ],
            [ 0, '', '' ],
        ] + all_colors(self_cmd, name)

def general_color(self_cmd, sele):
    return [[ 2, 'Color:'     ,''                        ]] + all_colors(self_cmd, sele)

def preset_ligand_sites(self_cmd, sele):
    return [[ 2, 'Ligand Sites:', ''],
              [ 1, 'cartoon'   , 'preset.ligand_cartoon("'+sele+'",_self=cmd)'          ],
              [ 0, '', ''],
              [ 1, 'solid surface'   , 'preset.ligand_sites("'+sele+'",_self=cmd)'          ],
              [ 1, 'solid (better)'   , 'preset.ligand_sites_hq("'+sele+'",_self=cmd)'          ],
              [ 0, '', ''],
              [ 1, 'transparent surface'   , 'preset.ligand_sites_trans("'+sele+'",_self=cmd)'          ],
              [ 1, 'transparent (better)'   , 'preset.ligand_sites_trans_hq("'+sele+'",_self=cmd)'          ],
              [ 0, '', ''],
              [ 1, 'dot surface'   , 'preset.ligand_sites_dots("'+sele+'",_self=cmd)'          ],
              [ 0, '', ''],
              [ 1, 'mesh surface'   , 'preset.ligand_sites_mesh("'+sele+'",_self=cmd)'          ]]

def presets(self_cmd, sele):
    return [[ 2, 'Preset:'       ,''                        ],
              [ 1, 'classified', 'preset.classified("'+sele+'",_self=cmd)' ],
              [ 0, '', '' ],
              [ 1, 'simple'   ,'preset.simple("'+sele+'",_self=cmd)'          ],
              [ 1, 'simple (no solvent)'   ,'preset.simple_no_solv("'+sele+'",_self=cmd)'          ],
              [ 1, 'ball and stick' , 'preset.ball_and_stick("'+sele+'",_self=cmd)' ],
              [ 1, 'b factor putty' , 'preset.b_factor_putty("'+sele+'",_self=cmd)' ],
              [ 1, 'technical'   , 'preset.technical("'+sele+'",_self=cmd)'          ],
              [ 1, 'ligands'   , 'preset.ligands("'+sele+'",_self=cmd)'          ],
              [ 1, 'ligand sites'   , preset_ligand_sites(self_cmd, sele)         ],
              [ 1, 'pretty ', 'preset.pretty("'+sele+'",_self=cmd)'          ],
              [ 1, 'pretty (with solvent)'     , 'preset.pretty_solv("'+sele+'",_self=cmd)'          ],
              [ 1, 'publication '   , 'preset.publication("'+sele+'",_self=cmd)'          ],
              [ 1, 'publication (with solvent)'   , 'preset.pub_solv("'+sele+'",_self=cmd)'          ],
              [ 0, ''               ,''                             ],
              [ 1, 'protein interface', 'preset.interface("'+sele+'",_self=cmd)' ],
              [ 0, '', '' ],
              [ 1, 'default'   ,'preset.default("'+sele+'",_self=cmd)'          ],
              ]

def hydrogens(self_cmd, sele):
   return [[ 2, 'Hydrogens:'       ,''                        ],
           [ 1, 'add'   ,'cmd.h_add("'+sele+'");cmd.sort("'+sele+' extend 1")' ],
           [ 1, 'add polar' , 'cmd.h_add("'+sele+' & (don.|acc.)");cmd.sort("'+sele+' extend 1")' ],
           [ 1, 'remove'   ,'cmd.remove("('+sele+') and hydro")'          ],
           [ 1, 'remove nonpolar' , 'cmd.remove("'+sele+' & hydro & not nbr. (don.|acc.)")' ],
           ]

def state(self_cmd, sele):
    return [[ 2, 'State:'       ,''                        ],
              [ 1, 'freeze'  ,'cmd.set("state",cmd.get_state(),"'+sele+'")'        ],
              [ 1, 'all states'  ,'cmd.set("state",0,"'+sele+'")' ],
              [ 1, 'thaw'  ,
                  'cmd.unset("all_states","'+sele+'");'
                  'cmd.unset("state","'+sele+'")'        ],
              [ 0, '', '' ],
              [ 1, 'split'   ,'cmd.split_states("'+sele+'")' ],
              ]

def movement(self_cmd, sele):
    return [[ 2, 'Movement:'       ,''                        ],
              [ 1, 'protect'   ,'cmd.protect("'+sele+'")'          ],
              [ 1, 'deprotect'   ,'cmd.deprotect("'+sele+'")'          ],
              ]

def sequence(self_cmd, sele):
    return [[ 2, 'Sequence:'       ,''                        ],
              [ 1, 'include'   ,'cmd.set("seq_view","on","'+sele+'")'          ],
              [ 1, 'exclude'   ,'cmd.set("seq_view","off","'+sele+'")'          ],
              [ 0, ''               ,''                             ],
              [ 1, 'default'   ,'cmd.unset("seq_view","'+sele+'")'          ],
              ]

def masking(self_cmd, sele):
    return [[ 2, 'Masking:'       ,''                        ],
              [ 1, 'mask'   ,'cmd.mask("'+sele+'")'          ],
              [ 1, 'unmask'   ,'cmd.unmask("'+sele+'")'          ],
              ]

def compute(self_cmd, sele):
    return [[ 2, 'Compute:', '' ],
            [ 1, 'atom count'   ,'cmd.count_atoms("'+sele+'",quiet=0)'          ],
            [ 1, 'charges',
              [[ 2, 'Charge:', ''],
               [ 1, 'formal charge sum'   ,'util.sum_formal_charges("'+sele+'",quiet=0,_self=cmd)'    ],
               [ 1, 'partial charges sum'   ,'util.sum_partial_charges("'+sele+'",quiet=0,_self=cmd)' ],
               ]],
            [ 1, 'surface area',
              [[ 2, 'Surface Area Type:', ''],
               [ 1, 'molecular', 'util.get_area("'+sele+'", -1, 0, quiet=0,_self=cmd)' ],
               [ 1, 'solvent accessible', 'util.get_sasa("'+sele+'",quiet=0,_self=cmd)' ],
               [ 1, 'per residue (relative sol. acc.)', 'cmd.get_sasa_relative("'+sele+'",quiet=0,_self=cmd)' ],
               ]],
            [ 1, 'molecular weight',
              [[ 2, 'Molecular Weight:', ''],
               [ 1, 'explicit', 'util.compute_mass("'+sele+'",implicit=False,quiet=0,_self=cmd)' ],
               [ 1, 'with missing hydrogens', 'util.compute_mass("'+sele+'",implicit=True,quiet=0,_self=cmd)' ],
             ]],
            ]

def vacuum(self_cmd, sele):
    return [[ 2, 'Vacuum Electrostatics:'       ,''                        ],
#              [ 2, '\\955WARNING:\\595 Unvalidated and experimental code!', '' ],
              [ 1, 'protein contact potential (local)', 'util.protein_vacuum_esp("'+sele+'",mode=2,quiet=0,_self=cmd)'          ],
#           [ 1, 'protein surface potential (absolute)', 'util.protein_vacuum_esp("'+sele+'",mode=0,quiet=0,_self=cmd)'          ],
#           [ 1, 'protein surface potential (relative)', 'util.protein_vacuum_esp("'+sele+'",mode=1,quiet=0,_self=cmd)'          ],
              [ 2, '\\955NOTE:\\559 Due to short cutoffs, truncations, and', ''],
              [ 2, '\\559lack of solvent "screening", these computed ', ''],
              [ 2, '\\559potentials are only qualitatively useful.', ''],
              [ 2, '\\559Please view with skepticism!', '' ],
              ]

def symmetry(self_cmd, sele):
    return [[ 2, 'Symmetry Mates:'       ,''                        ],
              [ 2, '\\955 +/- one unit cell and...', '' ],
              [ 1, 'within 4 A', 'cmd.symexp("'+sele+'_","'+sele+'","'+sele+'",cutoff=4,segi=1)'          ],
              [ 1, 'within 5 A', 'cmd.symexp("'+sele+'_","'+sele+'","'+sele+'",cutoff=5,segi=1)'          ],
              [ 1, 'within 6 A', 'cmd.symexp("'+sele+'_","'+sele+'","'+sele+'",cutoff=6,segi=1)'          ],
              [ 1, 'within 8 A', 'cmd.symexp("'+sele+'_","'+sele+'","'+sele+'",cutoff=8,segi=1)'          ],
              [ 1, 'within 12 A', 'cmd.symexp("'+sele+'_","'+sele+'","'+sele+'",cutoff=12,segi=1)'          ],
              [ 1, 'within 20 A', 'cmd.symexp("'+sele+'_","'+sele+'","'+sele+'",cutoff=20,segi=1)'          ],
              [ 1, 'within 50 A', 'cmd.symexp("'+sele+'_","'+sele+'","'+sele+'",cutoff=50,segi=1)'          ],
              [ 1, 'within 100 A', 'cmd.symexp("'+sele+'_","'+sele+'","'+sele+'",cutoff=100,segi=1)'          ],
              [ 1, 'within 250 A', 'cmd.symexp("'+sele+'_","'+sele+'","'+sele+'",cutoff=250,segi=1)'          ],
              [ 1, 'within 1000 A', 'cmd.symexp("'+sele+'_","'+sele+'","'+sele+'",cutoff=1000,segi=1)'          ]]

def mol_assign(self_cmd, sele):
    return [[ 2, 'Assign:'       ,''                        ],
              [ 1, 'Amber 99 atomic properties',  'util.assign_amber99("'+sele+'",_self=cmd)' ],
              ]

def selection(self_cmd, sele):
    return [[ 2, 'Selections:', '' ],
              [ 1, 'all', 'cmd.select("'+sele+'_all","'+sele+'")'],
              [ 1, 'polymer', 'cmd.select("'+sele+'_polymer","('+sele+') and polymer")'],
              [ 1, 'organic', 'cmd.select("'+sele+'_organic","('+sele+') and organic")'],
              [ 1, 'solvent', 'cmd.select("'+sele+'_solvent","('+sele+') and solvent")'],
              [ 1, 'polar hydrogens', 'cmd.select("'+sele+'_polar_h","('+sele+') and (e. H and bound_to e. S+O+N)")'],
              [ 1, 'non-polar hydrogens', 'cmd.select("'+sele+'_npolar_h","('+sele+') and (e. H and (not bound_to e. S+O+N))")'],
              [ 1, 'donors', 'cmd.select("'+sele+'_donors","('+sele+') and hbd")'],
              [ 1, 'acceptors', 'cmd.select("'+sele+'_acceptors","('+sele+') and hba")'],
              [ 1, 'surface atoms', 'util.find_surface_atoms(sele="'+sele+'", _self=cmd)' ],
              [ 1, 'C-alphas', 'cmd.select("'+sele+'_calpha","bycalpha ('+sele+')")'],
              ]

def mol_generate(self_cmd, sele):
    return [[ 2, 'Generate:'       ,''                        ],
              [ 1, 'selection', selection(self_cmd, sele) ],
              [ 1, 'symmetry mates', symmetry(self_cmd, sele) ],
              [ 1, 'vacuum electrostatics', vacuum(self_cmd, sele) ],
#           [ 1, 'assign', mol_assign(self_cmd, sele) ],
              ]

def invert(self_cmd, sele):
    return [[ 2, 'Invert:'       ,''                        ],
              [ 1, 'within object(s)'     ,'cmd.select("'+sele+'","((byobj '+sele+') and not '+sele+')",enable=1)'    ],
              [ 1, 'within segment(s)'     ,'cmd.select("'+sele+'","((byseg '+sele+') and not '+sele+')",enable=1)'    ],
              [ 1, 'within chain(s)'     ,'cmd.select("'+sele+'","((bychain '+sele+') and not '+sele+')",enable=1)'    ],
              [ 1, 'within residue(s)'   ,'cmd.select("'+sele+'","((byres '+sele+') and not '+sele+')",enable=1)'    ],
              [ 0, ''               ,''                             ],
              [ 1, 'within molecule(s)'     ,'cmd.select("'+sele+'","((bymol '+sele+') and not '+sele+')",enable=1)'    ],
              [ 0, ''               ,''                             ],
              [ 1, 'within any'     ,'cmd.select("'+sele+'","(not '+sele+')",enable=1)'    ],
              ]

def complete(self_cmd, sele):
    return [[ 2, 'Complete:'       ,''                        ],

              [ 1, 'residues'  ,'cmd.select("'+sele+'","(byres '+sele+')",enable=1)'      ],
              [ 1, 'chains'  ,'cmd.select("'+sele+'","(bychain '+sele+')",enable=1)'      ],
              [ 1, 'segments'  ,'cmd.select("'+sele+'","(byseg '+sele+')",enable=1)'      ],
              [ 1, 'objects'  ,'cmd.select("'+sele+'","(byobj '+sele+')",enable=1)'      ],
              [ 0, ''               ,''                             ],
              [ 1, 'molecules'  ,'cmd.select("'+sele+'","(bymol '+sele+')",enable=1)'      ],
              [ 0, ''               ,''                             ],
              [ 1, 'C-alphas'  ,'cmd.select("'+sele+'","(bycalpha '+sele+')",enable=1)'      ],
              ]

def modify_by_object(self_cmd, sele, op):
    list = self_cmd.get_names("public_objects",1)[0:25] # keep this practical
    list = [x for x in list if self_cmd.get_type(x)=="object:molecule"]
    result = [[ 2, 'Object:', '']]
    for a in list:
        if a!=sele:
            result.append([1,a,
                                'cmd.select("'+sele+'","('+sele+') '+op+' ('+a+')",enable=1)'])
    return result

def modify_by_sele(self_cmd, sele, op):
    list = self_cmd.get_names("public_selections",0)[0:25] # keep this practical
    result = [[ 2, 'Selection:', '']]
    for a in list:
        if a!=sele:
            result.append([1,a, 'cmd.select("'+sele+'","('+sele+') '+op+' ('+a+')",enable=1)'])
    return result

def restrict(self_cmd, sele):
    return [[ 2, 'Restrict:'       ,''                        ],
            [ 1, 'to object'   , modify_by_object(self_cmd, sele,'and') ],
            [ 1, 'to selection' , modify_by_sele(self_cmd, sele,'and') ],
            [ 0, ''               ,''                             ],
            [ 1, 'to visible'   , 'cmd.select("'+sele+'","('+sele+') and vis",enable=1)'],
            [ 0, ''               ,''                             ],
            [ 1, 'to polymer'   , 'cmd.select("'+sele+'","('+sele+') and polymer",enable=1)'],
            [ 1, 'to solvent'   , 'cmd.select("'+sele+'","('+sele+') and solvent",enable=1)'],
            [ 1, 'to organic'   , 'cmd.select("'+sele+'","('+sele+') and organic",enable=1)'],
            [ 1, 'to inorganic'   , 'cmd.select("'+sele+'","('+sele+') and inorganic",enable=1)'],
            ]

def include(self_cmd, sele):
    return [[ 2, 'Include:'       ,''                        ],
              [ 1, 'object'   , modify_by_object(self_cmd, sele,'or') ],
              [ 1, 'selection' , modify_by_sele(self_cmd, sele,'or') ],
              [ 0, ''               ,''                             ],
              [ 1, 'visible'   , 'cmd.select("'+sele+'","('+sele+') or vis",enable=1)'],
              ]

def exclude(self_cmd, sele):
    return [[ 2, 'Exclude:'       ,''                        ],
            [ 1, 'object'   , modify_by_object(self_cmd, sele,'and not') ],
            [ 1, 'selection' , modify_by_sele(self_cmd, sele,'and not') ],
            [ 0, ''               ,''                             ],
            [ 1, 'polymer'   , 'cmd.select("'+sele+'","('+sele+') and not organic",enable=1)'],
            [ 1, 'solvent'   , 'cmd.select("'+sele+'","('+sele+') and not solvent",enable=1)'],
            [ 1, 'organic'   , 'cmd.select("'+sele+'","('+sele+') and not organic",enable=1)'],
            [ 1, 'inorganic' , 'cmd.select("'+sele+'","('+sele+') and not organic",enable=1)'],
            ]

def expand(self_cmd, sele):
    return [[ 2, 'Expand:'       ,''                        ],
              [ 1, 'by 4 A'  ,'cmd.select("'+sele+'","('+sele+' expand 4)",enable=1)' ],
              [ 1, 'by 5 A'  ,'cmd.select("'+sele+'","('+sele+' expand 5)",enable=1)' ],
              [ 1, 'by 6 A'  ,'cmd.select("'+sele+'","('+sele+' expand 6)",enable=1)' ],
              [ 1, 'by 8 A'  ,'cmd.select("'+sele+'","('+sele+' expand 8)",enable=1)' ],
              [ 1, 'by 12 A'  ,'cmd.select("'+sele+'","('+sele+' expand 12)",enable=1)' ],
              [ 1, 'by 20 A'  ,'cmd.select("'+sele+'","('+sele+' expand 20)",enable=1)' ],
              [ 0, ''               ,''                             ],
              [ 1, 'by 4 A, residues'  ,'cmd.select("'+sele+'","(byres ('+sele+' expand 4))",enable=1)' ],
              [ 1, 'by 5 A, residues'  ,'cmd.select("'+sele+'","(byres ('+sele+' expand 5))",enable=1)' ],
              [ 1, 'by 6 A, residues'  ,'cmd.select("'+sele+'","(byres ('+sele+' expand 6))",enable=1)' ],
              [ 1, 'by 8 A, residues'   ,'cmd.select("'+sele+'","(byres ('+sele+' expand 8))",enable=1)' ],
              [ 1, 'by 12 A, residues'  ,'cmd.select("'+sele+'","(byres ('+sele+' expand 12))",enable=1)' ],
              [ 1, 'by 20 A, residues'  ,'cmd.select("'+sele+'","(byres ('+sele+' expand 20))",enable=1)' ],
              ]

def around(self_cmd, sele):
    return [[ 2, 'Around:'       ,''                        ],
              [ 1, 'atoms within 4 A'  ,'cmd.select("'+sele+'","('+sele+' around 4)",enable=1)' ],
              [ 1, 'atoms within 5 A'  ,'cmd.select("'+sele+'","('+sele+' around 5)",enable=1)' ],
              [ 1, 'atoms within 6 A'  ,'cmd.select("'+sele+'","('+sele+' around 6)",enable=1)' ],
              [ 1, 'atoms within 8 A'  ,'cmd.select("'+sele+'","('+sele+' around 8)",enable=1)' ],
              [ 1, 'atoms within 12 A'  ,'cmd.select("'+sele+'","('+sele+' around 12)",enable=1)' ],
              [ 1, 'atoms within 20 A'  ,'cmd.select("'+sele+'","('+sele+' around 20)",enable=1)' ],
              [ 0, ''               ,''                             ],
              [ 1, 'residues within 4 A'  ,'cmd.select("'+sele+'","(byres ('+sele+' around 4))",enable=1)' ],
              [ 1, 'residues within 5 A'  ,'cmd.select("'+sele+'","(byres ('+sele+' around 5))",enable=1)' ],
              [ 1, 'residues within 6 A'  ,'cmd.select("'+sele+'","(byres ('+sele+' around 6))",enable=1)' ],
              [ 1, 'residues within 8 A'  ,'cmd.select("'+sele+'","(byres ('+sele+' around 8))",enable=1)' ],
              [ 1, 'residues within 12 A'  ,'cmd.select("'+sele+'","(byres ('+sele+' around 12))",enable=1)' ],
              [ 1, 'residues within 20 A'  ,'cmd.select("'+sele+'","(byres ('+sele+' around 20))",enable=1)' ],
              ]

def extend(self_cmd, sele):
    return [[ 2, 'Extend:'       ,''                        ],
              [ 1, 'by 1 bond'  ,'cmd.select("'+sele+'","('+sele+' extend 1)",enable=1)' ],
              [ 1, 'by 2 bonds'  ,'cmd.select("'+sele+'","('+sele+' extend 2)",enable=1)' ],
              [ 1, 'by 3 bonds'  ,'cmd.select("'+sele+'","('+sele+' extend 3)",enable=1)' ],
              [ 1, 'by 4 bonds'  ,'cmd.select("'+sele+'","('+sele+' extend 4)",enable=1)' ],
              [ 1, 'by 5 bonds'  ,'cmd.select("'+sele+'","('+sele+' extend 5)",enable=1)' ],
              [ 1, 'by 6 bonds'  ,'cmd.select("'+sele+'","('+sele+' extend 6)",enable=1)' ],
              [ 0, ''               ,''                             ],
              [ 1, 'by 1 bond, residues'  ,'cmd.select("'+sele+'","(byres ('+sele+' extend 1))",enable=1)' ],
              [ 1, 'by 2 bonds, residues'  ,'cmd.select("'+sele+'","(byres ('+sele+' extend 2))",enable=1)' ],
              [ 1, 'by 3 bonds, residues'   ,'cmd.select("'+sele+'","(byres ('+sele+' extend 3))",enable=1)' ],
              [ 1, 'by 4 bonds, residues'  ,'cmd.select("'+sele+'","(byres ('+sele+' extend 4))",enable=1)' ],
              [ 1, 'by 5 bonds, residues'  ,'cmd.select("'+sele+'","(byres ('+sele+' extend 5))",enable=1)' ],
              [ 1, 'by 6 bonds, residues'  ,'cmd.select("'+sele+'","(byres ('+sele+' extend 6))",enable=1)' ],
              ]

def polar(self_cmd, sele):
    return [[ 2, 'Polar Contacts:', ''],
              [ 1, 'within selection'  ,
                 'cmd.dist("'+sele+'_polar_conts","'+sele+'","'+sele+'",quiet=1,mode=2,label=0,reset=1);cmd.enable("'+sele+'_polar_conts")'],
              [ 1, 'involving side chains'  ,
                 'cmd.dist("'+sele+'_polar_conts","('+sele+')","('+sele+
                 ') & sc.",quiet=1,mode=2,label=0,reset=1);cmd.enable("'+sele+'_polar_conts")'],
              [ 1, 'involving solvent'  ,
                 'cmd.dist("'+sele+'_polar_conts","('+sele+') and solvent","('+sele+
                 ') and not (solvent)",quiet=1,mode=2,label=0,reset=1);cmd.enable("'+sele+'_polar_conts")'],
              [ 1, 'excluding solvent'  ,
                 'cmd.dist("'+sele+'_polar_conts","('+sele+') and not (solvent)","('+sele+
                 ') and not (solvent)",quiet=1,mode=2,label=0,reset=1);cmd.enable("'+sele+'_polar_conts")'],
              [ 1, 'excluding main chain'  ,
                 'cmd.dist("'+sele+'_polar_conts","('+sele+') & !bb.","('+sele+') & !bb.",'
                 'quiet=1,mode=2,label=0,reset=1);cmd.enable("'+sele+'_polar_conts")'],
              [ 1, 'excluding intra-main chain'  ,
                 'cmd.dist("'+sele+'_polar_conts","('+sele+')","('+sele+
                 ') & !bb.",quiet=1,mode=2,label=0,reset=1);cmd.enable("'+sele+'_polar_conts")'],
              [ 1, 'just intra-side chain'  ,
                 'cmd.dist("'+sele+'_polar_conts","('+sele+') & sc.","('+sele+') & sc.",'
                 'quiet=1,mode=2,label=0,reset=1);cmd.enable("'+sele+
                 '_polar_conts")'],
              [ 1, 'just intra-main chain'  ,
                 'cmd.dist("'+sele+'_polar_conts","('+sele+') & bb.","('+sele+') & bb.",'
                 'quiet=1,mode=2,label=0,reset=1);cmd.enable("'+sele+
                 '_polar_conts")'],
              [ 0, '', '' ],
              [ 1, 'to other atoms in object',
                 'cmd.dist("'+sele+'_polar_conts","('+sele+')","(byobj ('+sele+')) and (not ('+sele+
                 '))",quiet=1,mode=2,label=0,reset=1);cmd.enable("'+sele+'_polar_conts")'],
              [ 1, 'to others excluding solvent',
                 'cmd.dist("'+sele+'_polar_conts","('+sele+') and not solvent","(byobj ('+sele+')) and (not ('+sele+
                 ')) and (not solvent)",quiet=1,mode=2,label=0,reset=1);cmd.enable("'+sele+'_polar_conts")'],
              [ 1, 'to any atoms',
                 'cmd.dist("'+sele+'_polar_conts","('+sele+')","(not '+sele+
                 ')",quiet=1,mode=2,label=0,reset=1);cmd.enable("'+sele+'_polar_conts")'],
              [ 1, 'to any excluding solvent',
                 'cmd.dist("'+sele+'_polar_conts","('+sele+') and not solvent","(not ('+sele+
                 ')) and not solvent",quiet=1,mode=2,label=0,reset=1);cmd.enable("'+sele+'_polar_conts")'],
              [ 0, '', '' ],
              [ 1, 'between chains',
                 'util.interchain_distances("'+sele+'_interchain_polar","'+sele+'",mode=2)'],
              ]


def find(self_cmd, sele):
    return [[ 2, 'Find:', ''],
              [ 1, 'polar contacts', polar(self_cmd, sele) ],
              [ 1, 'any contacts', [[ 2, 'Any Contacts:', '']] + [
                  [ 1, 'between chains within %.1fA' % d, 'util.interchain_distances("'+sele+'_interchain_any","'+sele+'",cutoff=%f)' % d]
                  for d in (3.0, 3.5, 4.0)
              ]],
              [ 1, 'pi interactions', [[ 2, 'Pi Interactions:', '']] + [
                  [1, 'all', 'cmd.pi_interactions("'+sele+'_pi_interactions","'+sele+'",reset=1)'],
                  [1, 'pi-pi', 'cmd.distance("'+sele+'_pi_pi","'+sele+'","same",reset=1,mode=6)'],
                  [1, 'pi-cation', 'cmd.distance("'+sele+'_pi_cation","'+sele+'","same",reset=1,mode=7)'],
              ]],
              ]

def align_to_object(self_cmd, sele):
    list = self_cmd.get_names("public_objects",1)[0:25] # keep this practical
    list = [x for x in list if self_cmd.get_type(x)=="object:molecule"]
    result = [[ 2, 'Object:', '']]
    for a in list:
        if a!=sele:
            result.append([1,a,
                           'cmd.align("polymer and name CA and ('+sele+')",'+
                           '"polymer and name CA and ('+a+')",quiet=0,'+
                           'object="aln_%s_to_%s",reset=1)'%(sele,a)])
    return result

def align_to_sele(self_cmd, sele):
    list = self_cmd.get_names("public_selections",0)[0:25] # keep this practical
    result = [[ 2, 'Selection:', '']]
    for a in list:
        if a!= sele:
            result.append([1,a,
                           'cmd.align("polymer and name CA and ('+sele+')",'+
                           '"polymer and name CA and ('+a+')",quiet=0,'+
                           'object="aln_%s_to_%s",reset=1)'%(sele,a)])
    return result

def mat_tran(self_cmd, sele, direction=0):
    list = self_cmd.get_names("public_objects",1)[0:25] # keep this practical
    list = [x for x in list if self_cmd.get_type(x)!="object:ramp"]
    result = [[ 2, 'Object:', '']]
    for a in list:
        if a!=sele:
            if direction:
                result.append([1,a,
                                    'cmd.matrix_transfer("'+a+'","'+sele+'");'])
            else:
                result.append([1,a,
                                    'cmd.matrix_transfer("'+sele+'","'+a+'");'])
    return result


def sele_align(self_cmd, sele):
    return [[ 2, 'Align:', ''],
              [ 1, 'to molecule (*/CA)', align_to_object(self_cmd, sele) ],
              [ 1, 'to selection (*/CA)', align_to_sele(self_cmd, sele) ],
              [ 0, '', None ],
              [ 1, 'enabled to this (*/CA)', 'util.mass_align("'+sele+'",1,_self=cmd)' ],
              [ 1, 'all to this (*/CA)', 'util.mass_align("'+sele+'",0,_self=cmd)' ],
              [ 0, '', None ],
              [ 1, 'states (*/CA)', 'cmd.intra_fit("('+sele+') and name CA")' ],
              [ 1, 'states', 'cmd.intra_fit("'+sele+'")' ],
              ]

def mol_align(self_cmd, sele):
    return sele_align(self_cmd, sele) + [
              [ 0, '', None ],
              [ 1, 'matrix from', mat_tran(self_cmd, sele,1) ],
              [ 1, 'matrix to', mat_tran(self_cmd, sele,0) ],
              [ 1, 'matrix reset', 'cmd.matrix_reset("'+sele+'")'],
              ]

def modify_sele(self_cmd, sele):
    return [[ 2, 'Modify:', ''],
              [ 1, 'around'         , around(self_cmd, sele)         ],
              [ 1, 'expand'         , expand(self_cmd, sele)         ],
              [ 1, 'extend'         , extend(self_cmd, sele)         ],
              [ 1, 'invert'         , invert(self_cmd, sele)         ],
              [ 1, 'complete'       , complete(self_cmd, sele)         ],
              [ 1, 'restrict'       , restrict(self_cmd, sele)       ],
              [ 1, 'include'        , include(self_cmd, sele)       ],
              [ 1, 'exclude'        , exclude(self_cmd, sele)       ]]


def copy_to(self_cmd, sele):
    onames = self_cmd.get_object_list('enabled')[:25]  # keep this practical
    selected = self_cmd.get_object_list(sele)
    return [[ 2, 'Copy To:', '' ],
              [ 1, 'new', 'cmd.create(None,"' + sele + '",zoom=0)' ],
              [ 0, '', '' ],
              ] + [
                  [ 1, oname, 'cmd.copy_to("' + oname + '","' + sele + '",zoom=0,quiet=0)' ]
                  for oname in onames
                  if oname not in selected
              ]

def move_to_group(self_cmd, sele):
    gnames = self_cmd.get_names_of_type('object:group')
    return [
        [ 2, 'Move to Group:', '' ],
        [ 1, 'new', 'cmd.group(cmd.get_unused_name("group"),"' + sele + '")' ],
        [ 1, 'ungroup', 'cmd.ungroup("' + sele + '")' ],
        [ 0, '', '' ],
    ] + [
        [ 1, gname, 'cmd.group("' + gname + '","' + sele + '",quiet=0)' ]
        for gname in gnames if gname != sele
    ]

def sele_action(self_cmd, sele):
    return [[ 2, 'Action:'       ,''                        ],
              [ 1, del_col + 'delete selection', 'cmd.delete("'+sele+'")'          ],
              [ 1, 'rename selection', 'cmd.wizard("renaming","'+sele+'")'          ],
              [ 0, ''               ,''                             ],
              [ 1, 'zoom'           ,'cmd.zoom("'+sele+'",animate=-1)'            ],
              [ 1, 'orient'         ,'cmd.orient("'+sele+'",animate=-1)'          ],
              [ 1, 'center'         ,'cmd.center("'+sele+'",animate=-1)'            ],
              [ 1, 'origin'         ,'cmd.origin("'+sele+'")'          ],
              [ 0, ''               ,''                             ],
              [ 1, 'drag coordinates'     , 'cmd.drag("'+sele+'")'    ],
              [ 1, 'clean'       , 'cmd.clean("'+sele+'")'    ],
              [ 0, ''               ,''                             ],
              [ 1, 'modify', modify_sele(self_cmd, sele) ],
              [ 1, 'preset'         ,presets(self_cmd, sele)         ],
              [ 1, 'find', find(self_cmd, sele) ],
              [ 1, 'align', sele_align(self_cmd, sele) ],
              [ 0, ''               ,''                             ],
              [ 1, rem_col + 'remove atoms'   ,'cmd.remove("'+sele+'");cmd.delete("'+sele+'")'          ],
              [ 1, 'hydrogens' , hydrogens(self_cmd, sele) ],
              [ 0, ''          ,''                                              ],
              [ 1, 'duplicate'      ,'cmd.select(None,"'+sele+'")'          ], # broken...
              [ 1, 'copy to object' , lambda: copy_to(self_cmd, sele) ],
              [ 1, 'extract object' ,'cmd.extract(None,"'+sele+'",zoom=0)' ],
              [ 0, ''          ,''                                  ],
              [ 1, 'masking'        , masking(self_cmd, sele)         ],
              [ 1, 'movement'       , movement(self_cmd, sele)         ],
              [ 1, 'compute'        , compute(self_cmd, sele) ],
              ]


def sele_action2(self_cmd, sele):
    return [[ 2, 'Action:'       ,''                        ],
              [ 1, del_col + 'delete selection', 'cmd.delete("'+sele+'")'          ],
              [ 1, 'rename selection', 'cmd.wizard("renaming","'+sele+'")'          ],
              [ 0, ''               ,''                             ],
              [ 1, 'preset'         ,presets(self_cmd, sele)         ],
              [ 1, 'find', find(self_cmd, sele) ],
              [ 0, ''               ,''                             ],
              [ 1, rem_col + 'remove atoms'   ,'cmd.remove("'+sele+'");cmd.delete("'+sele+'")'          ],
              [ 0, ''               ,''                             ],
              [ 1, 'around'         , around(self_cmd, sele)         ],
              [ 1, 'expand'         , expand(self_cmd, sele)         ],
              [ 1, 'extend'         , extend(self_cmd, sele)         ],
              [ 1, 'invert'         , invert(self_cmd, sele)         ],
              [ 1, 'complete'       , complete(self_cmd, sele)         ],
              [ 0, ''          ,''                                              ],
              [ 1, 'duplicate selection'      ,'cmd.select(None,"'+sele+'")'          ],
              [ 1, 'copy to object' , lambda: copy_to(self_cmd, sele) ],
              [ 1, 'extract object' ,'cmd.extract(None,"'+sele+'",zoom=0)' ],
            [ 0, ''          ,''                                  ],
              [ 1, 'masking'      , masking(self_cmd, sele)         ],
              [ 1, 'movement'       , movement(self_cmd, sele)         ],
              [ 1, 'compute'        , compute(self_cmd, sele)         ],
              ]


def group_action(self_cmd, sele):
    return [[ 2, 'Action:'     , ''                       ],
            [ 1, 'zoom'         , 'cmd.zoom("'+sele+'",animate=-1)'      ],
            [ 1, 'orient'       , 'cmd.orient("'+sele+'",animate=-1)'    ],
            [ 1, 'center'         ,'cmd.center("'+sele+'",animate=-1)'            ],
            [ 1, 'origin'       , 'cmd.origin("'+sele+'")'    ],
            [ 0, ''               ,''                             ],
            [ 1, 'drag' , 'cmd.drag("'+sele+'")' ],
            [ 1, 'reset' , 'cmd.reset(object="'+sele+'")' ],
            [ 0, ''               ,''                             ],
            [ 1, 'preset'  ,   presets(self_cmd, sele)       ],
            [ 1, 'find',     find(self_cmd, sele) ],
            [ 1, 'align',     mol_align(self_cmd, sele) ],
            [ 1, 'generate'  ,   mol_generate(self_cmd, sele)       ],
            [ 0, ''               ,''                             ],
            [ 1, 'assign sec. struc.'  ,'cmd.dss("'+sele+'")'        ],
            [ 0, ''             , ''                       ],
            [ 1, 'rename group', 'cmd.wizard("renaming","'+sele+'")'          ],
            [ 1, 'group' , lambda: move_to_group(self_cmd, sele) ],
            [ 1, del_col + 'delete group', 'cmd.delete("'+sele+'")'    ],
            [ 0, ''          ,''                                              ],
            [ 1, 'hydrogens' , hydrogens(self_cmd, sele)    ],
            [ 1, rem_col + 'remove waters', 'cmd.remove("(solvent and ('+sele+'))")'     ],
            [ 0, ''          ,''                                              ],
            [ 1, 'state'          , state(self_cmd, sele)         ],
            [ 1, 'masking'        , masking(self_cmd, sele)         ],
            [ 1, 'sequence'       , sequence(self_cmd, sele)         ],
            [ 1, 'movement'       , movement(self_cmd, sele)         ],
            [ 1, 'compute'        , compute(self_cmd, sele)         ],
            ]

def mol_action(self_cmd, sele):
    return [[ 2, 'Action:'     , ''                       ],
            [ 1, 'zoom'         , 'cmd.zoom("'+sele+'",animate=-1)'      ],
            [ 1, 'orient'       , 'cmd.orient("'+sele+'",animate=-1)'    ],
            [ 1, 'center'         ,'cmd.center("'+sele+'",animate=-1)'            ],
            [ 1, 'origin'       , 'cmd.origin("'+sele+'")'    ],
            [ 0, ''               ,''                             ],
            [ 1, 'drag matrix' , 'cmd.drag("'+sele+'")' ],
            [ 1, 'reset matrix' , 'cmd.reset(object="'+sele+'")' ],
            [ 0, ''               ,''                             ],
            [ 1, 'drag coordinates' , 'cmd.drag("('+sele+')")' ],
            [ 1, 'clean'       , 'cmd.clean("'+sele+'")'    ],
            [ 0, ''          ,''                                              ],
            [ 1, 'preset'  ,   presets(self_cmd, sele)       ],
            [ 1, 'find',     find(self_cmd, sele) ],
            [ 1, 'align',     mol_align(self_cmd, sele) ],
            [ 1, 'generate'  ,   mol_generate(self_cmd, sele)       ],
            [ 0, ''               ,''                             ],
            [ 1, 'assign sec. struc.'  ,'cmd.dss("'+sele+'")'        ],
            [ 0, ''             , ''                       ],
            [ 1, 'rename object', 'cmd.wizard("renaming","'+sele+'")'          ],
            [ 1, 'copy to object' , lambda: copy_to(self_cmd, sele) ],
            [ 1, 'group' , lambda: move_to_group(self_cmd, sele) ],
            [ 1, del_col + 'delete object', 'cmd.delete("'+sele+'")'    ],
            [ 0, ''          ,''                                              ],
            [ 1, 'hydrogens' , hydrogens(self_cmd, sele)    ],
            [ 1, rem_col + 'remove waters', 'cmd.remove("(solvent and ('+sele+'))")'     ],
            [ 0, ''          ,''                                              ],
              [ 1, 'state'          , state(self_cmd, sele)         ],
              [ 1, 'masking'        , masking(self_cmd, sele)         ],
              [ 1, 'sequence'       , sequence(self_cmd, sele)         ],
              [ 1, 'movement'       , movement(self_cmd, sele)         ],
              [ 1, 'compute'        , compute(self_cmd, sele)         ],
              ]

def slice_action(self_cmd, sele):
    return [[ 2, 'Action:'     , ''                       ],
              [ 1, 'zoom'         , 'cmd.zoom("'+sele+'",animate=-1)'      ],
              [ 1, 'center'       , 'cmd.center("'+sele+'",animate=-1)'    ],
              [ 1, 'origin'       , 'cmd.origin("'+sele+'")'    ],
              [ 0, ''             , ''                       ],
              [ 1, 'tracking on' , 'cmd.set("slice_track_camera",1,"'+sele+'")'      ],
              [ 1, 'tracking off' , 'cmd.set("slice_track_camera",0,"'+sele+'")'      ],
              [ 0, ''             , ''                       ],
              [ 1, 'height map on' , 'cmd.set("slice_height_map",1,"'+sele+'")'    ],
              [ 1, 'height map off', 'cmd.set("slice_height_map",0,"'+sele+'")'    ],
              [ 0, ''             , ''                       ],
              [ 1, 'dynamic grid on' , 'cmd.set("slice_dynamic_grid",1,"'+sele+'")'    ],
              [ 1, 'dynamic grid off', 'cmd.set("slice_dynamic_grid",0,"'+sele+'")'    ],
              [ 0, ''             , ''                       ],
              [ 1, 'rename'       , 'cmd.wizard("renaming","'+sele+'")'          ],
              [ 1, 'group' , lambda: move_to_group(self_cmd, sele) ],
              [ 0, ''             , ''                       ],
              [ 1, del_col + 'delete', 'cmd.delete("'+sele+'")'    ],
              ]

def simple_action(self_cmd, sele):
    return [[ 2, 'Action:'     , ''                       ],
            [ 1, 'zoom'         , 'cmd.zoom("'+sele+'",animate=-1)'      ],
            [ 1, 'center'       , 'cmd.center("'+sele+'",animate=-1)'    ],
            [ 1, 'origin'       , 'cmd.origin("'+sele+'")'    ],

            [ 0, ''             , ''                       ],
            [ 1, 'drag'       , 'cmd.drag("'+sele+'")'          ],
            [ 1, 'reset'       , 'cmd.reset(object="'+sele+'")'          ],
            [ 0, ''             , ''                       ],
            [ 1, 'rename'       , 'cmd.wizard("renaming","'+sele+'")'          ],
            [ 1, 'group' , lambda: move_to_group(self_cmd, sele) ],
            [ 0, ''             , ''                       ],
            [ 1, del_col + 'delete', 'cmd.delete("'+sele+'")'    ],
              ]

def iso_with_negative(mapname, suffix, rep, level=1, color='blue'):
    name = mapname + suffix
    return ('cmd.iso%s("%s","%s", %s);\n'
            'cmd.set("%s_negative_visible",1,"%s",quiet=0);\n'
            'cmd.color("%s","%s")' % (rep, name, mapname, level, rep, name, color, name))

def map_mesh(self_cmd, sele):
    return [[ 2, 'Mesh:',  '' ],
            [ 1, '@ level 1.0'         , 'cmd.isomesh("'+sele+'_mesh","'+sele+'",1.0)'      ],
            [ 0, ''             , ''                       ],
            [ 1, '@ level 2.0'         , 'cmd.isomesh("'+sele+'_mesh","'+sele+'",2.0)'      ],
            [ 1, '@ level 3.0'         , 'cmd.isomesh("'+sele+'_mesh","'+sele+'",3.0)'      ],
            [ 0, ''             , ''                       ],
            [ 1, '@ level +/-1.0'      , iso_with_negative(sele, '_mesh', 'mesh')],
            [ 1, '@ level +/-3.0'      , iso_with_negative(sele, '_mesh', 'mesh', 3, 'green')],
            [ 0, ''             , ''                       ],
            [ 1, '@ level 0.0'         , 'cmd.isomesh("'+sele+'_mesh","'+sele+'",0.0)'      ],
            [ 1, '@ level -1.0'         , 'cmd.isomesh("'+sele+'_mesh","'+sele+'",-1.0)'      ],
            [ 1, '@ level -2.0'         , 'cmd.isomesh("'+sele+'_mesh","'+sele+'",-2.0)'      ],
            [ 1, '@ level -3.0'         , 'cmd.isomesh("'+sele+'_mesh","'+sele+'",-3.0)'      ],
            ]

def map_volume(self_cmd, sele):
    from pymol.colorramping import namedramps
    return [[ 2, 'Volume:', ''],
            [ 1, 'default'              , 'cmd.volume("'+sele+'_volume","'+sele+'")'  ],
            [ 0, '', '' ],
        ] + [
            [ 1, p, 'cmd.volume("%s_volume","%s","%s")' % (sele, sele, p) ]
            for p in sorted(namedramps)
        ]

def map_surface(self_cmd, sele):
    return [[ 2, 'Surface:',  '' ],
            [ 1, '@ level 1.0'         , 'cmd.isosurface("'+sele+'_surf","'+sele+'",1.0)'      ],
            [ 0, ''             , ''                       ],
            [ 1, '@ level 2.0'         , 'cmd.isosurface("'+sele+'_surf","'+sele+'",2.0)'      ],
            [ 1, '@ level 3.0'         , 'cmd.isosurface("'+sele+'_surf","'+sele+'",3.0)'      ],
            [ 0, ''             , ''                       ],
            [ 1, '@ level +/-1.0'      , iso_with_negative(sele, '_surf', 'surface')],
            [ 1, '@ level +/-3.0'      , iso_with_negative(sele, '_surf', 'surface', 3, 'green')],
            [ 0, ''             , ''                       ],
            [ 1, '@ level 0.0'         , 'cmd.isosurface("'+sele+'_surf","'+sele+'",0.0)'      ],
            [ 1, '@ level -1.0'         , 'cmd.isosurface("'+sele+'_surf","'+sele+'",-1.0)'      ],
            [ 1, '@ level -2.0'         , 'cmd.isosurface("'+sele+'_surf","'+sele+'",-2.0)'      ],
            [ 1, '@ level -3.0'         , 'cmd.isosurface("'+sele+'_surf","'+sele+'",-3.0)'      ],
            ]

def map_gradient(self_cmd, sele):
    return [[ 2, 'Gradient:',  '' ],
            [ 1, 'default'         , 'cmd.gradient("'+sele+'_grad","'+sele+'");cmd.ramp_new("'+sele+
              '_grad_ramp","'+sele+'");cmd.color("'+sele+'_grad_ramp","'+sele+'_grad");' ]
            ]

def map_slice(self_cmd, sele):
    return [[ 2, 'Slice:',  '' ],
            [ 1, 'default'         , 'cmd.slice_new("'+sele+'_slice","'+sele+'");cmd.ramp_new("'+sele+
              '_slice_ramp","'+sele+'");cmd.color("'+sele+'_slice_ramp","'+sele+'_slice");'+
              'cmd.set("slice_track_camera",1,"'+sele+'_slice");'+
              'cmd.set("two_sided_lighting",1,"'+sele+'_slice");'+
              'cmd.set("ray_interior_color",0,"'+sele+'_slice");'+
              'cmd.set("slice_dynamic_grid",1,"'+sele+'_slice")'],
            ]

def map_action(self_cmd, sele):
    return [[ 2, 'Action:'     , ''                       ],
            [ 1, 'mesh'         , map_mesh(self_cmd, sele)  ],
            [ 1, 'surface'      , map_surface(self_cmd, sele)  ],
            [ 1, 'slice'        , map_slice(self_cmd, sele)  ],
            [ 1, 'gradient'     , map_gradient(self_cmd, sele)  ],
            [ 1, 'volume'       , map_volume(self_cmd, sele)  ],
            [ 0, ''             , ''                       ],
            [ 1, 'zoom'         , 'cmd.zoom("'+sele+'",animate=-1)'      ],
            [ 1, 'center'       , 'cmd.center("'+sele+'",animate=-1)'    ],
            [ 1, 'origin'       , 'cmd.origin("'+sele+'")'    ],
            [ 0, ''             , ''                       ],
            [ 1, 'drag' , 'cmd.drag("'+sele+'")' ],
            [ 1, 'reset'       , 'cmd.reset(object="'+sele+'")'          ],
            [ 0, ''             , ''                       ],
            [ 1, 'matrix_copy'  , mat_tran(self_cmd, sele, 1) ],
            [ 0, ''             , '' ],
            [ 1, 'rename'       , 'cmd.wizard("renaming","'+sele+'")'          ],
            [ 1, 'group' , lambda: move_to_group(self_cmd, sele) ],
            [ 0, ''             , ''                       ],
            [ 1, del_col + 'delete', 'cmd.delete("'+sele+'")'    ],
            ]

def level(self_cmd, sele):
    return [[ 2, 'Level',  '' ],
            [ 1, 'level 5.0'         , 'cmd.isolevel("'+sele+'",5.0)'      ],
            [ 1, 'level 4.0'         , 'cmd.isolevel("'+sele+'",4.0)'      ],
            [ 1, 'level 3.0'         , 'cmd.isolevel("'+sele+'",3.0)'      ],
            [ 1, 'level 2.0'         , 'cmd.isolevel("'+sele+'",2.0)'      ],
            [ 1, 'level 1.5'         , 'cmd.isolevel("'+sele+'",1.5)'      ],
            [ 1, 'level 1.0'         , 'cmd.isolevel("'+sele+'",1.0)'      ],
            [ 1, 'level 0.5'         , 'cmd.isolevel("'+sele+'",0.5)'      ],
            [ 1, 'level 0.0'         , 'cmd.isolevel("'+sele+'",0.0)'      ],
            [ 1, 'level -0.5'         , 'cmd.isolevel("'+sele+'",-0.5)'      ],
            [ 1, 'level -1.0'         , 'cmd.isolevel("'+sele+'",-1.0)'      ],
            [ 1, 'level -1.5'         , 'cmd.isolevel("'+sele+'",-1.5)'      ],
            [ 1, 'level -2.0'         , 'cmd.isolevel("'+sele+'",-2.0)'      ],
            [ 1, 'level -3.0'         , 'cmd.isolevel("'+sele+'",-3.0)'      ],
            [ 1, 'level -4.0'         , 'cmd.isolevel("'+sele+'",-4.0)'      ],
            [ 1, 'level -5.0'         , 'cmd.isolevel("'+sele+'",-5.0)'      ],
            ]

def surface_action(self_cmd, sele):
    return [[ 2, 'Action:'     , ''                       ],
            [ 1, 'level'         , level(self_cmd, sele)  ],
            [ 0, ''             , ''                       ],
            [ 1, 'zoom'         , 'cmd.zoom("'+sele+'",animate=-1)'      ],
            [ 1, 'center'       , 'cmd.center("'+sele+'",animate=-1)'    ],
            [ 1, 'origin'       , 'cmd.origin("'+sele+'")'    ],
            [ 0, ''             , ''                       ],
            [ 1, 'drag' , 'cmd.drag("'+sele+'")' ],
            [ 1, 'reset'       , 'cmd.reset(object="'+sele+'")'          ],
            [ 0, ''             , ''                       ],
            [ 1, 'rename'       , 'cmd.wizard("renaming","'+sele+'")'          ],
            [ 1, 'group' , lambda: move_to_group(self_cmd, sele) ],
            [ 0, ''             , ''                       ],
            [ 1, del_col + 'delete', 'cmd.delete("'+sele+'")'    ],
            ]

def mesh_action(self_cmd, sele):
    return [[ 2, 'Action:'     , ''                       ],
            [ 1, 'level'         , level(self_cmd, sele)  ],
            [ 0, ''             , ''                       ],
            [ 1, 'zoom'         , 'cmd.zoom("'+sele+'",animate=-1)'      ],
            [ 1, 'center'       , 'cmd.center("'+sele+'",animate=-1)'    ],
            [ 1, 'origin'       , 'cmd.origin("'+sele+'")'    ],
            [ 0, ''             , ''                       ],
            [ 1, 'drag' , 'cmd.drag("'+sele+'")' ],
            [ 1, 'reset'       , 'cmd.reset(object="'+sele+'")'          ],
            [ 0, ''             , ''                       ],
            [ 1, 'rename'       , 'cmd.wizard("renaming","'+sele+'")'          ],
            [ 1, 'group' , lambda: move_to_group(self_cmd, sele) ],
            [ 0, ''             , ''                       ],
            [ 1, del_col + 'delete', 'cmd.delete("'+sele+'")'    ],
            ]

def ramp_action(self_cmd, sele):
    return [[ 2, 'Action:'     , ''                       ],
              [ 1, 'levels', [
                  [ 1, 'Range +/- %.1f' % (L),
                      'cmd.ramp_update("%s", range=[%f, %f])' % (sele, -L, L) ]
                  for L in [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100]
                  ]],
              [ 0, '', '' ],
              [ 1, 'group' , lambda: move_to_group(self_cmd, sele) ],
              [ 0, '', '' ],
              [ 1, del_col + 'delete', 'cmd.delete("'+sele+'")'    ],
              ]

def ramp_color(self_cmd, sele):
    from pymol.creating import ramp_spectrum_dict
    return [[ 2, 'Color:', '' ],
            ] + [
                    [ 1, name, 'cmd.ramp_update("%s", color="%s")' % (sele, name) ]
                    for name in [
                        '[red, white, blue]',
                        ] + list(ramp_spectrum_dict)
            ]

def test1(self_cmd, sele):
        return [[ 2, 'Test1:'     , ''                      ],
              [ 1, 'zoom'         , 'cmd.zoom("all",animate=-1)'     ],
              [ 1, 'center'   , 'cmd.center("all",animate=-1)'   ],
              [ 1, 'origin'   , 'cmd.origin("all")'   ],
              ]

def test2(self_cmd, sele):
        return [[ 2, 'Test2:'     , ''                      ],
              [ 1, 'zoom'         , 'cmd.zoom("all",animate=-1)'     ],
              [ 1, 'center'   , 'cmd.center("all",animate=-1)'   ],
              [ 1, 'origin'   , 'cmd.origin("all")'   ],
              ]

def all_action(self_cmd, sele):
    return [[ 2, 'Action:'     , ''                      ],
              [ 1, 'zoom'         , 'cmd.zoom("all",animate=-1)'     ],
              [ 1, 'center'   , 'cmd.center("all",animate=-1)'   ],
              [ 1, 'origin'   , 'cmd.origin("all")'   ],
              [ 0, ''             , ''                      ],
              [ 1, 'preset'  , presets(self_cmd, "all")     ],
              [ 1, 'find', find(self_cmd, "all") ],
              [ 0, ''          ,''                                              ],
              [ 1, 'hydrogens' ,hydrogens(self_cmd, sele)     ],
#              [ 1, 'add hydrogens' ,'cmd.h_add("'+sele+'")'     ],
#              [ 1, 'remove hydrogens'  ,'cmd.remove("(hydro and ('+sele+'))")'     ],
              [ 1, rem_col + 'remove waters', 'cmd.remove("(solvent and ('+sele+'))")'     ],
              [ 0, ''             , ''                      ],
              [ 1, del_col + 'delete selections', 'map(cmd.delete,cmd.get_names("selections"))'     ],
              [ 0, ''          ,''                                              ],
              [ 1, del_col + 'delete everything', 'cmd.delete("all")'     ],
              [ 0, ''          ,''                                              ],
              [ 1, 'masking'      , masking(self_cmd, sele)         ],
              [ 1, 'movement'       , movement(self_cmd, sele)         ],
              [ 1, 'compute'        , compute(self_cmd, sele)         ],
              ]

def label_props(self_cmd, sele):
    return [[ 2, 'Other Properties:'       ,''                        ],

              [ 1, 'formal charge' ,
  'cmd.label("'+sele+'","(\'%+d\'%formal_charge) if formal_charge else \'\'")' ],
              [ 0, ''               , ''                                  ],
              [ 1, 'partial charge (0.00)' ,
  'cmd.label("'+sele+'","\'%.2f\'%partial_charge")'                      ],
              [ 1, 'partial charge (0.0000)' ,
  'cmd.label("'+sele+'","\'%.4f\'%partial_charge")'                      ],
              [ 0, ''               , ''                                  ],
              [ 1, 'elec. radius'       , 'cmd.label("'+sele+'","\'%1.2f\'%elec_radius")'  ],
              [ 0, ''               , ''                                  ],
              [ 1, 'text type'      , 'cmd.label("'+sele+'","text_type")'    ],
              [ 1, 'numeric type'   , 'cmd.label("'+sele+'","numeric_type")' ],
              [ 0, ''               , ''                                  ],
              [ 1, 'stereochemistry', 'cmd.label("'+sele+'","stereo")' ]
              ]

def label_ids(self_cmd, sele):
    return [[ 2, 'Atom Identifiers:'       ,''                        ],
              [ 1, 'rank'           , 'cmd.label("'+sele+'","rank")' ],
              [ 1, 'ID'             , 'cmd.label("'+sele+'","ID")' ],
              [ 1, 'index'          , 'cmd.label("'+sele+'","index")' ],
              ]

def mol_labels(self_cmd, sele):
    with menucontext(self_cmd, sele) as mc:
        return [[ 2, 'Label:'        , ''                                  ],
              [ 1, 'clear'          , 'cmd.label("'+sele+'","\'\'")'         ],
              [ 0, ''               , ''                                  ],
              [ 1, 'residues'       , """cmd.label('''(name """+self_cmd.get("label_anchor")+"""+C1*+C1' and (byres("""+sele+""")))''','''"%s-%s"%(resn,resi)''')"""  ],
              [ 1, 'residues (oneletter)', "cmd.label('''byca(" + sele + ")''', 'oneletter+resi')"],
              [ 1, 'chains'       ,   'util.label_chains("'+sele+'",_self=cmd)'  ],
              [ 1, 'segments'       ,   'util.label_segments("'+sele+'",_self=cmd)'  ],
              [ 0, ''               , ''                                  ],
              [ 1, 'atom name'      , 'cmd.label("'+sele+'","name")'         ],
              [ 1, 'element symbol' , 'cmd.label("'+sele+'","elem")'         ],
              [ 1, 'residue name'  , 'cmd.label("'+sele+'","resn")'         ],
              [ 1, 'one letter code', 'cmd.label("'+sele+'","oneletter")' ],
              [ 1, 'residue identifier'    , 'cmd.label("'+sele+'","resi")'         ],
              [ 1, 'chain identifier' , 'cmd.label("'+sele+'","chain")'         ],
              [ 1, 'segment identifier'       , 'cmd.label("'+sele+'","segi")'         ],
              [ 0, ''               , ''                                  ],
              [ 1, 'b-factor'       , 'cmd.label("'+sele+'","\'%1.2f\'%b")'  ],
              [ 1, 'occupancy'       , 'cmd.label("'+sele+'","\'%1.2f\'%q")'  ],
              [ 1, 'vdw radius'       , 'cmd.label("'+sele+'","\'%1.2f\'%vdw")'  ],
              [ 0, ''               , ''                                  ],
              [ 1, 'other properties' , label_props(self_cmd, sele) ],
              [ 0, ''               , ''                                  ],
              [ 1, 'atom identifiers' , label_ids(self_cmd, sele) ],
              ]

def mol_ss(self_cmd, sele):
    fmt = '''cmd.alter("{}","ss='{}'",quiet=0);cmd.rebuild()'''
    return [
        [2, 'Secondary Structure:', '' ],
        [1, 'H (helix)', fmt.format(sele, 'H')],
        [1, 'S (sheet)', fmt.format(sele, 'S')],
        [1, 'L (loop)' , fmt.format(sele, 'L')],
        [0, '', ''],
        [1, 'dss' , 'cmd.dss("{0}", context="byobject ({0})")'.format(sele) ],
    ]

def mol_view(self_cmd, sele):
    return [
        [ 1, 'zoom'           ,'cmd.zoom("'+sele+'",animate=-1)'            ],
        [ 1, 'orient'           ,'cmd.orient("'+sele+'",animate=-1)'            ],
        [ 1, 'center'           ,'cmd.center("'+sele+'",animate=-1)'            ],
        [ 1, 'origin'           ,'cmd.origin("'+sele+'")'            ],
        ]

def all_option(self_cmd, sele):
    return [
        [ 2, '(all)'      , '' ],
        [ 1, 'show'      , mol_show(self_cmd, sele) ],
        [ 1, 'hide'      , mol_hide(self_cmd, sele) ],
        [ 1, 'color'      , mol_color(self_cmd, sele) ],
#      [ 1, 'view'      , mol_view(self_cmd, sele) ],
        [ 1, 'preset'      , presets(self_cmd, sele) ],
        [ 0, ''             , ''                      ],
        [ 1, 'zoom'           ,'cmd.zoom("'+sele+'",animate=-1)'            ],
        [ 1, 'orient'           ,'cmd.orient("'+sele+'",animate=-1)'            ],
        [ 1, 'center'           ,'cmd.center("'+sele+'",animate=-1)'            ],
        [ 1, 'origin'           ,'cmd.origin("'+sele+'")'            ],
        [ 1, 'select'        ,'cmd.select("'+sele+'",enable=1,merge=2)'            ],
        [ 0, ''             , ''                      ],
        [ 1, 'label'      , mol_labels(self_cmd, sele) ],
        [ 0, '', '' ],
        [ 1, 'enable'         ,'cmd.enable("'+sele+'")'            ],
        [ 1, 'disable'        ,'cmd.disable("'+sele+'")'            ],
        ]

def enable_disable(self_cmd, enable):
    names_enabled = self_cmd.get_names('objects', enabled_only=1)
    if enable:
        result = [[ 2, 'Enable', '' ]]
        cmmd = 'enable '
        names = [ob for ob in self_cmd.get_names('objects')
                 if ob not in names_enabled]
    else:
        result = [[ 2, 'Disable', '']]
        cmmd = 'disable '
        names = names_enabled
    names = ['all'] + names
    result += [[1, ob, cmmd + ob] for ob in names]
    result.insert(2, [0, '', ''])
    if not enable:
        result.insert(2, [1, 'selections', 'deselect'])
    return result

def scene_buttons(self_cmd):
    return [[ 2, 'Buttons', '' ],
            [ 1, 'on', 'cmd.set("scene_buttons")'],
            [ 1, 'off', 'cmd.set("scene_buttons", 0)']]

def scene_main(self_cmd):
    list = self_cmd.get_scene_list()
    recall_list = [ [2, 'Scenes' , ''] ]
    for entry in list:
        recall_list.append([1,entry,'cmd.scene("""'+entry+'""")'])
    return [
        [ 2, 'Scene', '' ],
        [ 1, 'next' , 'cmd.scene()' ],
        [ 0, ''             , ''                      ],
        [ 1, 'append' , 'cmd.scene("new","append",quiet=0)' ],
        [ 1, 'update' , 'cmd.scene("auto","update",quiet=0)' ],
        [ 0, ''             , ''                      ],
        [ 1, 'recall' , recall_list ],
        [ 0, ''             , ''                      ],
        [ 1, 'buttons', scene_buttons(self_cmd)] ]

def main_pseudoatom_sub(self_cmd,pos, screenpos):
    return [
        [ 2, 'Pseudoatom' ,'' ],
        [ 1, 'label' ,'cmd.wizard("pseudoatom","label",pos=[%1.7f,%1.7f,%1.7f])'%pos ],
        [ 0, ''             , ''                      ],
        [ 1, 'single ', 'cmd.pseudoatom(pos=[%1.7f,%1.7f,%1.7f])'%pos ],
        ]

def main_pseudoatom(self_cmd,pos, screenpos):
    return [
        [ 2, 'New'  , '' ],
        [ 1, 'pseudoatom' , main_pseudoatom_sub(self_cmd,pos, screenpos) ],
        ]

def movie_panel(self_cmd):
    return [[ 2, 'Panel', '' ],
            [ 1, 'on', 'cmd.set("movie_panel")'],
            [ 1, 'off', 'cmd.set("movie_panel", 0)']]

def movie_main(self_cmd):
    return [
        [ 2, 'Movie', ''],
        [ 1, 'play', 'cmd.mplay()'],
        [ 1, 'stop', 'cmd.mstop()'],
        [ 0, '', '' ],
        [ 1, 'rewind', 'cmd.rewind()'],
        [ 0, '', '' ],
        [ 1, 'panel', movie_panel(self_cmd) ]
        ]

def main_menu(self_cmd,pos, screenpos):
    return [
        [ 2, 'Main Pop-Up'  , '' ],
        [ 1, 'new'             , main_pseudoatom(self_cmd,pos, screenpos) ],
        [ 0, ''             , ''                      ],
        [ 1, 'zoom (vis)'           ,'cmd.zoom("visible",animate=-1)'            ],
        [ 1, 'orient (vis)'           ,'cmd.orient("visible",animate=-1)'            ],
        [ 1, 'center (vis)'           ,'cmd.center("visible",animate=-1)'            ],
        [ 1, 'reset'           ,'cmd.reset()'            ],
        [ 0, ''             , ''                      ],
        [ 1, 'movie'           , movie_main(self_cmd) ],
        [ 1, 'scene'           , scene_main(self_cmd) ],
        [ 0, ''             , ''                      ],
        [ 1, 'enable', enable_disable(self_cmd, 1) ],
        [ 1, 'disable', enable_disable(self_cmd,0) ],
        [ 0, ''             , ''                      ],
        [ 1, '(all)'      , all_option(self_cmd,"all") ],
        [ 1, '(visible)'      , all_option(self_cmd,"visible") ],
        [ 0, ''             , ''                      ],
        [ 1, 'ray'           ,'cmd.ray()' ],
        [ 0, ''             , ''                      ],
        [ 1, del_col + 'delete all', 'cmd.delete("all")' ],
        [ 1, del_col + 'reinitialize', 'cmd.reinitialize()' ],
        [ 1, del_col + 'quit', 'cmd.quit()' ],
        ]


def pick_sele(self_cmd, sele, title):
    result = [
        [ 2, title, '' ],
        [ 1, 'disable'    , 'cmd.disable("'+sele+'")' ],
        [ 0, ''             , ''                      ],
        [ 1, 'actions', sele_action2(self_cmd, sele) ],
        [ 0, ''             , ''                      ],
        [ 1, 'color'      , mol_color(self_cmd, sele) ],
        [ 1, 'show'      , mol_show(self_cmd, sele) ],
        [ 1, 'hide'      , mol_hide(self_cmd, sele) ],
        [ 1, 'preset'  , presets(self_cmd, sele)       ],
        [ 1, 'label'       , mol_labels(self_cmd, sele) ],
        [ 1, 'ss'        , mol_ss(self_cmd, sele) ],
        [ 0, ''             , ''                      ],
        [ 1, 'zoom'           ,'cmd.zoom("'+sele+'",animate=-1)'            ],
        [ 1, 'orient'           ,'cmd.orient("'+sele+'",animate=-1)'            ],
        [ 1, 'center'           ,'cmd.center("'+sele+'",animate=-1)'            ],
        [ 1, 'origin'           ,'cmd.origin("'+sele+'")'            ],
        [ 0, ''               ,''                             ],
        [ 1, 'drag'             ,'cmd.drag("'+sele+'")'            ],
        [ 1, 'clean'             ,'cmd.clean("'+sele+'")'            ],
        [ 0, ''               ,''                             ],
        [ 1, rem_col + 'remove', 'cmd.remove("'+sele+'")'            ],
        ]
    return result

def pick_option(self_cmd, sele, title, object=0):
    if object:
        sele = self_cmd.identify(sele, 1)[0][0]
    result = [
        [ 2, title, '' ],
        [ 1, 'color'      , lambda: mol_color(self_cmd, sele) ],
        [ 1, 'show'      , mol_show(self_cmd, sele) ],
        [ 1, 'hide'      , mol_hide(self_cmd, sele) ],
        [ 1, 'preset'  , presets(self_cmd, sele)       ],
        [ 1, 'label'          , mol_labels(self_cmd, sele) ],
        [ 0, ''             , ''                      ],
        [ 1, 'zoom'           ,'cmd.zoom("'+sele+'",animate=-1)'            ],
        [ 1, 'orient'           ,'cmd.orient("'+sele+'",animate=-1)'            ],
        [ 1, 'center'           ,'cmd.center("'+sele+'",animate=-1)'            ],
        [ 1, 'origin'           ,'cmd.origin("'+sele+'")'            ],
        [ 1, 'select'        ,'cmd.select("'+sele+'",enable=1,merge=2)'            ],
        [ 0, ''               ,''                             ]]

    if object:
        result.append([ 1, 'drag'             ,   [[ 1, 'coordinates', 'cmd.drag("'+sele+'")'],
                                     [ 1, 'matrix', 'cmd.drag("'+sele+'",mode=1)']]])
    else:
        result.append([ 1, 'drag'   ,  'cmd.drag("'+sele+'")'])

    result.extend([
        [ 1, 'clean'             ,'cmd.clean("'+sele+'")'            ],
        [ 1, 'masking'        , masking(self_cmd, sele)         ],
        [ 1, 'movement'       , movement(self_cmd, sele)         ],
        ])

    if object:
        result.extend([
            [ 1, del_col + 'delete', 'cmd.delete("'+sele+'")'            ],
            [ 0, ''             , ''                      ],
            [ 1, 'disable'        ,'cmd.disable("'+sele+'")'            ],
            [ 1, 'disable others' ,'cmd.disable("*");cmd.enable("'+sele+'", 1)' ],
            ])
    else:
        result.extend([
            [ 1, rem_col + 'remove atoms', 'cmd.remove("'+sele+'")' ],
            [ 0, ''             , ''                      ],
            [ 1, 'copy to object' , lambda: copy_to(self_cmd, sele) ],
            [ 1, 'extract object' ,'cmd.extract(None,"'+sele+'",zoom=0)' ],
            ])
    return result


def pick_menu(self_cmd, title, sele2):
    with menucontext(self_cmd, sele2):
        return [[ 2, title     , '' ],
            [ 1, 'drag object matrix'      ,'cmd.drag("(byobj ('+sele2+'))",mode=1)'            ],
            [ 1, 'drag object coords'      ,'cmd.drag("(byobj ('+sele2+'))")'            ],
            [ 0, ''             , ''                      ],
            [ 1, 'atom'    , pick_option(self_cmd, sele2, "Atom") ],
            [ 1, 'residue' , pick_option(self_cmd, "(byres ("+sele2+"))", "Residue") ],
            [ 1, 'chain'   , pick_option(self_cmd, "(bychain ("+sele2+"))", "Chain") ],
            [ 1, 'segment' , pick_option(self_cmd, "(byseg ("+sele2+"))", "Segment") ],
            [ 1, 'object'  , pick_option(self_cmd, sele2, "Object",1) ],
            [ 0, ''             , ''                      ],
            [ 1, 'molecule', pick_option(self_cmd, "(bymol ("+sele2+"))", "Molecule") ],
            [ 0, ''             , ''                      ],
            [ 1, 'fragment', pick_option(self_cmd, "(byfrag ("+sele2+"))", "Fragment") ],
            [ 1, 'fragment+joint(s)', pick_option(self_cmd, "((byfrag ("+sele2+")) extend 1)", "Fragment") ],
            ]

def seq_option(self_cmd, sele, title, object=0):
    c=len(title)-1
    while title[c]!='/':
        c = c-1
    title = title[0:c+1]

    result = [
        [ 2, title, '' ],
        [ 1, 'color'     , mol_color(self_cmd, sele) ],
        [ 1, 'show'      , mol_show(self_cmd, sele) ],
        [ 1, 'hide'      , mol_hide(self_cmd, sele) ],
        [ 1, 'preset'    , presets(self_cmd, sele)       ],
        [ 1, 'label'     , mol_labels(self_cmd, sele) ],
        [ 1, 'ss'        , mol_ss(self_cmd, sele) ],
        [ 0, ''          , ''                      ],
        [ 1, 'zoom'      ,'cmd.zoom("'+sele+'",animate=-1)'            ],
        [ 1, 'orient'    ,'cmd.orient("'+sele+'",animate=-1)'            ],
        [ 1, 'center'    ,'cmd.center("'+sele+'",animate=-1)'            ],
        [ 1, 'origin'    ,'cmd.origin("'+sele+'")'            ],
        [ 1, 'select'    ,'cmd.select("'+sele+'",enable=1,merge=2)'            ],
        [ 0, ''               ,''                             ],
        [ 1, 'drag'      ,'cmd.drag("'+sele+'")'            ],
        [ 1, 'clean'      ,'cmd.clean("'+sele+'")'            ],
        ]

    if object:
        result.extend([
            [ 0, ''             , ''                      ],
            [ 1, 'disable'        ,'cmd.disable("'+sele+'")'            ],
            [ 0, ''             , ''                      ],
            [ 1, del_col + 'delete', 'cmd.delete("'+sele+'")'            ]
            ])
    else:
        result.extend([
            [ 0, ''             , ''                      ],
            [ 1, 'create object','cmd.create(None,"'+sele+'",zoom=0)'            ],
            [ 1, 'extract object' ,'cmd.extract(None,"'+sele+'",zoom=0)' ],
            [ 0, ''             , ''                      ],
            [ 1, rem_col + 'remove atoms', 'cmd.remove("'+sele+'")' ],
            ])
    return result

def scene_menu(self_cmd, name):
    safe_name = name.replace('"','\\"') # just in case
    return [[ 2, 'Scene '+name    , '' ],
            [ 1, 'rename', 'cmd.wizard("renaming","'+name+'",mode="scene")'          ],
            [ 0, ''             , ''                      ],
            [ 1, 'update', 'cmd.scene("'+safe_name+'","update")'],
            [ 0, ''             , ''                      ],
            [ 1, del_col + 'delete', 'cmd.scene("'+safe_name+'","delete")'],
            ]
