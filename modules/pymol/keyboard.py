import editor
import cmd

def get_special(self_cmd=cmd):
    return {
        1        : [ 'F1'        , None                   , () , {} ],
        2        : [ 'F2'        , None                   , () , {} ],
        3        : [ 'F3'        , None                   , () , {} ],
        4        : [ 'F4'        , None                   , () , {} ],
        5        : [ 'F5'        , None                   , () , {} ],
        6        : [ 'F6'        , None                   , () , {} ],
        7        : [ 'F7'        , None                   , () , {} ],
        8        : [ 'F8'        , None                   , () , {} ],
        9        : [ 'F9'        , None                   , () , {} ],
        10       : [ 'F10'       , None                   , () , {} ],
        11       : [ 'F11'       , None                   , () , {} ],
        12       : [ 'F12'       , None                   , () , {} ],
        100      : [ 'left'      , self_cmd.backward               , () , {} ],
        101      : [ 'up'        , None                   , () , {} ],
        102      : [ 'right'     , self_cmd.forward                , () , {} ],
        103      : [ 'down'      , None                   , () , {} ],
        104      : [ 'pgup'      , self_cmd.scene                  , ('','previous') , {} ],
        105      : [ 'pgdn'      , self_cmd.scene                  , ('','next') , {} ],
        106      : [ 'home'      , self_cmd.zoom                   , () ,  {'animate':-1} ],
        107      : [ 'end'       , self_cmd.mtoggle                , () , {} ],
        108      : [ 'insert'    , self_cmd.rock                   , () , {} ]   
        }

def get_shft_special(self_cmd=cmd):
    return {
        1        : [ 'F1'        , None                   , () , {} ],
        2        : [ 'F2'        , None                   , () , {} ],
        3        : [ 'F3'        , None                   , () , {} ],
        4        : [ 'F4'        , None                   , () , {} ],
        5        : [ 'F5'        , None                   , () , {} ],
        6        : [ 'F6'        , None                   , () , {} ],
        7        : [ 'F7'        , None                   , () , {} ],
        8        : [ 'F8'        , None                   , () , {} ],
        9        : [ 'F9'        , None                   , () , {} ],
        10       : [ 'F10'       , None                   , () , {} ],
        11       : [ 'F11'       , None                   , () , {} ],
        12       : [ 'F12'       , None                   , () , {} ],
        100      : [ 'left'      , self_cmd.backward               , () , {} ],
        101      : [ 'up'        , None                   , () , {} ],
        102      : [ 'right'     , self_cmd.forward                , () , {} ],
        103      : [ 'down'      , None                   , () , {} ],
        104      : [ 'pgup'      , self_cmd.scene                  , ('','previous') , {} ],
        105      : [ 'pgdn'      , self_cmd.scene                  , ('','next') , {} ],
        106      : [ 'home'      , self_cmd.rewind                 , () ,  {'animate':-1} ],
        107      : [ 'end'       , self_cmd.ending                 , () , {} ],
        108      : [ 'insert'    , self_cmd.rock                   , () , {} ]   
        }

def get_alt_special(self_cmd=cmd):
    return { # NOTE: some OSes/Windowing systems intercept ALT-Fn keys.
            1        : [ 'F1'        , None                   , () , {} ],
            2        : [ 'F2'        , None                   , () , {} ],
            3        : [ 'F3'        , None                   , () , {} ],
            4        : [ 'F4'        , None                   , () , {} ],
            5        : [ 'F5'        , None                   , () , {} ],
            6        : [ 'F6'        , None                   , () , {} ],
            7        : [ 'F7'        , None                   , () , {} ],
            8        : [ 'F8'        , None                   , () , {} ],
            9        : [ 'F9'        , None                   , () , {} ],
            10       : [ 'F10'       , None                   , () , {} ],
            11       : [ 'F11'       , None                   , () , {} ],
            12       : [ 'F12'       , None                   , () , {} ],
            100      : [ 'left'      , self_cmd.backward               , () , {} ],
            101      : [ 'up'        , None                   , () , {} ],
            102      : [ 'right'     , self_cmd.forward                , () , {} ],
            103      : [ 'down'      , None                   , () , {} ],
            104      : [ 'pgup'      , self_cmd.rewind                 , () , {} ],
            105      : [ 'pgdn'      , self_cmd.ending                 , () , {} ],
            106      : [ 'home'      , self_cmd.zoom                   , () ,  {'animate':-1} ],
            107      : [ 'end'       , self_cmd.ending                 , () , {} ],
            108      : [ 'insert'    , self_cmd.rock                   , () , {} ]   
        }

def get_ctrl_special(self_cmd=cmd):
    return { # NOTE: some OSes/Windowing systems intercept CTRL-Fn keys.
            1        : [ 'F1'        , self_cmd.scene  , ('F1','store') , {} ],
            2        : [ 'F2'        , self_cmd.scene,('F2','store')    , {} ],
            3        : [ 'F3'        , self_cmd.scene,('F3','store')    , {} ],
            4        : [ 'F4'        , self_cmd.scene,('F4','store')    , {} ],
            5        : [ 'F5'        , self_cmd.scene,('F5','store')    , {} ],
            6        : [ 'F6'        , self_cmd.scene,('F6','store')    , {} ],
            7        : [ 'F7'        , self_cmd.scene,('F7','store')    , {} ],
            8        : [ 'F8'        , self_cmd.scene,('F8','store')    , {} ],
            9        : [ 'F9'        , self_cmd.scene,('F9','store')    , {} ],
            10       : [ 'F10'       , self_cmd.scene,('F10','store')   , {} ],
            11       : [ 'F11'       , self_cmd.scene,('F11','store')   , {} ],
            12       : [ 'F12'       , self_cmd.scene,('F12','store')   , {} ],
            100      : [ 'left'      , self_cmd.backward               , () , {} ],
            101      : [ 'up'        , None                   , () , {} ],
            102      : [ 'right'     , self_cmd.forward                , () , {} ],
            103      : [ 'down'      , None                   , () , {} ],
            104      : [ 'pgup'      , self_cmd.scene                  , ('','insert_before') , {} ],
            105      : [ 'pgdn'      , self_cmd.scene                  , ('','insert_after') , {} ],
            106      : [ 'home'      , self_cmd.zoom                   , () , {'animate':-1} ],
            107      : [ 'end'       , self_cmd.scene                  , ('new','store') , {} ],
            108      : [ 'insert'    , self_cmd.scene                  , ('auto','store') , {} ]   
        }

def get_ctsh_special(self_cmd=cmd):
    return { # NOTE: some OSes/Windowing systems intercept CTRL-Fn keys.
        1        : [ 'F1'        , self_cmd.scene,('SHFT-F1','store') , {} ],
        2        : [ 'F2'        , self_cmd.scene,('SHFT-F2','store')    , {} ],
        3        : [ 'F3'        , self_cmd.scene,('SHFT-F3','store')    , {} ],
        4        : [ 'F4'        , self_cmd.scene,('SHFT-F4','store')    , {} ],
        5        : [ 'F5'        , self_cmd.scene,('SHFT-F5','store')    , {} ],
        6        : [ 'F6'        , self_cmd.scene,('SHFT-F6','store')    , {} ],
        7        : [ 'F7'        , self_cmd.scene,('SHFT-F7','store')    , {} ],
        8        : [ 'F8'        , self_cmd.scene,('SHFT-F8','store')    , {} ],
        9        : [ 'F9'        , self_cmd.scene,('SHFT-F9','store')    , {} ],
        10       : [ 'F10'       , self_cmd.scene,('SHFT-F10','store')   , {} ],
        11       : [ 'F11'       , self_cmd.scene,('SHFT-F11','store')   , {} ],
        12       : [ 'F12'       , self_cmd.scene,('SHFT-F12','store')   , {} ],
        100      : [ 'left'      , self_cmd.backward               , () , {} ],
        101      : [ 'up'        , None                   , () , {} ],
        102      : [ 'right'     , self_cmd.forward                , () , {} ],
        103      : [ 'down'      , None                   , () , {} ],
        104      : [ 'pgup'      , self_cmd.scene                  , ('','insert_before') , {} ],
        105      : [ 'pgdn'      , self_cmd.ending                 , ('','insert_after') , {} ],
        106      : [ 'home'      , self_cmd.zoom                 , () ,  {'animate':-1} ],
        107      : [ 'end'       , self_cmd.mtoggle                , () , {} ],
        108      : [ 'insert'    , self_cmd.rock                   , () , {} ]   
        }

def get_ctrl(self_cmd=cmd):
    return {
        'A' : [ self_cmd.redo                   , () , {}],
        'B' : [ self_cmd.replace                , ('Br',1,1), {} ],
        'C' : [ self_cmd.replace                , ('C',4,4), {} ],
        'D' : [ self_cmd.remove_picked          , () , {'quiet':0} ],
        'E' : [ self_cmd.invert                 , () , {'quiet':0} ],      
        'F' : [ self_cmd.replace                , ('F',1,1), {} ],   
        'G' : [ self_cmd.replace                , ('H',1,1), {} ],
        'I' : [ self_cmd.replace                , ('I',1,1), {} ],
        'J' : [ self_cmd.alter                  , ('pk1','formal_charge=-1.0'), {} ],
        'K' : [ self_cmd.alter                  , ('pk1','formal_charge =1.0'), {} ],
        'L' : [ self_cmd.replace                , ('Cl',1,1) , {}],
        'N' : [ self_cmd.replace                , ('N',4,3) , {}],
        'O' : [ self_cmd.replace                , ('O',4,2) , {}],   
        'P' : [ self_cmd.replace                , ('P',4,1) , {}],
        'Q' : [ None                   , () , {}],   
        'R' : [ self_cmd.h_fill                 , () , {} ],   
        'S' : [ self_cmd.replace                , ('S',4,2) , {}],
        'T' : [ lambda a,b,c=self_cmd.bond,u=self_cmd.unpick: (c(a,b),u()) , ('pk1','pk2') , {} ],   
        'U' : [ self_cmd.alter                  , ('pk1','formal_charge =0.0') , {}],
        'W' : [ self_cmd.cycle_valence          , () , {}],
        #         'X' : [ lambda a,b,c,d=dist,u=unpick:(d(a,b,c),u()), (None,'pk1','pk2') , {} ],
        'X' : [ self_cmd.auto_measure           , () , {} ],   
        'Y' : [ self_cmd.attach                 , ('H',1,1) , {} ],
        'Z' : [ self_cmd.undo                   , () , {} ],   
        }

def get_alt(self_cmd=cmd):
    return {
        '1' : [ editor.attach_fragment  , ("pk1","formamide",5,0), {}],
        '2' : [ editor.attach_fragment  , ("pk1","formamide",4,0), {}],
        '3' : [ editor.attach_fragment  , ("pk1","sulfone",3,1), {}],
        '4' : [ editor.attach_fragment  , ("pk1","cyclobutane",4,0), {}],
        '5' : [ editor.attach_fragment  , ("pk1","cyclopentane",5,0), {}],
        '6' : [ editor.attach_fragment  , ("pk1","cyclohexane",7,0), {}],
        '7' : [ editor.attach_fragment  , ("pk1","cycloheptane",8,0), {}],
        '8' : [ editor.attach_fragment  , ("pk1","cyclopentadiene",5,0), {}],
        '9' : [ editor.attach_fragment  , ("pk1","benzene",6,0), {}],
        '0' : [ editor.attach_fragment  , ("pk1","formaldehyde",2,0), {}],
        'a' : [ editor.attach_amino_acid, ("pk1","ala"), {}],
        'b' : [ editor.attach_amino_acid, ("pk1","ace"), {}],                                 
        'c' : [ editor.attach_amino_acid, ("pk1","cys"), {}],
        'd' : [ editor.attach_amino_acid, ("pk1","asp"), {}],
        'e' : [ editor.attach_amino_acid, ("pk1","glu"), {}],
        'f' : [ editor.attach_amino_acid, ("pk1","phe"), {}],
        
        'g' : [ editor.attach_amino_acid, ("pk1","gly"), {}],
        'h' : [ editor.attach_amino_acid, ("pk1","his"), {}],
        'i' : [ editor.attach_amino_acid, ("pk1","ile"), {}],
        
        'j' : [ editor.attach_fragment,   ("pk1","acetylene",2,0), {}],
        'k' : [ editor.attach_amino_acid, ("pk1","lys"), {}],
        'l' : [ editor.attach_amino_acid, ("pk1","leu"), {}],
        
        'm' : [ editor.attach_amino_acid, ("pk1","met"), {}],
        'n' : [ editor.attach_amino_acid, ("pk1","asn"), {}],
        'p' : [ editor.attach_amino_acid, ("pk1","pro"), {}],
        'q' : [ editor.attach_amino_acid, ("pk1","gln"), {}],
        'r' : [ editor.attach_amino_acid, ("pk1","arg"), {}],
        
        's' : [ editor.attach_amino_acid, ("pk1","ser"), {}],
        't' : [ editor.attach_amino_acid, ("pk1","thr"), {}],
        'v' : [ editor.attach_amino_acid, ("pk1","val"), {}],
        'w' : [ editor.attach_amino_acid, ("pk1","trp"), {}],
        'y' : [ editor.attach_amino_acid, ("pk1","tyr"), {}],
        'z' : [ editor.attach_amino_acid, ("pk1","nme"), {}],
        }

