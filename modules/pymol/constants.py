
# constant objects

from .parsing import QuietException
from .shortcut import Shortcut
from .constants_palette import palette_dict
import re

class _loadable:
    pdb = 0
    mol = 1
    molstr = 3
    mmod = 4
    mmodstr = 6
    xplor = 7
    model = 8
    pdbstr = 9
    brick = 10    # chempy.brick object
    chempymap = 11 # chempy.map object (only used by cctbx/examples/view_fft_map.py)
    callback = 12 # pymol callback obejct
    cgo = 13      # compiled graphic object
    xyz = 15      # xyz, tinker format
    ccp4 = 18     # CCP4 map, under development
    pmo = 19      # pmo, experimental molecular object format
    top = 21      # AMBER topology
    trj = 22      # AMBER trajectory
    crd = 23      # AMBER coordinate
    rst = 24      # AMBER restart
    xplorstr = 26 # XPLOR map as string
    phi = 27      # Delphi/Grasp
    fld = 28      # AVS field format (not yet general -- only uniform allowed)
    brix = 29     # BRIX/DSN6/O map format
    grd = 30      # Insight II Grid format
    pqr = 31      # PQR file (modified PDB file for APBS)
    dx = 32       # DX file (APBS)
    mol2 = 33     # MOL2 file (TRIPOS)
    mol2str = 34  # MOL2 file string (TRIPOS)
    p1m = 35      # P1M file (combined data & secure commands)
    ccp4str = 36  # CCP4 map string
    sdf = 37      # new default...
    sdf2 = 37     # SDF using C-based SDF parser (instead of Python)
    sdf2str = 38  # SDF ditto
    trj2 = 45     # trj trajectroy format (via plugin)
    xyzstr = 49   #
    phistr = 51   # electrostatic map as a string
    acnt = 56     # Tripos/Sybyl acnt grid file (proposed)
    cif = 60      # C++ based CIF parser
    cifstr = 61
    plugin = 64
    mae = 65
    maestr = 66
    pdbqt = 67
    vdbstr = 69
    vdb = 70
    mmtf = 71
    mmtfstr = 72
    map = 73      # unspecified CCP4 or MRC map
    mrc = 74
    dxstr = 75    # DX file (APBS)
    mapstr = 76   # unspecified CCP4 or MRC map
    mrcstr = 77

class loadable(_loadable):
    @classmethod
    def _reverse_lookup(cls, number):
        for name in dir(cls):
            if getattr(cls, name) == number:
                return name
        return ''

_load2str = { loadable.pdb : loadable.pdbstr,
              loadable.vdb: loadable.vdbstr,
              loadable.cif : loadable.cifstr,
              loadable.mmtf : loadable.mmtfstr,
              loadable.mae : loadable.maestr,
              loadable.mol : loadable.molstr,
              loadable.xplor : loadable.xplorstr,
              loadable.mol2 : loadable.mol2str,
              loadable.mmod : loadable.mmodstr,
              loadable.ccp4 : loadable.ccp4str,
              loadable.mrc : loadable.mrcstr,
              loadable.map : loadable.mapstr,
              loadable.dx : loadable.dxstr,
              loadable.xyz  : loadable.xyzstr,
              loadable.sdf2 : loadable.sdf2str}

sanitize_alpha_list_re = re.compile(r"[^a-zA-Z0-9_\'\"\.\-\[\]\,]+")
nt_hidden_path_re = re.compile(r"\$[\/\\]")

def safe_alpha_list_eval(st):
    '''Like `safe_eval` but removes most non-alpha-numeric characters.

    >>> safe_alpha_list_eval("[A B/C, D+E:F]")
    ['ABC', 'DEF']
    '''
    st = sanitize_alpha_list_re.sub('',st)
    return safe_list_eval(st)

class SafeEvalNS(object):
    def __getitem__(self, name):
        return name

def safe_eval(st):
    '''Safe version of "eval" which evaluates names to strings.

    # "foo" is a string
    >>> safe_eval('foo, 123, 4 + 5, "A B C", {}, "{}"')
    ('foo', 123, 9, 'A B C', {}, '{}')

    # no harmful code possible
    >>> safe_eval('__import__("os").unlink("foo.txt")')
    TypeError: 'str' object is not callable
    '''
    return eval(st, {}, SafeEvalNS())

safe_list_eval = safe_eval

DEFAULT_ERROR = -1
DEFAULT_SUCCESS = None

#--------------------------------------------------------------------
# shortcuts...

toggle_dict = {'on':1,'off':0,'1':1,'0':0,'toggle':-1, '-1':-1}
toggle_sc = Shortcut(toggle_dict.keys())

stereo_dict = {'on':-2,'off':0,'0':0,'1':-2,'swap':-1,
               'chromadepth': -3,
               'quadbuffer':1,'crosseye':2,
               'walleye':3,'geowall':4,'sidebyside':5,
               'byrow':6, 'bycolumn':7, 'checkerboard':8,
               'custom': 9, 'anaglyph' : 10,
               'dynamic' : 11, 'clonedynamic': 12,
               'openvr' : 13 }

stereo_sc = Shortcut(stereo_dict.keys())

space_dict = {
    'cmyk': "$PYMOL_DATA/pymol/cmyk.png",
    'pymol': 'pymol',
    'rgb': 'rgb',
    'greyscale': 'greyscale',
}

space_sc = Shortcut(space_dict)

window_dict = { 'show' : 1, 'hide' : 0, 'position' : 2, 'size' : 3,
                'box' : 4, 'maximize' : 5, 'fit' : 6, 'focus' : 7,
                'defocus' : 8 }
window_sc = Shortcut(window_dict.keys())

repres = {
    'everything'    : -1,
    'sticks'        : 0,
    'spheres'       : 1,
    'surface'       : 2,
    'labels'        : 3,
    'nb_spheres'    : 4,
    'cartoon'       : 5,
    'ribbon'        : 6,
    'lines'         : 7,
    'mesh'          : 8,
    'dots'          : 9,
    'dashes'        :10,
    'nonbonded'     :11,
    'cell'          :12,
    'cgo'           :13,
    'callback'      :14,
    'extent'        :15,
    'slice'         :16,
    'angles'        :17,
    'dihedrals'     :18,
    'ellipsoids'    :19,
    'volume'        :20,
}

repmasks = {
    'everything'    : 0b000111111111111111111111,
    'sticks'        : 0b000000000000000000000001,
    'spheres'       : 0b000000000000000000000010,
    'surface'       : 0b000000000000000000000100,
    'labels'        : 0b000000000000000000001000,
    'nb_spheres'    : 0b000000000000000000010000,
    'cartoon'       : 0b000000000000000000100000,
    'ribbon'        : 0b000000000000000001000000,
    'lines'         : 0b000000000000000010000000,
    'mesh'          : 0b000000000000000100000000,
    'dots'          : 0b000000000000001000000000,
    'dashes'        : 0b000000000000010000000000,
    'nonbonded'     : 0b000000000000100000000000,
    'cell'          : 0b000000000001000000000000,
    'cgo'           : 0b000000000010000000000000,
    'callback'      : 0b000000000100000000000000,
    'extent'        : 0b000000001000000000000000,
    'slice'         : 0b000000010000000000000000,
    'angles'        : 0b000000100000000000000000,
    'dihedrals'     : 0b000001000000000000000000,
    'ellipsoids'    : 0b000010000000000000000000,
    'volume'        : 0b000100000000000000000000,

    # combinations
    'licorice'      : 0b000000000000000000010001, # sticks | nb_spheres
    'wire'          : 0b000000000000100010000000, # lines | nonbonded
}

repmasks_sc = Shortcut(repmasks)
repres_sc = Shortcut(repres.keys())

boolean_dict = {
    'yes'           : 1,
    'no'            : 0,
    '1'             : 1,
    '0'             : 0,
    'on'            : 1,
    'off'           : 0
    }

boolean_sc = Shortcut(boolean_dict.keys())

palette_sc = Shortcut(palette_dict.keys())


location_code = {
    'first' : -1,
    'top' : -1,
    'upper' : -2,
    'current' : 0,
    'bottom' : 1,
    'last' : 1
    }
location_sc = Shortcut(location_code.keys())

class fb_action:
    set = 0
    enable = 1
    disable = 2
    push = 3
    pop = 4

class fb_module:

# This first set represents internal C systems

    all                       =0
    isomesh                   =1
    map                       =2
    matrix                    =3
    mypng                     =4
    triangle                  =5
    match                     =6
    raw                       =7
    isosurface                =8
    opengl                    =9

    color                     =10
    cgo                       =11
    feedback                  =12
    scene                     =13
    threads                   =14
    symmetry                  =15
    ray                       =16
    setting                   =17
    object                    =18
    ortho                     =19
    movie                     =20
    python                    =21
    extrude                   =22
    rep                       =23
    shaker                    =24

    coordset                  =25
    distset                   =26
    gadgetset                 =27

    objectmolecule            =30
    objectmap                 =31
    objectmesh                =32
    objectdist                =33
    objectcgo                 =34
    objectcallback            =35
    objectsurface             =36
    objectgadget              =37
    objectslice               =38
    objectvolume              =39

    repangle                  =43
    repdihederal              =44
    repwirebond               =45
    repcylbond                =46
    replabel                  =47
    repsphere                 =49
    repsurface                =50
    repmesh                   =51
    repdot                    =52
    repnonbonded              =53
    repnonbondedsphere        =54
    repdistdash               =55
    repdistlabel              =56
    repribbon                 =57
    repcartoon                =58
    sculpt                    =59
    vfont                     =60
    # in layer0
    shader                    =61
    shadermgr                 =62
    shaderprg                 =63
    session                   =64
    # in layer1
    property                  =65

    executive                 =70
    selector                  =71
    editor                    =72
    nag                       =73

    export                    =75
    ccmd                      =76
    api                       =77

    main                      =80

# This second set, with negative indices
# represent "python-only" subsystems

    parser                    =-1
    cmd                       =-2

class fb_mask:
    output =              0x01 # Python/text output
    results =             0x02
    errors =              0x04
    actions =             0x08
    warnings =            0x10
    details =             0x20
    blather =             0x40
    debugging =           0x80
    everything =          0xFF


# State count starts at 1. Values <=0 have special meaning.
ALL_STATES = 0
CURRENT_STATE = -1
