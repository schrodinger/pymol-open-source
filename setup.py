 #!/usr/bin/env python
#
# This script only applies if you are performing a Python Distutils-based
# installation of PyMOL.
#
# It may assume that all of PyMOL's external dependencies are
# pre-installed into the system.

from distutils.core import setup, Extension
import sys, os

#============================================================================
if sys.platform=='win32': 
    # NOTE: this branch not tested in years and may not work...
    inc_dirs=["ov/src",
              "layer0","layer1","layer2",
              "layer3","layer4","layer5",
              "win32/include"]
    libs=["opengl32","glu32","glut32","libpng","zlib"]
    pyogl_libs = ["opengl32","glu32","glut32"]
    lib_dirs=["win32/lib"]
    def_macros=[("_PYMOL_MODULE",None),
                ("WIN32",None),
                ("_PYMOL_LIBPNG",None),
                ]
    ext_comp_args=[]
    ext_link_args=['/NODEFAULTLIB:"LIBC"']
#============================================================================
elif sys.platform=='cygwin':
    # NOTE: this branch not tested in years and may not work...
    inc_dirs=["ov/src",
              "layer0","layer1","layer2",
              "layer3","layer4","layer5",
	      "modules/cealign/src", 
	      "modules/cealign/src/tnt" ]
    libs=["glut32","opengl32","glu32","png"]
    pyogl_libs = ["glut32","opengl32","glu32"]
    lib_dirs=["/usr/lib/w32api"]
    def_macros=[("_PYMOL_MODULE",None),
#                        ("_PYMOL_NUMPY",None),
                ("CYGWIN",None),
                ("_PYMOL_LIBPNG",None)]
    ext_comp_args=[]
    ext_link_args=[]
#============================================================================
elif sys.platform=='darwin':
    using_fink = "'/sw/" in str(sys.path)
    if using_fink:
        # NOTE: this code not tested in years and may not work...
        #
        # under Fink, with the following packages installed:
        #
        #  python24
        #  libpng3
        #  pmw-py24
        #  freetype2
        #  freetype2-dev
        #
        # REMEMBER to use Fink's Python!
        #
        inc_dirs=["ov/src",
                  "layer0","layer1","layer2",
                  "layer3","layer4","layer5", 
                  "/System/Library/Frameworks/OpenGL.framework/Headers",
                  "/System/Library/Frameworks/GLUT.framework/Headers",
                  "/System/Library/Frameworks/CoreFoundation.framework/Headers",
                  "/System/Library/Frameworks/AppKit.framework/Headers",
                  "/System/Library/Frameworks/ApplicationServices.framework/Headers",
                  "/System/Library/Frameworks/Cocoa.framework/Headers",
                  "/System/Library/Frameworks/IOKit.framework/Headers",
                  "/sw/lib/freetype2/include",
                  "/sw/lib/freetype2/include/freetype2",
                  "/sw/include",
		  "modules/cealign/src",
		  "modules/cealign/src/tnt",
		  #"contrib/uiuc/plugins/include/",
		  #"contrib/uiuc/plugins/molfile_plugin/src",
	]
        libs=[]
        pyogl_libs = []
        lib_dirs=[]
        def_macros=[("_PYMOL_MODULE",None),
                    ("_PYMOL_OSX",None),
                    ("_PYMOL_LIBPNG",None),
                    ("_PYMOL_FREETYPE",None),
                    ]
        ext_comp_args=[]
        ext_link_args=["-framework","OpenGL",
                       "-framework","AppKit",
                       "-framework","ApplicationServices",
                       "-framework","CoreFoundation",
                       "-framework","Cocoa",
                       "-framework","IOKit",
                       "-framework","GLUT",
                       "-lpng",
                       "-L/sw/lib/freetype2/lib", "-lfreetype" ]
    else:
        # Not using Fink -- building as if we are on Linux/X11 with
        # the external dependencies compiled into "./ext" in the
        # current working directory,
        #
        # REMEMEBER to use "./ext/bin/python ..."
        #
        EXT = os.getcwd()+"/ext"
        inc_dirs=["ov/src",
                  "layer0","layer1","layer2",
                  "layer3","layer4","layer5", 
                  "/usr/X11R6/include",
                  EXT+"/include",
                  EXT+"/include/GL",
                  EXT+"/include/freetype2",
		  "modules/cealign/src",
		  "modules/cealign/src/tnt",
                  ]
        libs=[]
        pyogl_libs = []
        lib_dirs=[]
        def_macros=[("_PYMOL_MODULE",None),
                    ("_PYMOL_LIBPNG",None),
                    ("_PYMOL_FREETYPE",None),
                    ]
        ext_comp_args=[]
        ext_link_args=[
	               "-L/usr/X11R6/lib64", "-lGL", "-lXxf86vm",
                   "-L"+EXT+"/lib", "-lpng", "-lglut", "-lfreetype"
                   ]
#============================================================================
else: # linux or other unix
    inc_dirs=["ov/src",
              "layer0","layer1","layer2",
              "layer3","layer4","layer5",
              "/usr/include/freetype2",
#                 "/users/warren/ext/include",
# VMD plugin support
#              "contrib/uiuc/plugins/include",
#              "contrib/uiuc/plugins/molfile_plugin/src",
              "modules/cealign/src",
	      "modules/cealign/src/tnt", ]
    libs=["GL","GLU","glut","png","z","freetype", "GLEW",
#          "Xxf86vm"
          ]	
    pyogl_libs = ["GL","GLU","glut", "GLEW"]
    lib_dirs=[
        "/usr/X11R6/lib64",
#        "/users/warren/pymol/ext/lib"
        ]
    def_macros=[("_PYMOL_MODULE",None),
                ("_PYMOL_INLINE",None),
                ("_PYMOL_FREETYPE",None),
                ("_PYMOL_LIBPNG",None),
# OpenGL shaders
               ("_PYMOL_OPENGL_SHADERS",None),
# Numeric Python support                    
#                ("_PYMOL_NUMPY",None),
# VMD plugin support           
#               ("_PYMOL_VMD_PLUGINS",None)
		("NO_MMLIBS",None),
                ]
    ext_comp_args=["-ffast-math","-funroll-loops","-O3"]
    ext_link_args=[]
  
setup ( # Distribution meta-data
    name = "pymol",
	version = "1.4.0", # see layer0/Version.h for updated version
	package_dir = {'' : 'modules'},
	packages = ['chempy',
                'chempy/bmin',
                'chempy/champ',
                'chempy/fast',
                'chempy/fragments',
                'chempy/tinker',
                'pmg_tk',
                'pmg_tk/startup',
                'pmg_tk/skins',
                'pmg_tk/skins/normal',                    
                'pmg_wx',
                'pymol',
                'pymol/contrib',                
                'pymol/opengl',
                'pymol/opengl/gl',
                'pymol/opengl/glu',
                'pymol/opengl/glut',
                'pymol/wizard',
		'web',
		'web/examples',
		'web/javascript', ],
	ext_modules = [
    Extension("pymol._cmd", [
    "modules/cealign/src/ccealignmodule.cpp",
    "ov/src/OVContext.c",
    "ov/src/OVHeapArray.c",
    "ov/src/OVHeap.c",
    "ov/src/OVLexicon.c",
    "ov/src/OVOneToOne.c",
    "ov/src/OVOneToAny.c",
    "ov/src/OVRandom.c",
    "ov/src/ov_utility.c",
    "layer0/Block.c",
    "layer0/Crystal.c",
    "layer0/Debug.c",
    "layer0/Deferred.c",
    "layer0/Err.c",
    "layer0/Feedback.c",
    "layer0/Field.c",
    "layer0/Isosurf.c",
    "layer0/Map.c",
    "layer0/Match.c",
    "layer0/Matrix.c",
    "layer0/MemoryDebug.c",
    "layer0/MemoryCache.c",
    "layer0/MyPNG.c",
    "layer0/Parse.c",
    "layer0/Pixmap.c",
    "layer0/Queue.c",
    "layer0/Raw.c",
    "layer0/Sphere.c",
    "layer0/ShaderMgr.c",
    "layer0/Tetsurf.c",
    "layer0/Texture.c",
    "layer0/Tracker.c",
    "layer0/Triangle.c",
    "layer0/Util.c",
    "layer0/Vector.c",
    "layer0/Word.c",
    "layer0/os_gl.c",
    "layer1/Basis.c",
    "layer1/ButMode.c",
    "layer1/Character.c",
    "layer1/CGO.c",
    "layer1/Color.c",
    "layer1/Control.c",
    "layer1/Extrude.c",
    "layer1/Font.c",
    "layer1/FontType.c",
    "layer1/FontGLUT.c",
    "layer1/FontGLUT8x13.c",
    "layer1/FontGLUT9x15.c",
    "layer1/FontGLUTHel10.c",
    "layer1/FontGLUTHel12.c",
    "layer1/FontGLUTHel18.c",
    "layer1/Movie.c",
    "layer1/Ortho.c",
    "layer1/P.c",
    "layer1/PConv.c",
    "layer1/Pop.c",
    "layer1/PyMOLObject.c",
    "layer1/Ray.c",
    "layer1/Rep.c",
    "layer1/Scene.c",
    "layer1/ScrollBar.c",
    "layer1/Seq.c",
    "layer1/Setting.c",
    "layer1/Shaker.c",
    "layer1/Symmetry.c",
    "layer1/Text.c",
    "layer1/TypeFace.c",
    "layer1/Wizard.c",
    "layer1/View.c",
    "layer2/AtomInfo.c",
    "layer2/CoordSet.c",
    "layer2/GadgetSet.c",    
    "layer2/DistSet.c",
    "layer2/ObjectAlignment.c",
    "layer2/ObjectCGO.c",
    "layer2/ObjectCallback.c",
    "layer2/ObjectDist.c",
    "layer2/ObjectMap.c",
    "layer2/ObjectMesh.c",
    "layer2/ObjectMolecule.c",
    "layer2/ObjectMolecule2.c",
    "layer2/ObjectSurface.c",
    "layer2/ObjectSlice.c",
    "layer2/ObjectVolume.c",
    "layer2/RepCartoon.c",
    "layer2/RepCylBond.c",
    "layer2/RepDistDash.c",
    "layer2/RepDistLabel.c",
    "layer2/RepDot.c",
    "layer2/RepLabel.c",
    "layer2/RepMesh.c",
    "layer2/ObjectGadget.c",
    "layer2/ObjectGadgetRamp.c",
    "layer2/ObjectGroup.c",
    "layer2/RepAngle.c",        
    "layer2/RepDihedral.c",    
    "layer2/RepNonbonded.c",
    "layer2/RepNonbondedSphere.c",
    "layer2/RepRibbon.c",
    "layer2/RepSphere.c",
    "layer2/RepEllipsoid.c",    
    "layer2/RepSurface.c",
    "layer2/RepWireBond.c",
    "layer2/Sculpt.c",
    "layer2/SculptCache.c",
    "layer2/VFont.c",    
    "layer3/PlugIOManager.c",
    "layer3/Editor.c",
    "layer3/Executive.c",
    "layer3/Seeker.c",
    "layer3/Selector.c",
    "layer4/Cmd.c",
    "layer4/Export.c",
    "layer4/Menu.c",
    "layer4/PopUp.c",
    "layer5/PyMOL.c",
    "layer5/TestPyMOL.c",
    "layer5/main.c"
# VMD plugin support
    # switch the 0 to 1 to activate the additional source code
    ] + 0 * [
# (incomplete support -- only TRJ, TRR, XTC, DCD so far...)
    'contrib/uiuc/plugins/molfile_plugin/src/PlugIOManagerInit.c',
    'contrib/uiuc/plugins/molfile_plugin/src/avsplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/basissetplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/bgfplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/binposplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/biomoccaplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/brixplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/carplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/ccp4plugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/corplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/cpmdlogplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/cpmdplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/crdplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/cubeplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/dcdplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/dlpolyplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/dsn6plugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/dtrplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/dxplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/edmplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/fs4plugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/gamessplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/gaussianplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/graspplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/grdplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/gridplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/gromacsplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/hash.c',
    'contrib/uiuc/plugins/molfile_plugin/src/jsplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/maeffplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/mapplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/mdfplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/mmcif.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/mol2plugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/moldenplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/mrcplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/msmsplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/namdbinplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/parm7plugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/parmplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/pbeqplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/pdbplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/phiplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/pltplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/pqrplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/psfplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/raster3dplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/rst7plugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/situsplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/spiderplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/stlplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/tinkerplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/uhbdplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/vaspchgcarplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/vaspoutcarplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/vaspparchgplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/vaspposcarplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/vaspxdatcarplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/vaspxmlplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/vtfplugin.c',
    'contrib/uiuc/plugins/molfile_plugin/src/xbgfplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/xsfplugin.cpp',
    'contrib/uiuc/plugins/molfile_plugin/src/xyzplugin.c'    
    ],
    include_dirs = inc_dirs,
    libraries = libs,
    library_dirs = lib_dirs,
    define_macros = def_macros,
    extra_link_args = ext_link_args,
    extra_compile_args = ext_comp_args,
              ),
# PyMOL no longers uses sglite directly
#
#    Extension("pymol.sglite", [
#    "contrib/sglite/runtests.c",
#    "contrib/sglite/sgcb.c",
#    "contrib/sglite/sgcharmx.c",
#    "contrib/sglite/sgfile.c",
#    "contrib/sglite/sggen.c",
#    "contrib/sglite/sgglobal.c",
#    "contrib/sglite/sghall.c",
#    "contrib/sglite/sghkl.c",
#    "contrib/sglite/sglitemodule.c",
#    "contrib/sglite/sgltr.c",
#    "contrib/sglite/sgmath.c",
#    "contrib/sglite/sgmetric.c",
#    "contrib/sglite/sgnorm.c",
#    "contrib/sglite/sgprop.c",
#    "contrib/sglite/sgss.c",
#    "contrib/sglite/sgstr.c",
#    "contrib/sglite/sgsymbols.c",
#    "contrib/sglite/sgtidy.c",
#    "contrib/sglite/sgtype.c",
#    "contrib/sglite/sgutil.c"
#    ],
#    define_macros=[("PythonTypes",None)],    
#    include_dirs=["contrib/sglite","contrib/modules"]
#                 ),
#    Extension("pymol.ExtensionClass",["contrib/modules/ExtensionClass.c"]),
    
    Extension("chempy.champ._champ", [
    "contrib/champ/champ.c",
    "contrib/champ/champ_module.c",
    "contrib/champ/chiral.c",
    "contrib/champ/err2.c",
    "contrib/champ/feedback2.c",
    "contrib/champ/list.c",
    "contrib/champ/os_memory.c",
    "contrib/champ/sort.c",
    "contrib/champ/strblock.c",
    "contrib/champ/vla.c",
    ],
    include_dirs=["contrib/champ"]
              ),
    Extension("pymol.opengl.glu._glu_num", ["contrib/pyopengl/_glu_nummodule.c"],
              include_dirs = inc_dirs,
              libraries = pyogl_libs,
              library_dirs = lib_dirs,
              define_macros = def_macros,
              ),
    Extension("pymol.opengl.glu._glu", ["contrib/pyopengl/_glumodule.c"],
              include_dirs = inc_dirs,
              libraries = pyogl_libs,
              library_dirs = lib_dirs,
              define_macros = def_macros
              ),
    Extension("pymol.opengl.glut._glut", ["contrib/pyopengl/_glutmodule.c"],
              include_dirs = inc_dirs,
              libraries = pyogl_libs,
              library_dirs = lib_dirs,
              define_macros = def_macros
              ),
    Extension("pymol.opengl.gl._opengl_num", ["contrib/pyopengl/_opengl_nummodule.c"],
              include_dirs = inc_dirs,
              libraries = pyogl_libs,
              library_dirs = lib_dirs,
              define_macros = def_macros
              ),
    Extension("pymol.opengl.gl._opengl", ["contrib/pyopengl/_openglmodule.c"],
              include_dirs = inc_dirs,
              libraries = pyogl_libs,
              library_dirs = lib_dirs,
              define_macros = def_macros
              ),
    Extension("pymol.opengl.gl.openglutil", ["contrib/pyopengl/openglutil.c"],
              include_dirs = inc_dirs,
              libraries = pyogl_libs,
              library_dirs = lib_dirs,
              define_macros = def_macros
              ),
    Extension("pymol.opengl.gl.openglutil_num", ["contrib/pyopengl/openglutil_num.c"],
              include_dirs = inc_dirs,
              libraries = pyogl_libs,
              library_dirs = lib_dirs,
              define_macros = def_macros
              )
])

print '''
 After running:

     python setup.py install
     
 Please run, to complete the installation:

     python setup2.py install

 To uninstall PyMOL, run:

     python setup2.py uninstall
'''
