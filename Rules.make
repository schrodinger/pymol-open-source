#---------------------------------------------------------------------
# PyMOL Makefile Rules 
#---------------------------------------------------------------------
#
#- Choose One Pair ---------------------------------------------------
#--- Build PyMOL with an embedded Python interpreter (3X performance)
#MODULE = 
#BUILD = -o pymol.exe
#--- Build PyMOL as an importable python module (can use system(),fork())
MODULE = -D_PYMOL_MODULE
BUILD = -shared -o modules/_pm.so
#---------------------------------------------------------------------
#
#- Choose One --------------------------------------------------------
#--- Workaround for XFree86/DRI linux dll problem for module build
BUGS = -D_DRI_WORKAROUND
#---
#BUGS =
#---------------------------------------------------------------------
#
#- Choose One --------------------------------------------------------
#--- Linux
CCOPT1 = -m486 -D__i686__ -ffast-math -Wall -ansi -Wmissing-prototypes
#--- SGI Irix 6.x
#CCOPT1 = -ansi -n32 -woff 1429,1204
#---------------------------------------------------------------------
#
#- Choose One --------------------------------------------------------
#--- Linux Optimized
#CCOPT2 = -O3 -funroll-loops -fomit-frame-pointer
#CCOPT2 = -pg -O3 -funroll-loops
#--- Irix Optimized
#CCOPT2 = -O2
#--- Debugging
CCOPT2 = -g
#---------------------------------------------------------------------
#
#- Choose One Pair ---------------------------------------------------
#--- Python tree in standard system locations
#PYLIB = -L/usr/local/python/lib/python1.5/config 
#PYINC = -I/usr/include/python1.5
#--- Python tree in local user files (shared python in ext/lib)
PYLIB = 
PYINC = -I../ext/include/python1.5
#
#---------------------------------------------------------------------
#
#- Choose One Pair----------------------------------------------------
#--- Libpng2 in system libraries (/usr/include,/usr/lib)
#PNG = -D_HAVE_LIBPNG 
#ZLIB = 
#--- Libpng2 in local user files (pymol/ext/...)
PNG = -D_HAVE_LIBPNG 
ZLIB = -lz
#--- Libpng2 not available
#PNG = 
#ZLIB = 
#---------------------------------------------------------------------
#
#---------------------------------------------------------------------
# No changes normally required below here
#---------------------------------------------------------------------

CFLAGS = $(CCOPT1) $(CCOPT2) -D_PYMOL_THREADS $(INC_DIRS) $(PNG) $(MODULE) $(BUGS)

CC = cc

INC_DIRS = -I../layer0 -I../layer1 -I../layer2 \
	-I../layer3 -I../layer4 -I../layer5 \
	-I../ext/include -I/usr/X11R6/include $(PYINC)

LIB_DIRS = -L./ext/lib -L/usr/X11R6/lib $(GLUT) $(PYLIB)

LIBS = -lpython1.5 -ltk8.0 -ltcl8.0 -lglut -lGL -lGLU -ldl  -lX11 -lXext -lXmu -lXi -lpng\
	$(ZLIB) -lm











