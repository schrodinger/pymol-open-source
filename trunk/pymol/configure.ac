# Process this file with autoconf to produce a configure script.
AC_INIT([pymol], [1.1.0], [bugs@delsci.com])
AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE
AC_PROG_CC
AC_PROG_INSTALL
AC_DISABLE_STATIC
AC_PROG_LIBTOOL

AC_PATH_XTRA

WLD_PATH_PYTHON( [2.0], [
    echo "$pythoninc"], [
    echo "Missing Python."
    exit -1
    ])

AC_CHECK_LIB(png, main, [AC_DEFINE(_PYMOL_LIBPNG)
                         LIBS="$LIBS -lpng"], [
    echo "Missing libpng."
    exit -1             
    ])

AC_CHECK_LIB(glut, main, [LIBS="$LIBS -lglut"], [
    echo "Missing libglut."
    exit -1             
])           

AC_CHECK_LIB(dl, main, [LIBS="$LIBS -ldl"], [
    echo "Missing libdl."
    exit -1             
    ])

AC_CHECK_FT2([2.0], [AC_DEFINE(_PYMOL_FREETYPE)], [
    echo "Missing freetype2."
    exit -1             
    ])

LIBS="$FT2_LIBS $LIBS"
CFLAGS="$FT2_CFLAGS $CFLAGS"

LIBS="$LIBS $X_LIBS"
INCS="$X_CFLAGS"

SIM_AC_CHECK_OPENGL([OPENGLLIB="true"], AC_MSG_ERROR([OpenGL library not found]))
if test x"$OPENGLLIB" = x"true" ; then
    SIM_AC_CHECK_GLU([GLULIB="true"], AC_MSG_ERROR([GLU library not found]))
fi
LIBS="$GL_LIBS $LIBS"
CPPFLAGS="-I$pythoninc $CPPFLAGS"

AC_CONFIG_HEADERS([config.h])

AC_DEFINE([_PYMOL_MODULE],[],[Build as a module])
AC_DEFINE([_PYMOL_VMD_PLUGINS],[],[Uses VMD Plugins])

AC_CONFIG_FILES([
    Makefile
    ov/src/Makefile
    layer0/Makefile
    layer1/Makefile
    layer2/Makefile
    layer3/Makefile
    layer4/Makefile
    layer5/Makefile
    contrib/champ/Makefile
    contrib/modules/Makefile
    contrib/sglite/Makefile
    contrib/pyopengl/Makefile
    contrib/uiuc/plugins/molfile_plugin/src/Makefile
    ])  
AC_OUTPUT
