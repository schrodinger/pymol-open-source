SUBDIRS = ov/src layer0 layer1 layer2 layer3 layer4 layer5 \
  contrib/uiuc/plugins/molfile_plugin/src \
  contrib/champ contrib/pyopengl 

# contrib/sglite contrib/modules

pymollibdir = $(PWD)/built
pymollib_LTLIBRARIES = libpymol.la
libpymol_la_SOURCES =
libpymol_la_LIBADD = \
contrib/uiuc/plugins/molfile_plugin/src/libmolfile_plugin0.la \
ov/src/libov.la \
layer0/liblayer0.la \
layer1/liblayer1.la \
layer2/liblayer2.la \
layer3/liblayer3.la \
layer4/liblayer4.la \
layer5/liblayer5.la 

champlibdir = $(PWD)/built
champlib_LTLIBRARIES = libchamp.la
libchamp_la_SOURCES =
libchamp_la_LIBADD = contrib/champ/libchamp0.la

#moduleslibdir = $(PWD)/built
#moduleslib_LTLIBRARIES =  libExtensionClass.la
#libExtensionClass_la_SOURCES =
#libExtensionClass_la_LIBADD =  contrib/modules/libExtensionClass0.la

#sglitelibdir = $(PWD)/built
#sglitelib_LTLIBRARIES = libsglite.la
#libsglite_la_SOURCES =
#libsglite_la_LIBADD = contrib/sglite/libsglite0.la

pyopengllibdir = $(PWD)/built
pyopengllib_LTLIBRARIES = \
    libopenglutil.la \
    libopenglutil_num.la \
    lib_openglmodule.la \
    lib_opengl_nummodule.la \
    lib_glutmodule.la \
    lib_glumodule.la \
    lib_glu_nummodule.la 
libopenglutil_num_la_SOURCES = 
libopenglutil_num_la_LIBADD = contrib/pyopengl/libopenglutil_num0.la
libopenglutil_la_SOURCES = 
libopenglutil_la_LIBADD = contrib/pyopengl/libopenglutil0.la
lib_openglmodule_la_SOURCES = 
lib_openglmodule_la_LIBADD = contrib/pyopengl/lib_openglmodule0.la
lib_opengl_nummodule_la_SOURCES = 
lib_opengl_nummodule_la_LIBADD = contrib/pyopengl/lib_opengl_nummodule0.la
lib_glutmodule_la_SOURCES = 
lib_glutmodule_la_LIBADD = contrib/pyopengl/lib_glutmodule0.la
lib_glumodule_la_SOURCES =  
lib_glumodule_la_LIBADD =  contrib/pyopengl/lib_glumodule0.la
lib_glu_nummodule_la_SOURCES = 
lib_glu_nummodule_la_LIBADD = contrib/pyopengl/lib_glu_nummodule0.la

install-exec-hook:
	touch $(srcdir)/modules/pymol/nofail.so $(srcdir)/modules/chempy/champ/nofail.so \
    $(srcdir)/modules/pymol/opengl/gl/nofail.so 
	/bin/rm -f $(srcdir)/modules/pymol/*.so $(srcdir)/modules/chempy/champ/*.so \
      $(srcdir)/modules/pymol/opengl/*/*.so
	/bin/ln -s $(DESTDIR)$(pymollibdir)/libpymol.so.0.0.0 $(srcdir)/modules/pymol/_cmd.so
#	/bin/ln -s $(DESTDIR)$(pymollibdir)/libsglite.so.0.0.0 $(srcdir)/modules/pymol/sglite.so
#	/bin/ln -s $(DESTDIR)$(pymollibdir)/libExtensionClass.so.0.0.0 $(srcdir)/modules/pymol/ExtensionClass.so
	/bin/ln -s $(DESTDIR)$(pymollibdir)/libchamp.so.0.0.0 $(srcdir)/modules/chempy/champ/_champ.so
	/bin/ln -s $(DESTDIR)$(pymollibdir)/libopenglutil.so.0.0.0 $(srcdir)/modules/pymol/opengl/gl/openglutil.so
	/bin/ln -s $(DESTDIR)$(pymollibdir)/libopenglutil_num.so.0.0.0 $(srcdir)/modules/pymol/opengl/gl/openglutil_num.so
	/bin/ln -s $(DESTDIR)$(pymollibdir)/lib_openglmodule.so.0.0.0 $(srcdir)/modules/pymol/opengl/gl/_openglmodule.so
	/bin/ln -s $(DESTDIR)$(pymollibdir)/lib_opengl_nummodule.so.0.0.0 $(srcdir)/modules/pymol/opengl/gl/_opengl_nummodule.so
	/bin/ln -s $(DESTDIR)$(pymollibdir)/lib_glumodule.so.0.0.0 $(srcdir)/modules/pymol/opengl/glu/_glumodule.so
	/bin/ln -s $(DESTDIR)$(pymollibdir)/lib_glu_nummodule.so.0.0.0 $(srcdir)/modules/pymol/opengl/glu/_glu_nummodule.so
	/bin/ln -s $(DESTDIR)$(pymollibdir)/lib_glutmodule.so.0.0.0 $(srcdir)/modules/pymol/opengl/glut/_glutmodule.so
	echo "#!/bin/sh" > ./pymol
	echo PYMOL_PATH=$(PWD) >> ./pymol
	echo export PYMOL_PATH >> ./pymol
	echo exec $(PYTHON) '$$'PYMOL_PATH/modules/launch_pymol.py '"$$@"' >> ./pymol
	chmod 755 ./pymol >> ./pymol
#	/bin/mv $(DESTDIR)$(pymollibdir)/libsglite.so.0.0.0 $(DESTDIR)$(pymollibdir)/sglite.so
#	/bin/rm -f $(DESTDIR)$(pymollibdir)/libsglite* $(DESTDIR)$(pymollibdir)/libpymol*
#	/bin/cp -r $(srcdir)/data $(pkgdatadir)/data
#	$(python_exe) setup.py install --prefix $(DESTDIR)$(exec_prefix)
#	/bin/rm -rf $(srcdir)/build

