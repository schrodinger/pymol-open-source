SUBDIRS = ov/src layer0 layer1 layer2 layer3 layer4 layer5 contrib/sglite contrib/modules contrib/champ

moduleslibdir = $(PWD)/built
moduleslib_LTLIBRARIES = libExtensionClass.la
libExtensionClass_la_SOURCES =
libExtensionClass_la_LIBADD = contrib/modules/libExtensionClass0.la

champlibdir = $(PWD)/built
champlib_LTLIBRARIES = libchamp.la
libchamp_la_SOURCES =
libchamp_la_LIBADD = contrib/champ/libchamp0.la

sglitelibdir = $(PWD)/built
sglitelib_LTLIBRARIES = libsglite.la
libsglite_la_SOURCES =
libsglite_la_LIBADD = contrib/sglite/libsglite0.la

pymollibdir = $(PWD)/built
pymollib_LTLIBRARIES = libpymol.la
libpymol_la_SOURCES =
libpymol_la_LIBADD = \
ov/src/libov.la \
layer0/liblayer0.la \
layer1/liblayer1.la \
layer2/liblayer2.la \
layer3/liblayer3.la \
layer4/liblayer4.la \
layer5/liblayer5.la

install-exec-hook:
	touch $(srcdir)/modules/pymol/nofail.so $(srcdir)/modules/chempy/champ/nofail.so
	/bin/rm -f $(srcdir)/modules/pymol/*.so $(srcdir)/modules/chempy/champ/*.so
	/bin/ln -s $(DESTDIR)$(pymollibdir)/libpymol.so.0.0.0 $(srcdir)/modules/pymol/_cmd.so
	/bin/ln -s $(DESTDIR)$(pymollibdir)/libsglite.so.0.0.0 $(srcdir)/modules/pymol/sglite.so
	/bin/ln -s $(DESTDIR)$(pymollibdir)/libExtensionClass.so.0.0.0 $(srcdir)/modules/pymol/ExtensionClass.so
	/bin/ln -s $(DESTDIR)$(pymollibdir)/libchamp.so.0.0.0 $(srcdir)/modules/chempy/champ/_champ.so
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

