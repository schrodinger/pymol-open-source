
include Rules.make

all: unix contrib

PRIME = ls *.c | sed 's/.c$$/.o/'| awk 'BEGIN{printf("OBJS=")}{printf("%s ",$$1)}END{print}'>.files;ls *.c | sed 's/.c$$/.p/'| awk 'BEGIN{printf("DEPS=")}{printf("%s ",$$1)}END{print}'>>.files; touch .depends; cat .files .depends > .includes

MINDEP=$(PYMOL_PATH)/products/unix-mindep
MDP=$(MINDEP)/pymol

.includes:
	cd ov/src;$(PRIME)
	cd layer0;$(PRIME)
	cd layer1;$(PRIME)
	cd layer2;$(PRIME)
	cd layer3;$(PRIME)
	cd layer4;$(PRIME)
	cd layer5;$(PRIME)
	touch .includes

.update:
	cd ov/src;$(MAKE)
	cd layer0;$(MAKE)
	cd layer1;$(MAKE)
	cd layer2;$(MAKE)
	cd layer3;$(MAKE)
	cd layer4;$(MAKE)
	cd layer5;$(MAKE)
	touch .update

ov: 
	cd ov/src;$(MAKE)

0:
	cd layer0;$(MAKE)

1:
	cd layer1;$(MAKE)

2:
	cd layer2;$(MAKE)

3:
	cd layer3;$(MAKE)

4:
	cd layer4;$(MAKE)

5:
	cd layer5;$(MAKE)

.depends: 
	/bin/rm -f .includes
	cd ov/src;$(MAKE) depends
	cd layer0;$(MAKE) depends
	cd layer1;$(MAKE) depends
	cd layer2;$(MAKE) depends
	cd layer3;$(MAKE) depends
	cd layer4;$(MAKE) depends
	cd layer5;$(MAKE) depends

.contrib:
	cd contrib;$(MAKE)
	touch .contrib

contrib: .contrib

lib:  .includes .depends .update 
	/bin/rm -f .update .includes
	ar crv libPyMOL.a */*.o ov/src/*.o 
	ranlib libPyMOL.a

unix: .includes .depends .update 
	/bin/rm -f .update .includes
	$(CC) $(BUILD) $(DEST) */*.o ov/src/*.o $(CFLAGS)  $(LIB_DIRS) $(LIBS)	


semistatic: .includes .depends .update
	/bin/rm -f .update .includes
	cd contrib;$(MAKE) static
	$(CC) $(BUILD) $(DEST) */*.o ov/src/*.o $(CFLAGS) $(LIB_DIRS) $(LIBS)	

unix-mindep-build: semistatic
	$(PYTHON_EXE) modules/compile_pymol.py
	/bin/rm -rf $(MINDEP)
	install -d $(MDP)/ext/lib
	cp -r modules $(MDP)
	cp -r test $(MDP)
	cp -r data $(MDP)	
	cp -r examples $(MDP)
	cp -r scripts $(MDP)
	cp -r pymol.exe $(MDP)
	cp -r ext/lib/python2.3 $(MDP)/ext/lib
	cp -r ext/lib/tcl8.4 $(MDP)/ext/lib
	cp -r ext/lib/tk8.4 $(MDP)/ext/lib
	/bin/rm -f $(MDP)/ext/lib/python2.3/config/libpython2.3.a
	/bin/rm -rf $(MDP)/ext/lib/python2.3/test
	cp LICENSE $(MDP)
	cp README $(MDP)
	cp setup/INSTALL.unix-mindep $(MDP)/INSTALL
	cp setup/setup.sh.unix-mindep $(MDP)/setup.sh
	cd $(MINDEP);chown -R nobody pymol
	cd $(MINDEP);chgrp -R nobody pymol

unix-mindep: unix-mindep-build
	cd $(MINDEP);tar -cvf - pymol | gzip > ../pymol-0_xx-bin-xxxxx-mindep.tgz

unix-mindep-beta: unix-mindep-build
	cp epymol/data/pymol/beta/splash.png $(MDP)/data/pymol/splash.png
	cd $(MINDEP);tar -cvf - pymol | gzip > ../pymol-0_xx-bin-xxxxx-mindep.tgz

unix-helper: unix-mindep-build
	cp setup/setup.sh.unix-helper $(MDP)/setup.sh
	cd $(MINDEP);tar -cvf - pymol | gzip > ../helperpymol-0_xx-bin-xxxxx-mindep.tgz

irix-mindep: semistatic
	$(PYTHON_EXE) modules/compile_pymol.py
	/bin/rm -rf $(MINDEP)
	mkdir products/unix-mindep
	mkdir products/unix-mindep/pymol
	mkdir products/unix-mindep/pymol/ext
	mkdir products/unix-mindep/pymol/ext/lib
	cp -r modules $(MDP)
	cp -r test $(MDP)
	cp -r data $(MDP)	
	cp -r examples $(MDP)
	cp -r pymol.exe $(MDP)
	cp -r ext/lib/python2.3 $(MDP)/ext/lib
	cp -r ext/lib/tcl8.4 $(MDP)/ext/lib
	cp -r ext/lib/tk8.4 $(MDP)/ext/lib
	/bin/rm -f $(MDP)/ext/lib/python2.3/config/libpython2.3.a
	/bin/rm -rf $(MDP)/ext/lib/python2.3/test
	cp LICENSE $(MDP)
	cp README $(MDP)
	cp setup/INSTALL.unix-mindep $(MDP)/INSTALL
	cp setup/setup.sh.unix-mindep $(MDP)/setup.sh
	cd $(MINDEP);chown -R nobody.nobody pymol
	cd $(MINDEP);tar -cvf - pymol |gzip > ../pymol-0_xx-bin-xxxxx-mindep.tgz

# Sharp3D display running under linux

unix-s3d: .includes .depends .update 
	cd sharp3d/src; $(MAKE)
	$(CC) $(BUILD) $(DEST) */*.o ov/src/*.o sharp3d/src/*.o $(CFLAGS) $(LIB_DIRS) $(LIBS) -Lsharp3d/lib -lsgl.1.3 -lspl.1.2 

unix-s3d-semistatic: .includes .depends .update
	/bin/rm -f .update .includes
	cd contrib;$(MAKE) static
	cd sharp3d/src; $(MAKE)
	$(CC) $(BUILD) $(DEST) */*.o ov/src/*.o sharp3d/src/*.o $(CFLAGS) $(LIB_DIRS) $(LIBS) -Lsharp3d/lib -lsgl.1.3 -lspl.1.2 

unix-s3d-build: unix-s3d-semistatic
	$(PYTHON_EXE) modules/compile_pymol.py
	/bin/rm -rf $(MINDEP)
	install -d $(MDP)/ext/lib
	cp -r modules $(MDP)
	cp -r test $(MDP)
	cp -r data $(MDP)	
	cp -r examples $(MDP)
	cp -r scripts $(MDP)
	cp -r pymol.exe $(MDP)
	cp -r ext/lib/python2.3 $(MDP)/ext/lib
	cp -r ext/lib/tcl8.4 $(MDP)/ext/lib
	cp -r ext/lib/tk8.4 $(MDP)/ext/lib
	cp sharp3d/lib/* $(MDP)/ext/lib/
	/bin/rm -f $(MDP)/ext/lib/python2.3/config/libpython2.3.a
	/bin/rm -rf $(MDP)/ext/lib/python2.3/test
	cp LICENSE $(MDP)
	cp README $(MDP)
	cp setup/INSTALL.unix-mindep $(MDP)/INSTALL
	cp setup/setup.sh.unix-s3d $(MDP)/setup.sh
	cd $(MINDEP);chown -R nobody pymol
	cd $(MINDEP);chgrp -R nobody pymol

unix-s3d-product: unix-s3d-build
	cd $(MINDEP);tar -cvf - pymol | gzip > ../pymol-0_xx-bin-sharp3d.tgz

#windows: .includes .depends .update 
#
#	echo "EXPORTS" > _cmd.def
#	nm --demangle --defined-only */*.o | grep ' T ' | sed 's/.* T //' >> _cmd.def 
#	dllwrap --dllname _cmd.pyd --driver-name gcc $(BUILD) --def _cmd.def -o _cmd.pyd */*.o -s --entry _DllMain@12 --target=i386-mingw32 $(LIB_DIRS) $(LIBS) 
#	/bin/rm -f .update .includes

#cygwin: .includes .depends .update
#	$(CC) -shared -W1,--enable-auto-image-base */*.o $(LIB_DIRS) $(LIBS) $(DEST)

fast: .update
	/bin/rm -f .update 
	$(CC) $(BUILD) */*.o ov/src/*.o $(CFLAGS) $(LIB_DIRS) $(LIBS)

depends: 
	/bin/rm -f */*.p
	$(MAKE) .depends

partial:
	touch layer5/main.c
	touch layer1/P.c
	touch layer4/Cmd.c
	/bin/rm -f modules/pymol/_cmd.so pymol.exe
	$(MAKE)

clean: 
	touch .no_fail
	/bin/rm -f layer*/*.o ov/src/*.o layer*/*.p modules/*/*.pyc modules/*/*/*.pyc \
	layer*/.files layer*/.depends layer*/.includes \
	*.log core */core game.* log.* _cmd.def .update .contrib .no_fail*
	cd contrib;$(MAKE) clean

distclean: clean
	touch .no_fail
	/bin/rm -f modules/*.pyc modules/*.so modules/*/*.so modules/*/*/*.so \
	modules/*/*/*/*.so pymol.exe \
	modules/*/*.pyc modules/*/*/*.pyc modules/*/*/*/*.pyc .no_fail* test/cmp/*
	/bin/rm -rf build
	/bin/rm -rf products/*.tgz products/unix-mindep
	cd contrib;$(MAKE) distclean

pyclean: clean
	/bin/rm -rf build
	/bin/rm -rf ext/lib/python2.1/site-packages/pymol
	/bin/rm -rf ext/lib/python2.1/site-packages/chempy
	/bin/rm -rf ext/lib/python2.1/site-packages/pmg_tk
	/bin/rm -rf ext/lib/python2.1/site-packages/pmg_wx

dist: distclean
	cd ..;tar -cvf - pymol | gzip > pymol.tgz

pmw: 
	cd modules; gunzip < ./pmg_tk/pmw.tgz | tar xvf -

compileall:
	$(PYTHON_EXE) modules/compile_pymol.py

# Everything below here is for the MacPyMOL Incentive Product
# Compilation of MacPyMOL requires layerOSX source code (closed source)

OSXPROD=products/MacPyMOL.app
OSXVIEWER=products/PyMOLViewer.app
OSXHYBRID=products/PyMOLX11Hybrid.app
OSXHELPER=products/HelperPyMOL.app
OSXFRWK=products/FrameworkPyMOL.app
OSXDEMO=products/MacPyMOL\ Demos
OSXPYMOL=$(OSXPROD)/pymol
OSXEXE=$(OSXPROD)/Contents/MacOS/PyMOL
OSXPY=$(OSXPROD)/py23

osx-wrap:
	/bin/rm -rf $(OSXPYMOL) $(OSXEXE) $(OSXPY)
	$(PYMOL_PATH)/layerOSX/tar -czvf layerOSX/bundle/app.hfstar $(OSXPROD)

osx-unwrap:
	/bin/rm -rf $(OSXPROD)
	$(PYMOL_PATH)/layerOSX/tar -xzvf layerOSX/bundle/app.hfstar

osx-python-framework:
	cc layerOSX/bundle/python.c -o $(OSXEXE) $(DEFS)\
$(PYTHON_INC_DIR) \
-framework CoreFoundation -framework Python -lc -Wno-long-double

osx-python-standalone:
	cc layerOSX/bundle/python.c -o $(OSXEXE) $(DEFS)\
$(PYTHON_INC_DIR) -Lext/lib -Lext/lib/python2.3/config -lpython2.3 \
-framework CoreFoundation -lc -Wno-long-double -D_PYMOL_OSX_PYTHONHOME


osx: 
	cd layerOSX/src; $(MAKE)
	$(MAKE) 

osx-dev: osx
	cp modules/pymol/_cmd.so $(OSXPYMOL)/modules/pymol

osx-pdev:
	/bin/rm -rf $(OSXPYMOL)/modules/pymol
	cp -R modules/pymol $(OSXPYMOL)/modules/pymol
	cp -R layerOSX/plugin/pmg_aqua $(OSXPYMOL)/modules/

osx-product: osx 
	$(PYTHON_EXE) modules/compile_pymol.py
	/bin/rm -rf $(OSXPYMOL)
	install -d $(OSXPYMOL)
	cp -R modules $(OSXPYMOL)/
	cp -R layerOSX/plugin/pmg_aqua $(OSXPYMOL)/modules/
	cp -R test $(OSXPYMOL)/
	cp -R data $(OSXPYMOL)/	
	cp -R scripts $(OSXPYMOL)/	
	cp -R examples $(OSXPYMOL)/
	cp LICENSE $(OSXPYMOL)/
	cp README $(OSXPYMOL)/


osx-standalone: osx-unwrap osx-python-standalone osx-product
	/bin/rm -rf $(OSXPY)
	install -d $(OSXPY)/lib
	cp -R ext/lib/python2.3 $(OSXPY)/lib/

mac-wrap-demos:
	/usr/local/bin/tar -czvf layerOSX/applescript.hfstar $(OSXDEMO)

mac-unwrap-demos:
	/usr/local/bin/tar -xzvf layerOSX/applescript.hfstar

mac-demo-data:
	install -d $(OSXPYMOL)/data/demo
	cp -R demo_data/* $(OSXPYMOL)/data/demo/

mac-demo: osx-standalone mac-demo-data mac-unwrap-demos 

mac-framework: osx-unwrap osx-python-framework osx-product
	/bin/rm -rf $(OSXFRWK)
	/bin/cp -R $(OSXPROD) $(OSXFRWK)
	sed 's/MacPyMOL/FrameworkPyMOL/' $(OSXFRWK)/Contents/Info.plist > $(OSXFRWK)/Contents/Info.plist.tmp
	mv $(OSXFRWK)/Contents/Info.plist.tmp $(OSXFRWK)/Contents/Info.plist
	/bin/rm -r $(OSXFRWK)/Contents/Resources/English.lproj/MainMenu.nib
	/bin/rm -r $(OSXFRWK)/Contents/Resources/English.lproj/MainMenu~.nib

mac: osx-standalone
	/bin/cp layerOSX/bundle/splash.png $(OSXPYMOL)/data/pymol/

mac-helper: mac
	/bin/rm -rf $(OSXHELPER)
	/bin/cp -R $(OSXPROD) $(OSXHELPER)
	sed 's/MacPyMOL/HelperPyMOL/' $(OSXHELPER)/Contents/Info.plist > $(OSXHELPER)/Contents/Info.plist.tmp
	mv $(OSXHELPER)/Contents/Info.plist.tmp $(OSXHELPER)/Contents/Info.plist
	/bin/cp data/pymol/splash.png $(OSXHELPER)/pymol/data/pymol/
	/bin/rm -r $(OSXHELPER)/Contents/Resources/English.lproj/MainMenu.nib
	/bin/rm -r $(OSXHELPER)/Contents/Resources/English.lproj/MainMenu~.nib

mac-viewer: mac
	/bin/rm -rf $(OSXVIEWER)
	/bin/cp -R $(OSXPROD) $(OSXVIEWER)
	sed 's/MacPyMOL/PyMOLViewer/' $(OSXVIEWER)/Contents/Info.plist > $(OSXVIEWER)/Contents/Info.plist.tmp
	mv $(OSXVIEWER)/Contents/Info.plist.tmp $(OSXVIEWER)/Contents/Info.plist
	/bin/cp data/pymol/splash.png $(OSXVIEWER)/pymol/data/pymol/
	/bin/rm -r $(OSXVIEWER)/Contents/Resources/English.lproj/MainMenu.nib
	/bin/rm -r $(OSXVIEWER)/Contents/Resources/English.lproj/MainMenu~.nib

mac-hybrid: mac
	/bin/rm -rf $(OSXHYBRID)
	/bin/cp -R $(OSXPROD) $(OSXHYBRID)
	sed 's/MacPyMOL/PyMOLX11Hybrid/' $(OSXHYBRID)/Contents/Info.plist > $(OSXHYBRID)/Contents/Info.plist.tmp
	mv $(OSXHYBRID)/Contents/Info.plist.tmp $(OSXHYBRID)/Contents/Info.plist
	/bin/cp data/pymol/splash.png $(OSXHYBRID)/pymol/data/pymol/
	/bin/rm -r $(OSXHYBRID)/Contents/Resources/English.lproj/MainMenu.nib
	/bin/rm -r $(OSXHYBRID)/Contents/Resources/English.lproj/MainMenu~.nib
	/bin/cp -r ext/lib/tcl8.4 $(OSXHYBRID)
	/bin/cp -r ext/lib/tk8.4 $(OSXHYBRID)

mac-beta:
	make distclean
	make mac
