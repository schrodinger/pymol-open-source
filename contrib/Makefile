include ../Rules.make

all: make-modules make-sglite make-pyopengl make-champ

static:
	cd modules;$(MAKE) -f Makefile.static
	cd sglite;$(MAKE) static
	cd pyopengl;$(MAKE) static
	touch champ/.includes
	cd champ;$(MAKE) prime
	cd champ;$(MAKE) static

make-modules:
	cd modules;echo installdir=$(PYMOL_PATH)/modules > Makefile.pre.in 
	cd modules;echo PYTHON=$(PYTHON_EXE) >> Makefile.pre.in 
	cd modules;echo PYTHONFRAMEWORK=Python >> Makefile.pre.in
	cd modules;cat Makefile.pre.in.src >> Makefile.pre.in 
	cd modules;$(MAKE) -f Makefile.pre.in boot
	cd modules;$(MAKE)
	mv modules/*.so ../modules/pymol
	cd modules;$(MAKE) distclean
	cd modules;/bin/rm Makefile.pre.in

make-sglite:
	cd sglite;$(MAKE) 
	mv sglite/*.so ../modules/pymol

make-pyopengl:
	cd pyopengl;$(MAKE) 
	mv pyopengl/_glut*.so ../modules/pymol/opengl/glut
	mv pyopengl/_opengl*.so ../modules/pymol/opengl/gl
	mv pyopengl/_glu*.so ../modules/pymol/opengl/glu
	mv pyopengl/opengl*.so ../modules/pymol/opengl/gl

make-champ:
	cd champ;$(MAKE) 
#	mv champ/*.so ../modules/pymol

clean: 
	cd sglite;$(MAKE) clean
	cd champ;$(MAKE) clean
	/bin/rm -f */*.o

distclean: clean
	cd sglite;$(MAKE) distclean
