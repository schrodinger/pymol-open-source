include ../Rules.make

all: make-modules make-sglite make-pyopengl

make-modules:
	cd modules;echo installdir=$(PYMOL_PATH)/modules > Makefile.pre.in 
	cd modules;echo PYTHON=$(PYTHON_EXE) >> Makefile.pre.in 
	cd modules;cat Makefile.pre.in.src >> Makefile.pre.in 
	cd modules;$(MAKE) -f Makefile.pre.in boot
	cd modules;$(MAKE)
	mv modules/*.so ../modules
	cd modules;$(MAKE) distclean
	cd modules;/bin/rm Makefile.pre.in

make-sglite:
	cd sglite;$(MAKE) 
	mv sglite/*.so ../modules

make-pyopengl:
	cd pyopengl;$(MAKE) 
	mv pyopengl/*.so ../modules

clean: 
	cd sglite;$(MAKE) clean

distclean: clean
	cd sglite;$(MAKE) distclean



