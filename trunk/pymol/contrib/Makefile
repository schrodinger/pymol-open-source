
all: make-modules make-sglite make-pyopengl

make-modules:
	cd modules;$(MAKE) -f Makefile.pre.in boot
	cd modules;$(MAKE)
	mv modules/*.so ../modules
	cd modules;$(MAKE) distclean

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



