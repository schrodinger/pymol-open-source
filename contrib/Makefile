
all: make-modules make-sglite

make-modules:
	cd modules;$(MAKE) -f Makefile.pre.in boot
	cd modules;$(MAKE)
	mv modules/*.so ../modules
	cd modules;$(MAKE) distclean

make-sglite:
	cd sglite;$(MAKE) 
	mv sglite/*.so ../modules

clean: 
	cd sglite;$(MAKE) clean

distclean: clean
	cd sglite;$(MAKE) distclean



