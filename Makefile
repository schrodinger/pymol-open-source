include Rules.make

all: unix contrib

PRIME = ls *.c | sed 's/.c$$/.o/'| awk 'BEGIN{printf("OBJS=")}{printf("%s ",$$1)}END{print}'>.files;ls *.c | sed 's/.c$$/.p/'| awk 'BEGIN{printf("DEPS=")}{printf("%s ",$$1)}END{print}'>>.files; touch .depends; cat .files .depends > .includes

.includes:
	cd layer0;$(PRIME)
	cd layer1;$(PRIME)
	cd layer2;$(PRIME)
	cd layer3;$(PRIME)
	cd layer4;$(PRIME)
	cd layer5;$(PRIME)
	touch .includes

.update:
	cd layer0;$(MAKE)
	cd layer1;$(MAKE)
	cd layer2;$(MAKE)
	cd layer3;$(MAKE)
	cd layer4;$(MAKE)
	cd layer5;$(MAKE)
	touch .update

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

unix: .includes .depends .update 
	/bin/rm -f .update .includes
	cc $(BUILD) */*.o $(CFLAGS)  $(LIB_DIRS) $(LIBS)

windows: .includes .depends .update 
	echo "EXPORTS" > _cmd.def
	nm --demangle --defined-only */*.o | grep ' T ' | sed 's/.* T //' >> _cmd.def 
	dllwrap --dllname _cmd.pyd --driver-name gcc $(BUILD) --def _cmd.def -s  */*.o $(CFLAGS) $(LIB_DIRS) $(LIBS)
	/bin/rm -f .update .includes

fast: .update
	/bin/rm -f .update 
	cc $(BUILD) */*.o $(CFLAGS) $(LIB_DIRS) $(LIBS)

depends: 
	/bin/rm -f */*.p
	$(MAKE) .depends

partial:
	touch layer5/main.c
	touch layer1/P.c
	touch layer4/Cmd.c
	/bin/rm -f modules/_cmd.so pymol.exe
	$(MAKE)

clean: 
	touch .no_fail
	/bin/rm -f layer*/*.o layer*/*.p modules/*/*.pyc modules/*/*/*.pyc \
	layer*/.files layer*/.depends layer*/.includes \
	*.log core */core game.* log.* _cmd.def .update .contrib .no_fail*
	cd contrib;$(MAKE) clean

distclean: clean
	touch .no_fail
	/bin/rm -f modules/*.pyc modules/*.so modules/pymol/*.so pymol.exe \
	modules/*/*.pyc modules/*/*/*.pyc modules/*/*/*/*.pyc .no_fail* test/cmp/*
	cd contrib;$(MAKE) distclean

dist: distclean
	cd ..;tar -cvf - pymol | gzip > pymol.tgz


