include Rules.make

all: pymol

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

.depends: 
	/bin/rm -f .includes
	cd layer0;$(MAKE) depends
	cd layer1;$(MAKE) depends
	cd layer2;$(MAKE) depends
	cd layer3;$(MAKE) depends
	cd layer4;$(MAKE) depends
	cd layer5;$(MAKE) depends

pymol: .includes .depends .update 
	/bin/rm -f .update .includes
	cc $(BUILD) */*.o $(CFLAGS)  $(LIB_DIRS) $(LIBS)

fast: .update
	/bin/rm -f .update 
	cc */*.o $(CFLAGS) -o pymol $(LIB_DIRS) $(LIBS)

depends: 
	/bin/rm -f */*.p
	$(MAKE) .depends

clean: 
	/bin/rm -f *.log core */core game.* log.* layer*/*.o layer*/*.p .update layer*/.files layer*/.depends layer*/.includes 

distclean: clean
	/bin/rm -f modules/*.pyc modules/_pm.so pymol.exe
	/bin/rm -f modules/Pmw/*.pyc modules/Pmw/*/*.pyc modules/Pmw/*/*/*.pyc

dist: distclean
	cd ..;tar -cvf - pymol | gzip > pymol.tgz


