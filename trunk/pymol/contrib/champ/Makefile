
all: module

include ../../Rules.make
include .includes

PRIME = ls *.c | sed 's/.c$$/.o/'| awk 'BEGIN{printf("OBJS=")}{printf("%s ",$$1)}END{print}'>.files;ls *.c | sed 's/.c$$/.p/'| awk 'BEGIN{printf("DEPS=")}{printf("%s ",$$1)}END{print}'>>.files; touch .depends; cat .files .depends > .includes

.includes:
	$(PRIME)

.c.p:
	$(CC) -M $(CFLAGS) $< | sed 's/$*\\.o[ :]*/& $@/g;s/[^ ]*\/usr\/[^ ]*//g;s/[^ ]*\/ext\/[^ ]*//g' > $@

module: $(DEPS) $(OBJS)
	$(CC) $(CCOPT1) -shared *.o -o _champ.so $(GLLIBS)

depends: $(DEPS)
	cat *.p > .depends;cat .files .depends > .includes;

clean:
	/bin/rm -f *.o *.p *.so .depends .inclues .files








