
all: layer

include ../Rules.make
include .includes

.c.p:
	$(CC) -M $(CFLAGS) $< | sed 's/$*\\.o[ :]*/& $@/g;s/[^ ]*\/usr\/[^ ]*//g;s/[^ ]*\/ext\/[^ ]*//g' > $@

layer: $(DEPS) $(OBJS)

depends: $(DEPS)
	cat *.p > .depends;cat .files .depends > .includes;










