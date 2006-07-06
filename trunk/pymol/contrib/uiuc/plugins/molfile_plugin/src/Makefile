
all: library

include ../../../../../Rules.make
include objects.make

CFLAGS=-I../../include $(CCOPT1) $(CCOPT2)
CPPFLAGS=-I../../include  $(CCOPT1) $(CCOPT2)

library: $(OBJS) PlugIOManagerInit.o

static: $(DEPS) $(OBJS)

clean:
	/bin/rm -f *.o 









