include ../../Rules.make

CFLAGS = $(CCOPT1) $(CCOPT2) -I../../ext/include/python1.5 -I../../ext/include -DHAVE_ARRAYOBJECT_H $(DEFS) -I/usr/X11R6/include

GLLIBS = -L../../ext/lib $(XLIB_DIR) -lglut -lGL -lGLU -lXmu -lXi
all: \
_glu_nummodule.so \
_glumodule.so \
_glutmodule.so \
_opengl_nummodule.so \
_openglmodule.so \
openglutil.so \
openglutil_num.so                                                                  

.c.o: 
	$(CC) $(CFLAGS) -c $*.c -DPythonTypes

_glu_nummodule.so: _glu_nummodule.o
	$(CC) $(CCOPT1) -shared _glu_nummodule.o -o _glu_nummodule.so $(GLLIBS)

_glumodule.so: _glumodule.o
	$(CC) $(CCOPT1) -shared _glumodule.o -o _glumodule.so $(GLLIBS)

_glutmodule.so: _glutmodule.o
	$(CC) $(CCOPT1) -shared _glutmodule.o -o _glutmodule.so $(GLLIBS) 

_opengl_nummodule.so: _opengl_nummodule.o
	$(CC) $(CCOPT1) -shared _opengl_nummodule.o -o _opengl_nummodule.so $(GLLIBS)

_openglmodule.so: _openglmodule.o
	$(CC) $(CCOPT1) -shared _openglmodule.o -o _openglmodule.so $(GLLIBS)

openglutil.so: openglutil.o
	$(CC) $(CCOPT1) -shared openglutil.o -o openglutil.so $(GLLIBS)

openglutil_num.so: openglutil_num.o
	$(CC) $(CCOPT1) -shared openglutil_num.o -o openglutil_num.so $(GLLIBS)      
                                                          
clean:
	touch .no_fail
	/bin/rm -f *.so *.o .no_fail*

distclean: clean


