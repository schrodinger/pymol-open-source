#pragma once

#include "os_gl.h"

// Generic Error Handing
bool glCheckOkay();

void gl_debug_proc(GLenum, GLenum, GLuint, GLenum, GLsizei, const GLchar *, const void *);

size_t gl_sizeof(GLenum type);
