/*
 * Copyright (c) Schrodinger, LLC.
 *
 * Basic file IO.
 */

#ifndef _H_File
#define _H_File

#ifdef _WIN32
FILE * pymol_fopen(const char * filename, const char * mode);
#else
#define pymol_fopen fopen
#endif

char * FileGetContents(const char *filename, long *size);

#endif
