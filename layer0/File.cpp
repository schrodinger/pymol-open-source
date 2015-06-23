/*
 * Copyright (c) Schrodinger, LLC.
 *
 * Basic file IO.
 */

#include <stdio.h>
#include <stdlib.h>

#include "File.h"
#include "MemoryDebug.h"

/*
 * Get the size from the current file pointer to the end of the file
 */
static long fgetsize(FILE *fp) {
  long filesize, current = ftell(fp);
  fseek(fp, 0, SEEK_END);
  filesize = ftell(fp);
  fseek(fp, current, SEEK_SET);
  return filesize;
}

/*
 * Allocate memory and read the entire file from the given file pointer.
 * The file size is stored into the size pointer if not NULL.
 */
static char * fgetcontents(FILE *fp, long *size) {
  long filesize = fgetsize(fp);

  char *contents = (char*) mmalloc(filesize + 255);
  if (!contents)
    return NULL;

  if (1 != fread(contents, filesize, 1, fp)) {
    mfree(contents);
    return NULL;
  }

  if (size)
    *size = filesize;

  contents[filesize] = '\0';
  return contents;
}

/*
 * Allocate memory and read the entire file for the given filename.
 * The file size is stored into the size pointer if not NULL.
 */
char * FileGetContents(const char *filename, long *size) {
  char *contents;
  FILE *fp = fopen(filename, "rb");

  if (!fp)
    return NULL;

  contents = fgetcontents(fp, size);
  fclose(fp);
  return contents;
}
