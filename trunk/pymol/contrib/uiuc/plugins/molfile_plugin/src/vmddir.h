/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2006 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: vmddir.h,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.7 $       $Date: 2006/01/05 00:05:55 $
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#if defined(_MSC_VER)
#include <windows.h>

typedef struct {
  HANDLE h;
  LPWIN32_FIND_DATAA fd;
} VMDDIR;

#else
#include <dirent.h>

typedef struct {
  DIR * d;
} VMDDIR;
#endif


static VMDDIR * vmd_opendir(const char *);
static char * vmd_readdir(VMDDIR *);
static void vmd_closedir(VMDDIR *);
static int vmd_file_is_executable(const char * filename);


#define VMD_FILENAME_MAX 1024

#if defined(_MSC_VER) 
/* Windows version */

static VMDDIR * vmd_opendir(const char * filename) {
  VMDDIR * d;
 char dirname[VMD_FILENAME_MAX];

  strcpy(dirname, filename);
  strcat(dirname, "\\*");
  d = (VMDDIR *) malloc(sizeof(VMDDIR));
  if (d != NULL) {
    d->h = FindFirstFileA((char*)dirname, d->fd);
    if (d->h == ((HANDLE)(-1))) {
      free(d);
      return NULL;
    }
  }
  return d;
}

static char * vmd_readdir(VMDDIR * d) {
  if (FindNextFileA(d->h, d->fd)) {
    return d->fd->cFileName; 
  }
  return NULL;     
}

static void vmd_closedir(VMDDIR * d) {
  if (d->h != NULL) {
    FindClose(d->h);
  }
  free(d);
}


static int vmd_file_is_executable(const char * filename) {
  FILE * fp;
  if ((fp=fopen(filename, "rb")) != NULL) {
    fclose(fp);
    return 1;
  }

  return 0;
} 

#else

/* Unix version */

#include <sys/types.h>
#include <sys/stat.h>

static VMDDIR * vmd_opendir(const char * filename) {
  VMDDIR * d;

  d = (VMDDIR *) malloc(sizeof(VMDDIR));
  if (d != NULL) {
    d->d = opendir(filename);
    if (d->d == NULL) {
      free(d);
      return NULL;
    }
  }

  return d;
}

static char * vmd_readdir(VMDDIR * d) {
  struct dirent * p;
  if ((p = readdir(d->d)) != NULL) {
    return p->d_name;
  }

  return NULL;     
}

static void vmd_closedir(VMDDIR * d) {
  if (d->d != NULL) {
    closedir(d->d);
  }
  free(d);
}


static int vmd_file_is_executable(const char * filename) {
  struct stat buf;
  if (!stat(filename, &buf)) {
    if (buf.st_mode & S_IXUSR || 
        buf.st_mode & S_IXGRP ||
        buf.st_mode & S_IXOTH) {
      return 1;
    }
  }
  return 0;
} 

#endif




