/* $Id$ */

/* The source code contained in this file is            */
/* Copyright (C) 1994-2000 by Ralf W. Grosse-Kunstleve. */
/* Please see the LICENSE file for more information.    */

#ifdef _PYMOL_WIN32
#include"os_predef.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SG_GLOBAL
#include "sglite.h"
#include "sgconst.h"
#include "sgrefset.h"
#undef SG_GLOBAL


void SetSgError(const char *msg)
{
  if (SgError == NULL) SgError = msg;
}

void ClrSgError(void)
{
  SgError = NULL;
}


int SetSg_InternalError(int status, const char *file, const int line)
{
  if (SgError == NULL)
  {
    sprintf(SgErrorBuffer, "Internal Error: file \"%s\", line %d",
      file, line);

    SetSgError(SgErrorBuffer);
  }

  return status;
}


const void *pSetSg_InternalError(const void *ptr,
                                 const char *file, const int line)
{
  if (SgError == NULL)
  {
    sprintf(SgErrorBuffer, "Internal Error: file \"%s\", line %d",
      file, line);

    SetSgError(SgErrorBuffer);
  }

  return ptr;
}


int SetSg_NotEnoughCore(int status, const char *file, const int line)
{
  if (SgError == NULL)
  {
    sprintf(SgErrorBuffer, "Error: Not enough core: file \"%s\", line %d",
      file, line);

    SetSgError(SgErrorBuffer);
  }

  return status;
}


const void *pSetSg_NotEnoughCore(const void *ptr,
                                 const char *file, const int line)
{
  if (SgError == NULL)
  {
    sprintf(SgErrorBuffer, "Error: Not enough core: file \"%s\", line %d",
      file, line);

    SetSgError(SgErrorBuffer);
  }

  return ptr;
}
