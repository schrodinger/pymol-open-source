

/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/

#include"os_predef.h"
#include"os_std.h"

#include"Raw.h"
#include"Version.h"
#include"Base.h"
#include"Feedback.h"
#include"MemoryDebug.h"
#include"OOMac.h"

#define cRaw_file_stream 0

#define cRaw_header_size 16

#define swap_bytes(a) { \
char u; \
u=*(((char*)(a))); \
*(((char*)(a)))=*(((char*)(a))+3); \
*(((char*)(a))+3)=u; \
u=*(((char*)(a))+1); \
*(((char*)(a))+1)=*(((char*)(a))+2); \
*(((char*)(a))+2)=u; }

CRaw *RawOpenRead(PyMOLGlobals * G, char *fname)
{
  int target = 0x04030201;
  int reverse = 0x01020304;
  int actual;
  int ok = true;

  OOAlloc(G, CRaw);
  I->bufVLA = NULL;
  I->G = G;
  I->f = fopen(fname, "rb");
  if(!I->f) {
    ok = false;
  } else {
    if(feof(I->f))
      ok = false;
    else if(fread(&actual, 4, 1, I->f) != 1)
      ok = false;
    else if(actual == target)
      I->swap = false;
    else if(actual == reverse)
      I->swap = true;
    else {
      PRINTFB(G, FB_Raw, FB_Errors)
        "Error-RawOpenRead: Unrecognized byte ordering. This may not a PyMOL file.\n"
        ENDFB(G);
      ok = false;
    }
  }
  if(!ok) {
    if(I->f)
      fclose(I->f);
    OOFreeP(I);
    PRINTFB(G, FB_Raw, FB_Errors)
      "Error-RawOpenRead: Unable to open '%s'.\n", fname ENDFB(G);

  } else {
    I->mode = cRaw_file_stream;
  }
  return (I);
}

CRaw *RawOpenWrite(PyMOLGlobals * G, char *fname)
{
  int target = 0x04030201;
  int ok = true;
  OOAlloc(G, CRaw);
  I->bufVLA = NULL;
  I->G = G;
  I->f = fopen(fname, "wb");
  if(!I->f) {
    ok = false;
  } else {
    fwrite(&target, 4, 1, I->f);
  }
  if(!ok) {
    if(I->f)
      fclose(I->f);
    OOFreeP(I);
  } else {
    I->mode = cRaw_file_stream;
  }
  return (I);
}

CRaw *RawOpenAppend(PyMOLGlobals * G, char *fname)
{
  int target = 0x04030201;
  int ok = true;
  OOAlloc(G, CRaw);
  I->bufVLA = NULL;
  I->G = G;
  I->f = fopen(fname, "wba");
  if(!I->f) {
    ok = false;
  } else {
    if(!ftell(I->f))            /* write magic if this is a new file */
      fwrite(&target, 4, 1, I->f);
  }
  if(!ok) {
    if(I->f)
      fclose(I->f);
    OOFreeP(I);
    PRINTFB(G, FB_Raw, FB_Errors)
      "Error-RawOpenAppend: Unable to open '%s'.\n", fname ENDFB(G);
  } else {
    I->mode = cRaw_file_stream;
  }
  return (I);
}

void RawFree(CRaw * I)
{
  switch (I->mode) {
  case cRaw_file_stream:
    if(I->f) {
      fclose(I->f);
      I->f = NULL;
    }
    break;
  }
  VLAFreeP(I->bufVLA);
  OOFreeP(I);
}

int RawGetNext(CRaw * I, int *size, int *version)
{
  PyMOLGlobals *G = I->G;
  int result = cRaw_EOF;
  switch (I->mode) {
  case cRaw_file_stream:
    if(I->f) {
      if(!feof(I->f)) {
        if(fread((char *) I->header, cRaw_header_size, 1, I->f) != 1) {
          PRINTFD(G, FB_Raw)
            " RawGetNextType-Debug: Couldn't read header.\n" ENDFD;
        } else {
          if(I->swap) {
            swap_bytes(I->header);
            swap_bytes(I->header + 1);
            swap_bytes(I->header + 2);
            swap_bytes(I->header + 3);
          }
          fseek(I->f, -cRaw_header_size, SEEK_CUR);
          *size = I->header[0];
          result = I->header[1];
          *version = I->header[2];
        }
      }
    }
  }
  return (result);
}

int RawReadSkip(CRaw * I)
{
  PyMOLGlobals *G = I->G;
  int result = false;
  switch (I->mode) {
  case cRaw_file_stream:
    if(I->f) {
      if(!feof(I->f)) {
        if(fread((char *) I->header, cRaw_header_size, 1, I->f) != 1) {
          PRINTFB(G, FB_Raw, FB_Errors)
            "Error-Raw: Error reading header.\n" ENDFB(G);
        } else {
          if(I->swap) {
            swap_bytes(I->header);
            swap_bytes(I->header + 1);
            swap_bytes(I->header + 2);
            swap_bytes(I->header + 3);
          }
          fseek(I->f, I->header[0], SEEK_CUR);
          result = true;
        }
      }
    }
    break;
  }
  return (result);
}

char *RawRead(CRaw * I, int *type, unsigned int *size, int *serial)
{
  PyMOLGlobals *G = I->G;
  char *result = NULL;
  switch (I->mode) {
  case cRaw_file_stream:
    if(I->f) {
      if(!feof(I->f)) {
        if(fread((char *) I->header, cRaw_header_size, 1, I->f) != 1) {
          PRINTFB(G, FB_Raw, FB_Errors)
            "Error-Raw: Error reading header.\n" ENDFB(G);
        } else {
          if(I->swap) {
            swap_bytes(I->header);
            swap_bytes(I->header + 1);
            swap_bytes(I->header + 2);
            swap_bytes(I->header + 3);
          }
          VLACheck(I->bufVLA, char, I->header[0]);
          if(fread(I->bufVLA, I->header[0], 1, I->f) == 1) {
            result = I->bufVLA;
            *size = I->header[0];
            *type = I->header[1];       /* record type */
            *serial = I->header[3];
          } else {
            PRINTFB(G, FB_Raw, FB_Errors)
              "Error-RawRead: Data read error.\n" ENDFB(G);

          }
        }
      } else {
        *type = cRaw_EOF;
      }
    }
    break;
  }
  return (result);
}

char *RawReadPtr(CRaw * I, int type, int *size)
{
  PyMOLGlobals *G = I->G;
  char *result = NULL;
  switch (I->mode) {
  case cRaw_file_stream:
    if(I->f) {
      if(!feof(I->f)) {
        if(fread((char *) I->header, cRaw_header_size, 1, I->f) != 1) {
          PRINTFB(G, FB_Raw, FB_Errors)
            "Error-Raw: Error reading header.\n" ENDFB(G);
        } else {
          if(I->swap) {
            swap_bytes(I->header);
            swap_bytes(I->header + 1);
            swap_bytes(I->header + 2);
            swap_bytes(I->header + 3);
          }
          if((I->header[1]) != type) {
            fseek(I->f, -cRaw_header_size, SEEK_CUR);
            PRINTFD(G, FB_Raw)
              " RawReadPtr-Debug: Type mismatch.\n" ENDFD;
          } else {
            result = mmalloc(I->header[0]);
            if(fread(result, I->header[0], 1, I->f) != 1) {
              FreeP(result);
              PRINTFB(G, FB_Raw, FB_Errors)
                "Error-RawReadVLA: Data read error.\n" ENDFB(G);

            } else {
              *size = I->header[0];
            }
          }
        }
      }
    }
    break;
  }
  return (result);
}

char *RawReadVLA(CRaw * I, int type, unsigned int rec_size, int grow_factor,
                 int auto_zero)
{
  PyMOLGlobals *G = I->G;
  char *result = NULL;
  switch (I->mode) {
  case cRaw_file_stream:
    if(I->f) {
      if(!feof(I->f)) {
        if(fread((char *) I->header, cRaw_header_size, 1, I->f) != 1) {
          PRINTFB(G, FB_Raw, FB_Errors)
            "Error-Raw: Error reading header.\n" ENDFB(G);
        } else {
          if(I->swap) {
            swap_bytes(I->header);
            swap_bytes(I->header + 1);
            swap_bytes(I->header + 2);
            swap_bytes(I->header + 3);
          }
          if((I->header[1]) != type) {
            fseek(I->f, -cRaw_header_size, SEEK_CUR);
            PRINTFD(G, FB_Raw)
              " RawReadVLA-Debug: Type mismatch %d != %d.\n", I->header[1], type ENDFD;

          } else {
            result =
              VLAMalloc((I->header[0] / rec_size), rec_size, grow_factor, auto_zero);
            if(fread(result, I->header[0], 1, I->f) != 1) {
              VLAFreeP(result);
              PRINTFB(G, FB_Raw, FB_Errors)
                "Error-RawReadVLA: Data read error.\n" ENDFB(G);

            } else {
              result = (char *) VLASetSize(result, I->header[0] / rec_size);
            }
          }
        }
      }
    }
    break;
  }
  return (result);
}

int RawReadInto(CRaw * I, int type, unsigned int size, char *buffer)
{
  PyMOLGlobals *G = I->G;
  int ok = false;
  switch (I->mode) {
  case cRaw_file_stream:
    if(I->f) {
      if(!feof(I->f)) {
        if(fread((char *) I->header, cRaw_header_size, 1, I->f) != 1) {
          PRINTFB(G, FB_Raw, FB_Errors)
            "Error-RawReadInfo: Error reading header.\n" ENDFB(G);
        } else {
          if(I->swap) {
            swap_bytes(I->header);
            swap_bytes(I->header + 1);
            swap_bytes(I->header + 2);
            swap_bytes(I->header + 3);
          }
          if((I->header[1]) != type) {
            fseek(I->f, -cRaw_header_size, SEEK_CUR);
            PRINTFD(G, FB_Raw)
              " RawReadPtr-Debug: Type mismatch.\n" ENDFD;
          } else if(I->header[0] != (signed) size) {
            PRINTFB(G, FB_Raw, FB_Errors)
              "Error-RawReadInfo: Size mismatch %d!=%d (disk/RAM).\n", I->header[0], size
              ENDFB(G);
          } else {
            if(fread(buffer, size, 1, I->f) != 1) {
              PRINTFB(G, FB_Raw, FB_Errors)
                "Error-RawReadInfo: Data read error.\n" ENDFB(G);
            } else {
              ok = true;
            }
          }
        }
      }
    }
    break;
  }
  return (ok);
}

int RawWrite(CRaw * I, int type, unsigned int size, int serial, char *bytes)
{
  PyMOLGlobals *G = I->G;
  int header[4];
  int ok = false;
  PRINTFD(G, FB_Raw)
    " RawWrite-Debug: type %d size %d %p\n", type, size, bytes ENDFD;
  switch (I->mode) {
  case cRaw_file_stream:
    if(I->f) {
      header[0] = size;
      header[1] = type;
      header[2] = _PyMOL_VERSION_int;
      header[3] = serial;
      if(fwrite((char *) header, cRaw_header_size, 1, I->f) != 1) {
        PRINTFB(G, FB_Raw, FB_Errors)
          "Error-RawWrite: can't write header.\n" ENDFB(G);
      } else if(fwrite((char *) bytes, size, 1, I->f) != 1) {
        PRINTFB(G, FB_Raw, FB_Errors)
          "Error-RawWrite: can't write data.\n" ENDFB(G);
      } else {
        ok = true;
      }
    }
  }
  PRINTFD(G, FB_Raw)
    " RawWrite-Debug: leaving... %d\n", ok ENDFD;

  return (ok);
}
