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
#include <ctype.h>


#undef SG_GLOBAL
#include "sglite.h"


int ParseStrXYZ(const char *StrXYZ, int StopChar,
                T_RTMx *RTMx, int FacR, int FacT)
{
  unsigned int  P_mode;
  int           Row, Column, Sign, Mult, i, n;
  double        ValR[3], ValT, Value, V;
  T_RTMx        Mx[1];
  int           iStrXYZ;

#define cStrXYZ StrXYZ[iStrXYZ]
#define ReturnErr return -(++iStrXYZ)


  for (i = 0; i < 12; i++) RTMx->a[i] = 0;
  for (i = 0; i < 12; i++)   Mx->a[i] = 0;

#define P_Add    0x01u
#define P_Mult   0x02u
#define P_Value  0x04u
#define P_XYZ    0x08u
#define P_Comma  0x10u

  Row    = 0;
  Column = -1;
  Sign   = 1;
  Mult   = 0;
  for (i = 0; i < 3; i++) ValR[i] = 0.;
  ValT   = 0.;
  Value  = 0.;
  P_mode = P_Add | P_Value | P_XYZ;

  for (iStrXYZ = 0;; iStrXYZ++)
  {
    if (cStrXYZ == StopChar || cStrXYZ == '\0' || ! isspace(cStrXYZ))
    {
      switch (cStrXYZ == StopChar ? '\0' : cStrXYZ)
      {
        case '_':
          break;
        case '+': Sign =  1; goto ProcessAdd;
        case '-': Sign = -1;
         ProcessAdd:
          if ((P_mode & P_Add) == 0) ReturnErr;
          if (Column >= 0) ValR[Column] += Value;
          else             ValT         += Value;
          Value = 0.;
          Column = -1;
          Mult = 0;
          P_mode = P_Value | P_XYZ;
          break;
        case '*':
          if ((P_mode & P_Mult) == 0) ReturnErr;
          Mult = 1;
          P_mode = P_Value | P_XYZ;
          break;
        case '/':
        case ':':
          if ((P_mode & P_Mult) == 0) ReturnErr;
          Mult = -1;
          P_mode = P_Value;
          break;
        case '.':
        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
          if ((P_mode & P_Value) == 0) ReturnErr;
          i = 1;
          n = sscanf(&cStrXYZ, "%lf%n", &V, &i);
          iStrXYZ += (i - 1);
          if (n != 1) ReturnErr;
          if (Sign == -1) { V = -V; Sign = 1; }
          if      (Mult ==  1)
            Value *= V;
          else if (Mult == -1) {
            if      (V     != 0.) Value /= V;
            else if (Value != 0.) ReturnErr;
          }
          else
            Value = V;
          P_mode = P_Comma | P_Add | P_Mult | P_XYZ;
          break;
        case 'X':
        case 'x': Column = 0; goto Process_XYZ;
        case 'Y':
        case 'y': Column = 1; goto Process_XYZ;
        case 'Z':
        case 'z': Column = 2;
         Process_XYZ:
          if ((P_mode & P_XYZ) == 0) ReturnErr;
          if (Value == 0.) { Value = Sign; Sign = 1; }
          P_mode = P_Comma | P_Add | P_Mult;
          break;
        case ',':
        case ';':
          if (Row == 2) ReturnErr;
        case '\0':
          if ((P_mode & P_Comma) == 0) ReturnErr;
          if (Column >= 0) ValR[Column] += Value;
          else             ValT         += Value;
          for (i = 0; i < 3; i++)
            if (Discretize(ValR[i], &Mx->s.R[Row * 3 + i], FacR) != 0)
              ReturnErr;
          if (Discretize(ValT, &Mx->s.T[Row], FacT) != 0)
            ReturnErr;
          Row++;
          Column = -1;
          Sign   = 1;
          Mult   = 0;
          for (i = 0; i < 3; i++) ValR[i] = 0.;
          ValT   = 0.;
          Value = 0.;
          P_mode = P_Add | P_Value | P_XYZ;
          break;
        default:
          ReturnErr;
      }
    }

    if (cStrXYZ == StopChar || cStrXYZ == '\0')
      break;
  }

  if (Row != 3) ReturnErr;

  for (i = 0; i < 12; i++) RTMx->a[i] = Mx->a[i];

#undef P_Add
#undef P_Mult
#undef P_Value
#undef P_XYZ
#undef P_Comma

#undef cStrXYZ
#undef ReturnErr

  return ++iStrXYZ;
}


const char *FormatFraction(int nume, int deno, int Decimal,
                           char *Buffer, int SizeBuffer)
{
  int          n, d;
  char         *cp, *cpp;
  static char  StaticBuffer[40];


  if (NULL == Buffer) {
              Buffer =        StaticBuffer;
          SizeBuffer = sizeof StaticBuffer / sizeof (*StaticBuffer);
  }

  Buffer[SizeBuffer - 1] = '\0';

  if (nume == 0)
  {
    Buffer[0] = '0';
    Buffer[1] = '\0';
  }
  if (Decimal)
  {
    (void) sprintf(Buffer, "%.6g", (double) nume / deno);

         cp = Buffer;
    if (*cp == '-') cp++;
    if (*cp == '0') {
      cpp = cp + 1; while (*cp) *cp++ = *cpp++;
    }
  }
  else
  {
    SimplifyFraction(nume, deno, &n, &d);

    if (d == 1)
      (void) sprintf(Buffer, "%d", n);
    else
      (void) sprintf(Buffer, "%d/%d", n, d);
  }

  if (Buffer[SizeBuffer - 1] != '\0') {
      Buffer[SizeBuffer - 1] =  '\0';
    SetSgError("Internal Error: FormatFraction(): Buffer too small");
    return NULL;
  }

  return Buffer;
}


const char *RTMx2XYZ(const T_RTMx *RTMx, int RBF, int TBF,
                     int Decimal, int TrFirst, int Low,
                     const char *Separator,
                     char *BufferXYZ, int SizeBufferXYZ)
{
  static const char *UpperXYZ = "XYZ";
  static const char *LowerXYZ = "xyz";

  int         i, j, p, iRo, iTr;
  char        *xyz, buf_tr[32];
  const char  *sep, *LetterXYZ, *ro, *tr;

  static char  StaticBufferXYZ[80];


  if (NULL == BufferXYZ) {
              BufferXYZ  =        StaticBufferXYZ;
          SizeBufferXYZ  = sizeof StaticBufferXYZ / sizeof (*StaticBufferXYZ);
  }

  BufferXYZ[SizeBufferXYZ - 1] = '\0';

  if (Low)
    LetterXYZ = LowerXYZ;
  else
    LetterXYZ = UpperXYZ;

  if (Separator == NULL)
      Separator = ",";

  xyz = BufferXYZ;

  for (i = 0; i < 3; i++)
  {
    if (i != 0)
      for (sep = Separator; *sep; sep++) *xyz++ = *sep;

    sep = xyz;
                        iTr = RTMx->s.T[i];
    tr = FormatFraction(iTr, TBF, Decimal,
                        buf_tr, sizeof buf_tr / sizeof (*buf_tr));
    if (tr == NULL)
      return NULL;

    p = 0;

    if (  TrFirst && iTr) {
      if (*tr) p = 1;
      while (*tr) *xyz++ = *tr++;
    }

    for (j = 0; j < 3; j++)
    {
          iRo = RTMx->s.R[i * 3 + j];
      if (iRo)
      {
            ro = FormatFraction(iRo, RBF, Decimal, NULL, 0);
        if (ro == NULL)
          return NULL;

        if      (*ro == '-')
          *xyz++ = *ro++;
        else if (*ro && p)
          *xyz++ = '+';

        if (ro[0] != '1' || ro[1] != '\0') {
          while (*ro) *xyz++ = *ro++;
          *xyz++ = '*';
        }

        *xyz++ = LetterXYZ[j];

        p = 1;
      }
    }

    if (! TrFirst && iTr)
    {
      if (*tr && *tr != '-' && p)
        *xyz++ = '+';

      while (*tr) *xyz++ = *tr++;
    }

    if (xyz == sep)
      *xyz++ = '0';
  }

  *xyz = '\0';

  if (BufferXYZ[SizeBufferXYZ - 1] != '\0') {
      BufferXYZ[SizeBufferXYZ - 1] =  '\0';
    SetSgError("Internal Error: RTMx2XYZ(): BufferXYZ too small");
    return NULL;
  }

  return BufferXYZ;
}
