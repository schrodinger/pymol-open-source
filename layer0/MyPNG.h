

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
#ifndef _H_MyPNG
#define _H_MyPNG

#include <memory>
#include <vector>

#include "Image.h"

#include"PyMOLGlobals.h"

#define cMyPNG_FormatPNG 0
#define cMyPNG_FormatPPM 1

using png_outbuf_t = std::vector</* png_byte */ unsigned char>;

int MyPNGWrite(const char* file_name, const pymol::Image& img, const float dpi,
    const int format, const int quiet, const float screen_gamma,
    const float file_gamma, png_outbuf_t* io_ptr = nullptr);

std::unique_ptr<pymol::Image> MyPNGRead(const char *file_name);

#endif
