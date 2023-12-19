

/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* --------------------------------------------------\-----------------
I* Additional authors of this source file include:
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/
#ifndef _H_Texture
#define _H_Texture

#include "PyMOLGlobals.h"

int TextureInit(PyMOLGlobals * G);
void TextureFree(PyMOLGlobals * G);
/**
 * Determines whether the given character is textured.
 * @param char_id character id
 * @param[out] extent texture extent
 */
bool TextureIsCharTextured(PyMOLGlobals* G, int char_id, float* extent);
void TextureInvalidateTextTexture(PyMOLGlobals * G);
void TextureInitTextTexture(PyMOLGlobals * G);
/**
 * Binds the global Text Texture
 */
void TextureBindTexture(PyMOLGlobals* G);
void TextureGetPlacementForNewSubtexture(PyMOLGlobals * G, int new_texture_width, int new_texture_height, int *new_texture_posx, int *new_texture_posy);
int TextureGetTextTextureSize(PyMOLGlobals * G);

#endif
