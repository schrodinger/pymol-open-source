/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2005 by Warren Lyford Delano of DeLano Scientific. 
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

#include "MemoryDebug.h"
#include "Base.h"

#include "TypeFace.h"

#ifdef _PYMOL_FREETYPE

#include <ft2build.h>
#include FT_FREETYPE_H 

struct _CType { 
  FT_Library library;
};

struct _CTypeFace {
  PyMOLGlobals *G;
  FT_Face Face;
  float LastSize;
};

int TypeFaceCharacterNew(CTypeFace *I,CharFngrprnt *fprnt,float size)
{
  FT_GlyphSlot slot = I->Face->glyph; /* a small shortcut */ 
  if(I->LastSize!=size) {
    I->LastSize = size;
    FT_Set_Char_Size( I->Face, /* handle to face object */
                              0, /* char_width in 1/64th of points */
                              (int)(size*64), /* char_height in 1/64th of points */
                              72, /* horizontal device resolution */
                              72 ); /* vertical device resolution */
    
  }
  if(!FT_Load_Char( I->Face, fprnt->u.i.ch, FT_LOAD_RENDER )) 
    return CharacterNewFromBytemap(I->G,
                                   slot->bitmap.width,
                                   slot->bitmap.rows,
                                   -slot->bitmap.pitch,
                                   slot->bitmap.buffer + ((slot->bitmap.rows-1) * slot->bitmap.pitch),
                                   (float)-slot->bitmap_left,
                                   (float)slot->bitmap.rows-slot->bitmap_top, 
                                   slot->advance.x / 64.0F,
                                   fprnt);
  else
    return 0;
}

CTypeFace *TypeFaceLoad(PyMOLGlobals *G,unsigned char *dat,unsigned int len)
{
  register CType *I = G->Type;
  int ok=true;
  CTypeFace *result = Calloc(CTypeFace,1);
  if(result) {
    FT_Error error = FT_New_Memory_Face( I->library, dat, len, 0, &result->Face );
    result->G = G;
    if ( error ) 
      ok = false;
    else {
      result->LastSize = 12.0F; /* default size */
      error = FT_Set_Char_Size( result->Face, /* handle to face object */
                                0, /* char_width in 1/64th of points */
                                (int)(result->LastSize*64), /* char_height in 1/64th of points */
                                72, /* horizontal device resolution */
                                72 ); /* vertical device resolution */

      if(error) {
        ok = false;
      } else {
        error = FT_Select_Charmap( result->Face, FT_ENCODING_UNICODE );
      }
    }
  }
  if(!ok) {
    FreeP(result);
  }
  return result;
}

float TypeFaceGetKerning(CTypeFace *I,unsigned int last, unsigned int current, float size)
{
  float result = 0.0F;
  FT_UInt glyph_index, previous;
  /*  FT_Bool use_kerning = FT_HAS_KERNING( I->Face );*/
  FT_Bool use_kerning = 1;
  if(I->LastSize!=size) {
    I->LastSize = size;
    FT_Set_Char_Size( I->Face, /* handle to face object */
                      0, /* char_width in 1/64th of points */
                      (int)(size*64), /* char_height in 1/64th of points */
                      72, /* horizontal device resolution */
                      72 ); /* vertical device resolution */
    
  }
  if(use_kerning) {
    previous = FT_Get_Char_Index( I->Face, last );
    glyph_index = FT_Get_Char_Index( I->Face, current);
    if ( previous && glyph_index ) { 
      FT_Vector delta; 
      FT_Get_Kerning( I->Face, previous, glyph_index, FT_KERNING_DEFAULT, &delta ); 
      result = delta.x / 64.0F;
    }
  }
  return result;
}

void TypeFaceFree(CTypeFace *I)
{
  FT_Done_Face( I->Face );
  FreeP(I);
}


int TypeInit(PyMOLGlobals *G)
{
  register CType *I;
  if( (I=(G->Type=Calloc(CType,1)))) {
    FT_Error error = FT_Init_FreeType( &I->library ); 
    return !error;
  }
  return 0;
}

void TypeFree(PyMOLGlobals *G)
{
  register CType *I = G->Type;
  FT_Done_FreeType( I->library );
  FreeP(G->Type);
}

#else

int TypeFaceCharacterNew(CTypeFace *I,
                         CharFngrprnt *fprnt,
                         float size)
{
  return 0;
}

float TypeFaceGetKerning(CTypeFace *I,
                         unsigned int last, 
                         unsigned int current,
                         float size)
{
  return 0.0F;
}

CTypeFace *TypeFaceLoad(PyMOLGlobals *G,unsigned char *dat, unsigned int len)
{
  return NULL;
}
void TypeFaceFree(CTypeFace *face)
{}

int TypeInit(PyMOLGlobals *G) { return 1;}
void TypeFree(PyMOLGlobals *G) {}

#endif

