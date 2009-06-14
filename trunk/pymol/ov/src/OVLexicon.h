
#ifndef _H_OVLexicon
#define _H_OVLexicon

#include "OVHeap.h"
#include "OVreturns.h"


/* 
   OVLexicon -- a collection of strings and their identifiers
*/

#ifndef OVLexicon_DEFINED
typedef struct _OVLexicon OVLexicon;
#define OVLexicon_DEFINED
#endif

OVLexicon *OVLexicon_New(OVHeap * heap);
void OVLexicon_Del(OVLexicon * I);

#define OVLexicon_DEL_AUTO_NULL(I) { if(I) { OVLexicon_Del(I); I=OV_NULL; }}

OVreturn_word OVLexicon_GetFromCString(OVLexicon * uk, ov_char8 * str);

OVstatus OVLexicon_IncRef(OVLexicon * uk, ov_word id);
OVstatus OVLexicon_DecRef(OVLexicon * uk, ov_word id);

ov_char8 *OVLexicon_FetchCString(OVLexicon * uk, ov_word id);
OVreturn_word OVLexicon_BorrowFromCString(OVLexicon * uk, ov_char8 * str);

ov_uword OVLexicon_GetNActive(OVLexicon * uk);

OVstatus OVLexicon_Pack(OVLexicon * uk);

#if 0
ov_word OVLexicon_GetCStringHash(ov_char8 * str);
#endif

#endif
