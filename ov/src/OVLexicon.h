
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

OVreturn_word OVLexicon_GetFromCString(OVLexicon * uk, const ov_char8 * str);

OVstatus OVLexicon_IncRef(OVLexicon * uk, ov_word id);
OVstatus OVLexicon_DecRef(OVLexicon * uk, ov_word id);

ov_char8 *OVLexicon_FetchCString(OVLexicon * uk, ov_word id);
OVreturn_word OVLexicon_BorrowFromCString(OVLexicon * uk, const ov_char8 * str);

ov_uword OVLexicon_GetNActive(OVLexicon * uk);

OVstatus OVLexicon_Pack(OVLexicon * uk);

#if 0
ov_word OVLexicon_GetCStringHash(ov_char8 * str);
#endif

static const char EMPTY_CSTR[1] = "";

#define LexStr(G, i) (i ? OVLexicon_FetchCString(G->Lexicon, i) : EMPTY_CSTR)
#define LexIdx(G, s) ((s && s[0]) ? OVLexicon_GetFromCString(G->Lexicon, s).word : 0)
#define LexDec(G, i) OVLexicon_DecRef(G->Lexicon, i)
#define LexInc(G, i) OVLexicon_IncRef(G->Lexicon, i)

#endif
