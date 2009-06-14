#ifndef _H_OVOneToAny
#define _H_OVOneToAny


/* OneToAny: A bidirectional associative array for ov_words
 * 
 * associates a unique ov_word to a unique ov_word in both forward and
 * reverse directions.
 *
*/

#include "OVHeap.h"
#include "OVreturns.h"

typedef struct _OVOneToAny OVOneToAny;

OVOneToAny *OVOneToAny_New(OVHeap * heap);
OVstatus OVOneToAny_Init(OVOneToAny * o2o, OVHeap * heap);
void OVOneToAny_Purge(OVOneToAny * o2o);
void OVOneToAny_Del(OVOneToAny * o2o);
void OVOneToAny_Reset(OVOneToAny * up);

OVreturn_word OVOneToAny_GetKey(OVOneToAny * o2o, ov_word forward_value);
OVstatus OVOneToAny_SetKey(OVOneToAny * o2o, ov_word forward_value,
                           ov_word reverse_value);

OVstatus OVOneToAny_Pack(OVOneToAny * o2o);
OVreturn_size OVOneToAny_GetSize(OVOneToAny * o2o);
OVstatus OVOneToAny_DelKey(OVOneToAny * o2o, ov_word forward_value);
void OVOneToAny_Stats(OVOneToAny * o2o);
void OVOneToAny_Dump(OVOneToAny * o2o);

#define OVOneToAny_DEL_AUTO_NULL(I) { if(I) { OVOneToAny_Del(I); I=OV_NULL; }}

#endif
