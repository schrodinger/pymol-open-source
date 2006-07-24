#ifndef _H_OVOneToOne
#define _H_OVOneToOne

/* OneToOne: A bidirectional associative array for ov_words
 * 
 * associates a unique ov_word to a unique ov_word in both forward and
 * reverse directions.
 *
*/

#include "OVHeap.h"
#include "OVreturns.h"

typedef struct _OVOneToOne OVOneToOne;

OVOneToOne *OVOneToOne_New(OVHeap *heap);
OVstatus OVOneToOne_Init(OVOneToOne *o2o,OVHeap *heap);
void OVOneToOne_Purge(OVOneToOne *o2o);
void OVOneToOne_Del(OVOneToOne *o2o);
void OVOneToOne_Reset(OVOneToOne *up);

OVreturn_word OVOneToOne_IterateForward(OVOneToOne *o2o,ov_word *hidden);

OVreturn_word OVOneToOne_GetForward(OVOneToOne *o2o,ov_word forward_value);
OVreturn_word OVOneToOne_GetReverse(OVOneToOne *o2o,ov_word reverse_value);
OVstatus OVOneToOne_Set(OVOneToOne *o2o, ov_word forward_value, ov_word reverse_value);

OVstatus OVOneToOne_Pack(OVOneToOne *o2o);
OVreturn_size OVOneToOne_GetSize(OVOneToOne *o2o);
OVstatus OVOneToOne_DelReverse(OVOneToOne *o2o,ov_word reverse_value);
OVstatus OVOneToOne_DelForward(OVOneToOne *o2o,ov_word forward_value);
void OVOneToOne_Stats(OVOneToOne *o2o);
void OVOneToOne_Dump(OVOneToOne *o2o);

#define OVOneToOne_DEL_AUTO_NULL(I) { if(I) { OVOneToOne_Del(I); I=OV_NULL; }}

#endif
