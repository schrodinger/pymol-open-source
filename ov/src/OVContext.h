#ifndef _H_OVContext
#define _H_OVContext

#ifndef OVCONTEXT_DEFINED
typedef struct _OVContext OVContext;
#define OVCONTEXT_DEFINED
#endif

OVContext *OVContext_New(void);
void OVContext_Del(OVContext *I);

#include "OVHeap.h"

/* should only be accessed by special methods */

struct _OVContext {
  OVHeap *heap;
};

#endif
