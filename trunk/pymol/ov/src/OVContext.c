
#include "OVContext.h"

OVContext *OVContext_New(void)
{
  OVContext *I = NULL;
  OVHeap *heap = OVHeap_New();
  if(heap) {
    I = OVHeap_ALLOC(heap, OVContext);
    if(I) {
      I->heap = heap;
    } else
      OVHeap_Del(heap);
  }
  return I;
}

void OVContext_Del(OVContext * I)
{
  if(I) {
    OVHeap *heap = I->heap;
    OVHeap_Free(heap, I);
    OVHeap_Del(heap);
    if(0) {
      heap = (OVHeap *) 0;      /* suppress warning message */
    }
  }
}
