#ifndef _H_OVRandom
#define _H_OVRandom

#include "OVHeap.h"

typedef struct _OVRandom OVRandom;

OVRandom *OVRandom_NewBySeed(OVHeap *heap, ov_uint32 seed);
OVRandom *OVRandom_NewByArray(OVHeap *heap, ov_uint32 init_key[],ov_int32 key_length);

void OVRandom_Del(OVRandom *I);

/* generates a random number on [0,0xffffffff]-interval */
ov_uint32 OVRandom_Get_int32(OVRandom *I);

/* generates a random number on [0,0x7fffffff]-interval */
ov_int32 OVRandom_Get_int31(OVRandom *I);

/* generates a random number on [0,1]-real-interval */
ov_float64 OVRandom_Get_float64_inc1(OVRandom *I);

/* generates a random number on [0,1)-real-interval */
ov_float64 OVRandom_Get_float64_exc1(OVRandom *I);

/* generates a random number on (0,1)-real-interval */
ov_float64 OVRandom_Get_float64_exc01(OVRandom *I);

#endif
