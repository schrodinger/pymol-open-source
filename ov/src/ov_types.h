#ifndef _H_ov_types
#define _H_ov_types


/* Primitive types */

typedef signed int           ov_boolean;

typedef char                 ov_char8;
typedef unsigned char        ov_uchar8;
typedef signed short int     ov_int16;
typedef unsigned short int   ov_uint16;
typedef signed int           ov_int32;
typedef unsigned int         ov_uint32;

#ifdef WIN32
typedef __int64              ov_int64;
typedef unsigned __int64     ov_uint64;
#else
typedef signed long long     ov_int64;
typedef unsigned long long   ov_uint64;

#endif
typedef float                ov_float32;
typedef double               ov_float64;   

/* default: 32 bit words */

typedef ov_int32             ov_word;
typedef ov_uint32            ov_uword;
typedef ov_uint32            ov_size;

/* machine width pointer */

typedef void                *ov_pointer;

/* ASSUMPTION: ov_pointer must be equal or larger than ov_word */

#endif
