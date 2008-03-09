/* 
 * COPYRIGHT NOTICE: This file contains original source code from the
 * Jenarix (TM) Library, Copyright (C) 2007-8 by Warren L. Delano of
 * DeLano Scientific LLC, Palo Alto, California, United States.
 * Please see the accompanying LICENSE file for further information.
 * All rights not explicitly granted in that LICENSE file are
 * reserved.  It is unlawful to modify or remove this notice.
 * TRADEMARK NOTICE: Jenarix is a Trademark of DeLano Scientific LLC.
*/
#ifndef _H_ov_types
#define _H_ov_types

#ifdef __cplusplus
extern "C" {
#endif

  /* automatically detect 64-bit machines */

#ifndef OV_32_BIT
#if (((size_t)-1) != 0xFFFFFFFF)
#ifndef OV_64_BIT
#define OV_64_BIT
#endif
#endif
#endif

  /* Jenarix native structs */

  typedef struct ov__object ov_object;
  typedef struct ov__gc_head ov_gc_head;
  typedef struct ov__cell ov_cell;
  typedef struct ov__block ov_block;
  typedef struct ov__opaque_head ov_opaque_head;
  typedef struct ov__triplet ov_triplet;
  typedef struct ov__node ov_node; 
  
  /* Jenarix prototypes for native functions and object dispatch
   * NOTE: parameter "void" below is a workaround -- you must cast as
   * appropriate to enable parameter type checking */
  
  typedef ov_object (*ov__fn)(void);
  typedef ov_object (*ov__dispatch_fn)(void);

  /* owned primitive types */
  
  typedef               char ov_char8;
  typedef unsigned      char ov_uchar8;
  typedef               char ov_int8;
  typedef unsigned      char ov_uint8;
  typedef          short int ov_int16;
  typedef unsigned short int ov_uint16;
  typedef                int ov_int32;
  typedef unsigned       int ov_uint32;
  typedef              float ov_float32;
  typedef             double ov_float64;
#ifdef WIN32
  typedef            __int64 ov_int64;
  typedef unsigned   __int64 ov_uint64;
#else
  typedef          long long ov_int64;
  typedef unsigned long long ov_uint64;
#endif
  typedef             size_t ov_size; 
  typedef          ptrdiff_t ov_diff;
  
#ifdef OV_64_BIT

  /* 64 bit mode: Jenarix objects are 96 bits wide and can hold a
     double-precision floating point number without owning any
     external resources */

  typedef ov_float64 ov_float;

#define OV_FLOAT_ZERO 0.0
#define OV_UWORD_MAX 0xFFFFFFFFFFFFFFFF
  
#else

  /* 32 bit mode: Jenarix objects are 64 bits wide and can only hold a
     single-precision floating point number */

  typedef ov_float32 ov_float;
  
#define OV_FLOAT_ZERO 0.0F
#define OV_UWORD_MAX 0xFFFFFFFF
  
#endif

  /* machine word types */

  typedef ov_size    ov_uword;
  typedef ov_diff    ov_word;

  /* additional derived / convenience types */

  typedef ov_uword   ov_boolean;
  
  typedef ov_uint32  ov_meta; /* object "meta" bits */
  typedef ov_int32   ov_status;
  
  typedef ov_word    ov_int; /* C-like synonyms */
  typedef ov_uword   ov_uint;
  typedef ov_char8   ov_char;
  typedef ov_uchar8  ov_uchar;

  typedef ov_word    ov_selector; /* Jenarix method selector */
  
  /* generic untyped free pointer/resource method */
  
  typedef void ov_void_free(void *ptr);
  
  /* Jenarix object struct/union */
  
  struct ov__object {
    ov_uint32 meta;
    
    union {
      /* machine word */
      
      ov_word word;
      ov_uword uword;
      
      /* method selector */
      ov_selector selector;
      
      /* functions */
      
      ov__fn fn; /* cast as ov_fn for parameter checking */
      ov__dispatch_fn dispatch_fn; /* cast as ov_dispatch_fn for
                                      parameter checking */
      
      /* C string constants (i.e. for hardcoded use with method import) */
      
      const char *str_const;
      
      /* cells & blocks, requiring external storage */
      
      ov_word *ref_cnt; 
      ov_gc_head *gc; /* garbage collection and method dispatch */
      ov_block *block;
      ov_cell *cell;
      ov_opaque_head *opaque;
      ov_triplet *triplet;
      
      /* various in-place particle payloads */
      
      ov_word boolean;
      ov_char8 char8;
      ov_uchar8 uchar8;      
      ov_int8 int8;
      ov_uint8 uint8;
      ov_int16 int16;
      ov_uint16 uint16;
      ov_int32 int32;
      ov_uint32 uint32;
      ov_float32 float32;
      
#ifdef OV_64_BIT
      ov_float64 float64;
      ov_int64 int64;
      ov_uint64 uint64;
#endif
      
      /* data pointer types */
      
      union {
        ov_word *word;
        ov_uword *uword;
        ov_word *boolean;
        ov_char8 *char8;
        ov_uchar8 *uchar8;      
        ov_int8 *int8;
        ov_uint8 *uint8;
        ov_int16 *int16;
        ov_uint16 *uint16;
        ov_int32 *int32;
        ov_uint32 *uint32;
        ov_float32 *float32;
        ov_float64 *float64;
        ov_int64 *int64;
        ov_uint64 *uint64;
      } ptr;

    } data;
  };
  
  /* unfortunately, due to circularity of definition, we can't get
   * compile-time parameter checking with :
   *    ov_object.fn(...)  
   * so must use: 
   *    ((ov_fn)ov_object.fn)(...) instead */
  
  typedef ov_object (*ov_fn)(ov_object, ov_object);
  typedef ov_object (*ov_dispatch_fn)(ov_object, ov_object, ov_object);
  
#ifdef __cplusplus
}
#endif

#endif
