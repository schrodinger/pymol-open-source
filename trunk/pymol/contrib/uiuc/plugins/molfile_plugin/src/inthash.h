/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2016 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: inthash.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.4 $      $Date: 2016/11/28 05:01:54 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   A simple hash table implementation for ints, contributed by John Stone,
 *   derived from his ray tracer code.
 * NOTES: - this can only used for _positive_ data values (HASH_FAIL is -1)
 *        - this code is slightly modified from the version in VMD
 *          so that both, the string hash and the int hash can be used.
 ***************************************************************************/
#ifndef INTHASH_H
#define INTHASH_H

/** hash table top level data structure */
typedef struct inthash_t {
  struct inthash_node_t **bucket;        /* array of hash nodes */
  int size;                           /* size of the array */
  int entries;                        /* number of entries in table */
  int downshift;                      /* shift cound, used in hash function */
  int mask;                           /* used to select bits for hashing */
} inthash_t;

#define HASH_FAIL -1

#if defined(VMDPLUGIN_STATIC)
#define VMDEXTERNSTATIC static
#include "inthash.c"
#else

#define VMDEXTERNSTATIC 

#ifdef __cplusplus
extern "C" {
#endif

/** initialize hash table for first use */
void inthash_init(inthash_t *, int);

/** lookup a string key in the hash table returning its integer key */
int inthash_lookup(const inthash_t *, int);

/** insert a string into the hash table, along with an integer key */
int inthash_insert(inthash_t *, int, int);

/** delete an string from the hash table, given its string name */
int inthash_delete(inthash_t *, int);

/** return the number of entries in the hash table */
int inthash_entries(inthash_t *);

/** destroy the hash table completely, deallocate memory */
void inthash_destroy(inthash_t *);

/** print hash table vital stats */
char *inthash_stats(inthash_t *);

#ifdef __cplusplus
}
#endif

#endif

#endif
