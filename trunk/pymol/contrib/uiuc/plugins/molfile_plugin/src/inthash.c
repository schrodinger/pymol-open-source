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
 *      $RCSfile: inthash.c,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.3 $      $Date: 2016/11/28 05:01:54 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   A simple hash table implementation for ints, contributed by John Stone,
 *   derived from his ray tracer code.
 * NOTES: - this can only used for _positive_ data values (HASH_FAIL is -1)
 *        - this code is slightly modified from the version in VMD
 *          so that both, the string hash and the int hash can be used.
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inthash.h"

#define HASH_LIMIT 0.5

/** hash table node data structure */
typedef struct inthash_node_t {
  int data;                           /* data in hash node */
  int key;                            /* key for hash lookup */
  struct inthash_node_t *next;        /* next node in hash chain */
} inthash_node_t;

/*
 *  inthash() - Hash function returns a hash number for a given key.
 *
 *  tptr: Pointer to a hash table
 *  key: The key to create a hash number for
 */
static int inthash(const inthash_t *tptr, int key) {
  int hashvalue;

  hashvalue = (((key*1103515249)>>tptr->downshift) & tptr->mask);
  if (hashvalue < 0) {
    hashvalue = 0;
  }    

  return hashvalue;
}

/*
 *  inthash_init() - Initialize a new hash table.
 *
 *  tptr: Pointer to the hash table to initialize
 *  buckets: The number of initial buckets to create
 */
VMDEXTERNSTATIC void inthash_init(inthash_t *tptr, int buckets) {

  /* make sure we allocate something */
  if (buckets==0)
    buckets=16;

  /* initialize the table */
  tptr->entries=0;
  tptr->size=2;
  tptr->mask=1;
  tptr->downshift=29;

  /* ensure buckets is a power of 2 */
  while (tptr->size<buckets) {
    tptr->size<<=1;
    tptr->mask=(tptr->mask<<1)+1;
    tptr->downshift--;
  } /* while */

  /* allocate memory for table */
  tptr->bucket=(inthash_node_t **) calloc(tptr->size, sizeof(inthash_node_t *));

  return;
}

/*
 *  rebuild_table_int() - Create new hash table when old one fills up.
 *
 *  tptr: Pointer to a hash table
 */
static void rebuild_table_int(inthash_t *tptr) {
  inthash_node_t **old_bucket, *old_hash, *tmp;
  int old_size, h, i;

  old_bucket=tptr->bucket;
  old_size=tptr->size;

  /* create a new table and rehash old buckets */
  inthash_init(tptr, old_size<<1);
  for (i=0; i<old_size; i++) {
    old_hash=old_bucket[i];
    while(old_hash) {
      tmp=old_hash;
      old_hash=old_hash->next;
      h=inthash(tptr, tmp->key);
      tmp->next=tptr->bucket[h];
      tptr->bucket[h]=tmp;
      tptr->entries++;
    } /* while */
  } /* for */

  /* free memory used by old table */
  free(old_bucket);

  return;
}

/*
 *  inthash_lookup() - Lookup an entry in the hash table and return a pointer to
 *    it or HASH_FAIL if it wasn't found.
 *
 *  tptr: Pointer to the hash table
 *  key: The key to lookup
 */
VMDEXTERNSTATIC int inthash_lookup(const inthash_t *tptr, int key) {
  int h;
  inthash_node_t *node;


  /* find the entry in the hash table */
  h=inthash(tptr, key);
  for (node=tptr->bucket[h]; node!=NULL; node=node->next) {
    if (node->key == key)
      break;
  }

  /* return the entry if it exists, or HASH_FAIL */
  return(node ? node->data : HASH_FAIL);
}

/*
 *  inthash_insert() - Insert an entry into the hash table.  If the entry already
 *  exists return a pointer to it, otherwise return HASH_FAIL.
 *
 *  tptr: A pointer to the hash table
 *  key: The key to insert into the hash table
 *  data: A pointer to the data to insert into the hash table
 */
VMDEXTERNSTATIC int inthash_insert(inthash_t *tptr, int key, int data) {
  int tmp;
  inthash_node_t *node;
  int h;

  /* check to see if the entry exists */
  if ((tmp=inthash_lookup(tptr, key)) != HASH_FAIL)
    return(tmp);

  /* expand the table if needed */
  while (tptr->entries>=HASH_LIMIT*tptr->size)
    rebuild_table_int(tptr);

  /* insert the new entry */
  h=inthash(tptr, key);
  node=(struct inthash_node_t *) malloc(sizeof(inthash_node_t));
  node->data=data;
  node->key=key;
  node->next=tptr->bucket[h];
  tptr->bucket[h]=node;
  tptr->entries++;

  return HASH_FAIL;
}

/*
 *  inthash_delete() - Remove an entry from a hash table and return a pointer
 *  to its data or HASH_FAIL if it wasn't found.
 *
 *  tptr: A pointer to the hash table
 *  key: The key to remove from the hash table
 */
VMDEXTERNSTATIC int inthash_delete(inthash_t *tptr, int key) {
  inthash_node_t *node, *last;
  int data;
  int h;

  /* find the node to remove */
  h=inthash(tptr, key);
  for (node=tptr->bucket[h]; node; node=node->next) {
    if (node->key == key)
      break;
  }

  /* Didn't find anything, return HASH_FAIL */
  if (node==NULL)
    return HASH_FAIL;

  /* if node is at head of bucket, we have it easy */
  if (node==tptr->bucket[h])
    tptr->bucket[h]=node->next;
  else {
    /* find the node before the node we want to remove */
    for (last=tptr->bucket[h]; last && last->next; last=last->next) {
      if (last->next==node)
        break;
    }
    last->next=node->next;
  }

  /* free memory and return the data */
  data=node->data;
  free(node);

  return(data);
}


/*
 * inthash_entries() - return the number of hash table entries.
 *
 */
VMDEXTERNSTATIC int inthash_entries(inthash_t *tptr) {
  return tptr->entries;
}


/*
 * inthash_destroy() - Delete the entire table, and all remaining entries.
 * 
 */
VMDEXTERNSTATIC void inthash_destroy(inthash_t *tptr) {
  inthash_node_t *node, *last;
  int i;

  for (i=0; i<tptr->size; i++) {
    node = tptr->bucket[i];
    while (node != NULL) { 
      last = node;   
      node = node->next;
      free(last);
    }
  }     

  /* free the entire array of buckets */
  if (tptr->bucket != NULL) {
    free(tptr->bucket);
    memset(tptr, 0, sizeof(inthash_t));
  }
}

/*
 *  alos_int() - Find the average length of search.
 *
 *  tptr: Pointer to a hash table
 */
static float alos_int(inthash_t *tptr) {
  int i,j;
  float alos_int=0;
  inthash_node_t *node;


  for (i=0; i<tptr->size; i++) {
    for (node=tptr->bucket[i], j=0; node!=NULL; node=node->next, j++);
    if (j)
      alos_int+=((j*(j+1))>>1);
  } /* for */

  return(tptr->entries ? alos_int/tptr->entries : 0);
}


/*
 *  inthash_stats() - Return a string with stats about a hash table.
 *
 *  tptr: A pointer to the hash table
 */
VMDEXTERNSTATIC char * inthash_stats(inthash_t *tptr) {
  static char buf[1024];

  sprintf(buf, "%u slots, %u entries, and %1.2f ALOS",
    (int)tptr->size, (int)tptr->entries, alos_int(tptr));

  return(buf);
}




