

/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/

#ifndef _H_Tracker
#define _H_Tracker

#include "PyMOLGlobals.h"


/* Tracker is a generic data structure for tracking lists of
   candidates in a manner that:

   (1) is robust with respect to creation and deletion of candidates and lists
   (2) enables efficient iteration over the lists to which a candidate belongs
   (3) enables efficient iteration over the candidates to which a list belongs
   (4) is double-linked so that insertion and deletion take constant time
   (5) keeps running counts so that length determination takes constant time

*/

typedef struct _CTracker CTracker;

typedef void *TrackerRef;

CTracker *TrackerNew(PyMOLGlobals * G);
void TrackerFree(CTracker * I);


/* creating and deleting candidates, lists, and iterators */

int TrackerNewCand(CTracker * I, TrackerRef * ref);
int TrackerDelCand(CTracker * I, int cand_id);

int TrackerNewList(CTracker * I, TrackerRef * ref);
int TrackerDelList(CTracker * I, int list_id);
int TrackerNewListCopy(CTracker * I, int list_id, TrackerRef * ref);

int TrackerNewIter(CTracker * I, int cand_id, int list_id);
int TrackerDelIter(CTracker * I, int iter_id);


/* creating and destroying links */

int TrackerLink(CTracker * I, int cand_id, int list_id, int priority);
int TrackerUnlink(CTracker * I, int cand_id, int list_id);


/* querying */

int TrackerGetNList(CTracker * I);
int TrackerGetNCand(CTracker * I);
int TrackerGetNLink(CTracker * I);
int TrackerGetNIter(CTracker * I);
int TrackerGetNListForCand(CTracker * I, int cand_id);
int TrackerGetNCandForList(CTracker * I, int list_id);

int TrackerGetCandRef(CTracker * I, int cand_id, TrackerRef ** ref_ret);


/* iterating */

int TrackerIterNextCandInList(CTracker * I, int iter_id, TrackerRef ** ref_ret);
int TrackerIterNextListInCand(CTracker * I, int iter_id, TrackerRef ** ref_ret);

#ifdef TRACKER_UNIT_TEST
int TrackerUnitTest(PyMOLGlobals * G);
#endif

#endif
