

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

#include"os_predef.h"
#include"os_limits.h"
#include"os_std.h"

#include <unordered_map>
#include"OVContext.h"
#include"Base.h"
#include"MemoryDebug.h"
#include"Tracker.h"
#include"Util.h"

#define CAND_INFO 1
#define LIST_INFO 2
#define ITER_INFO 3


/* double-linked throughout */

struct TrackerInfo {
  int id;
  int type;
  int first, last;
  TrackerRef *ref;
  int length;
  int next, prev;
};

struct TrackerMember {
  int cand_id, cand_index;
  int cand_next, cand_prev;
  int list_id, list_index;
  int list_next, list_prev;
  int hash_next, hash_prev;
  int priority;
};

struct TrackerIter {
  int id;
  int member;
};

struct CTracker {
  int next_id = 1;
  int next_free_info = 0;
  int next_free_member = 0;
  int n_cand{}, n_list{};
  int n_info = 0;
  int n_member = 0;
  int n_link{};
  int n_iter{};
  int cand_start{};
  int list_start{};
  int iter_start{};
  std::vector<TrackerInfo> info;
  std::unordered_map<int, int> id2info;
  std::unordered_map<int, int> hash2member;
  std::vector<TrackerMember> member;
};

#define TRACKER_HASH_KEY(a,b) (a^b)

CTracker *TrackerNew(PyMOLGlobals * G)
{
  auto I = new CTracker();

  I->info.push_back(TrackerInfo{});
  I->member.push_back(TrackerMember{});
  return I;
}

static int GetNewInfo(CTracker * I)
{
  int result = 0;
  if(!I->next_free_info) {
    I->n_info++;
    result = I->n_info;
    I->info.push_back(TrackerInfo{});
  } else {
    result = I->next_free_info;
    I->next_free_info = I->info[result].next;
    I->info[result] = TrackerInfo{};
  }
  return result;
}

static int GetNewMember(CTracker * I)
{
  int result = 0;
  if(!(result = I->next_free_member)) {
    I->n_member++;
    result = I->n_member;
    I->member.push_back(TrackerMember{});
  } else {
    I->next_free_member = I->member[result].hash_next;
    I->member[result] = TrackerMember{};
  }
  I->n_link++;
  return result;
}

static void ReleaseInfo(CTracker * I, int index)
{
  I->info[index].next = I->next_free_info;
  I->next_free_info = index;
}

static void ReleaseMember(CTracker * I, int index)
{
  I->member[index].hash_next = I->next_free_member;
  I->next_free_member = index;
  I->n_link--;
}

static int GetUniqueValidID(CTracker * I)
{
  int result = I->next_id;
  for (auto id2info : I->id2info)
  {
    auto result = id2info.second;
    result = (result + 1) & 0x7FFFFFFF;
    if(!result)
      result = 1;
  }
  if(!(I->next_id = (result + 1) & 0x7FFFFFFF))
    I->next_id = 1;

  return result;
}

static void ProtectIterators(CTracker * I, int member_index)
{
  auto I_info = I->info.data();
  int iter_index;
  if((iter_index = I->iter_start) && member_index) {
    while(iter_index) {
      TrackerInfo *info = I_info + iter_index;
      if(info->first == member_index) {
        auto member = &I->member[member_index];
        switch (info->length) {
        case CAND_INFO:
          info->first = member->cand_next;
          break;
        case LIST_INFO:
          info->first = member->list_next;
          break;
        default:
          info->first = 0;
          break;
        }
      } else if(info->last == member_index) {
        auto member = &I->member[member_index];
        switch (info->length) {
        case CAND_INFO:
          info->last = member->cand_prev;
          break;
        case LIST_INFO:
          info->last = member->list_prev;
          break;
        default:
          info->last = 0;
          break;
        }
      }
      iter_index = info->next;
    }
  }
}

int TrackerNewCand(CTracker * I, TrackerRef * ref)
{
  int result = 0;
  int index = GetNewInfo(I);
  int id;
  auto I_info = I->info.data();
  if(index) {
    TrackerInfo *info = I_info + index;
    info->ref = ref;

    info->next = I->cand_start;
    if(info->next)
      I_info[info->next].prev = index;

    I->cand_start = index;
    id = GetUniqueValidID(I);
    I->id2info[id] = index;
    info->id = (result = id);
    info->type = CAND_INFO;
    I->n_cand++;
  }
  return result;
}

int TrackerNewList(CTracker * I, TrackerRef * ref)
{
  int result = 0;
  int index = GetNewInfo(I);
  int id;
  auto I_info = I->info.data();
  if(index) {
    TrackerInfo *info = I_info + index;
    info->ref = ref;
    info->next = I->list_start;
    if(info->next)
      I_info[info->next].prev = index;
    I->list_start = index;
    id = GetUniqueValidID(I);
    I->id2info[id] = index;
    info->id = (result = id);
    info->type = LIST_INFO;
    I->n_list++;
  }
  return result;
}

int TrackerNewListCopy(CTracker * I, int list_id, TrackerRef * ref)
{
  int new_list_id = TrackerNewList(I, ref);
  int iter_id = TrackerNewIter(I, 0, list_id);
  if(iter_id) {
    int cand_id;
    while((cand_id = TrackerIterNextCandInList(I, iter_id, NULL))) {
      TrackerLink(I, cand_id, new_list_id, 1);
    }
    TrackerDelIter(I, iter_id);
  }
  return new_list_id;
}

int TrackerNewIter(CTracker * I, int cand_id, int list_id)
{
  int result = 0;
  if((cand_id >= 0) || (list_id >= 0)) {
    int index = GetNewInfo(I);
    int id;
    auto I_info = I->info.data();
    if(index) {
      TrackerInfo *info = I_info + index;
      info->next = I->iter_start;
      if(info->next)
        I_info[info->next].prev = index;
      I->iter_start = index;
      id = GetUniqueValidID(I);
      I->id2info[id] = index;
      info->id = (result = id);
      info->type = ITER_INFO;
      I->n_iter++;
      if(cand_id && list_id) {        /* seeking a specific member */
        int hash_key = TRACKER_HASH_KEY(cand_id, list_id);
        auto hashStartIt = I->hash2member.find(hash_key);
        if (hashStartIt != I->hash2member.end()) {
          int member_index = hashStartIt->second;
          TrackerMember *I_member = I->member.data();
          TrackerMember *member;
          while(member_index) {
            member = I_member + member_index;
            if((member->cand_id == cand_id) && (member->list_id == list_id)) {
              info->first = member_index;
              break;
            }
            member_index = member->hash_next;
          }
        }
      } else if(list_id) {    /* for iterating over cands in a list */
        auto listIndexIt = I->id2info.find(list_id);
        if(listIndexIt != I->id2info.end()) {
          TrackerInfo *list_info = I_info + listIndexIt->second;
          info->first = list_info->first;
        }
      } else if(cand_id) {    /* for iterating over lists in a cand */
        auto candIndexIt = I->id2info.find(cand_id);
        if (candIndexIt != I->id2info.end()) {
          TrackerInfo *cand_info = I_info + candIndexIt->second;
          info->first = cand_info->first;
        }
      }
    }
  }
  return result;
}

int TrackerDelCand(CTracker * I, int cand_id)
{
  int result = false;
  if(cand_id >= 0) {
    auto candIndexIt = I->id2info.find(cand_id);
    auto I_info = I->info.data();

    if (candIndexIt != I->id2info.end()) {
      auto cand_index = candIndexIt->second;
      TrackerInfo *cand_info = I_info + cand_index;

      if(cand_info->type == CAND_INFO) {
        result = true;

        {                       /* first release all the members */
          int iter_start = I->iter_start;
          auto I_member = I->member.data();
          int member_index = cand_info->first;

          while(member_index) {
            TrackerMember *member = I_member + member_index;
            TrackerInfo *list_info = I_info + member->list_index;
            int cand_id = member->cand_id, list_id = member->list_id;

            /* if any iterators exist, then make sure they don't point at
               this member */
            if(iter_start) {
              ProtectIterators(I, member_index);
            }

            {
              /* extract from hash chain */
              int prev = member->hash_prev;
              int next = member->hash_next;
              int hash_key = TRACKER_HASH_KEY(cand_id, list_id);

              if(prev) {
                I_member[prev].hash_next = next;
              } else {
                I->hash2member.erase(hash_key);
                if(member->hash_next) {
                  I->hash2member[hash_key] = member->hash_next;
                  /* ASSUMING SUCCESS -- NOT TESTING FOR ERROR */
                }
              }
              if(next) {
                I_member[next].hash_prev = prev;
              }
            }

            {
              /* extract from list chain */
              int prev = member->list_prev;
              int next = member->list_next;

              if(prev) {
                I_member[prev].list_next = next;
              } else {
                list_info->first = next;
              }

              if(next) {
                I_member[next].list_prev = prev;
              } else {
                list_info->last = prev;
              }

              /* shorten length of the list */
              list_info->length--;
            }

            {                   /* continue along the chain */
              int member_index_copy = member_index;
              member_index = member->cand_next;
              ReleaseMember(I, member_index_copy);
            }
          }
        }

        /* delete the cand id */
        I->id2info.erase(cand_id);

        /* remove the cand from the cand chain */
        {
          int prev = cand_info->prev;
          int next = cand_info->next;

          if(prev) {
            I->info[prev].next = next;
          } else {
            I->cand_start = next;
          }

          if(next) {
            I->info[next].prev = prev;
          }

          /* shorten length of the list */
          I->n_cand--;
        }

        /* and release the info record */

        ReleaseInfo(I, cand_index);
      }
    }
  }
  return result;

}

int TrackerIterNextCandInList(CTracker * I, int iter_id, TrackerRef ** ref_ret)
{
  /* returns the next cand_id in the list */
  int result = 0;
  if(iter_id >= 0) {
    auto iterIndexIt = I->id2info.find(iter_id);
    auto I_info = I->info.data();
    if (iterIndexIt != I->id2info.end()) {
      TrackerInfo *iter_info = I_info + iterIndexIt->second;
      int member_index;
      if((member_index = iter_info->first)) {
        auto member = &I->member[member_index];
        result = member->cand_id;
        if(ref_ret) {
          TrackerInfo *cand_info = I_info + member->cand_index;
          *ref_ret = cand_info->ref;
        }
        iter_info->last = iter_info->first;
        iter_info->first = member->list_next;
      } else if((member_index = iter_info->last)) {     /* first is zero, so try last */
        auto member = &I->member[member_index];
        if(member->list_next) {
          member = &I->member[member->list_next];
          result = member->cand_id;
          if(ref_ret) {
            TrackerInfo *cand_info = I_info + member->cand_index;
            *ref_ret = cand_info->ref;
          }
          iter_info->last = iter_info->first;
          iter_info->first = member->list_next;
        }
      }
      iter_info->length = LIST_INFO;
    }
  }
  return result;
}

int TrackerIterNextListInCand(CTracker * I, int iter_id, TrackerRef ** ref_ret)
{
  /* returns the next cand_id in the list */
  int result = 0;
  if(iter_id >= 0) {
    auto iterIndexIt = I->id2info.find(iter_id);
    auto I_info = I->info.data();
    if (iterIndexIt != I->id2info.end()) {
      TrackerInfo *iter_info = I_info + iterIndexIt->second;
      int member_index;
      if((member_index = iter_info->first)) {
        auto member = &I->member[member_index];
        result = member->list_id;
        if(ref_ret) {
          TrackerInfo *list_info = I_info + member->list_index;
          *ref_ret = list_info->ref;
        }
        iter_info->last = member_index;
        iter_info->first = member->cand_next;
      } else if((member_index = iter_info->last)) {     /* first is zero, so try last */
        auto member = &I->member[member_index];
        if(member->cand_next) {
          member = &I->member[member->cand_next];
          result = member->list_id;
          if(ref_ret) {
            TrackerInfo *list_info = I_info + member->list_index;
            *ref_ret = list_info->ref;
          }
          iter_info->last = member_index;
          iter_info->first = member->cand_next;
        }
      }
      iter_info->length = CAND_INFO;
    }
  }
  return result;
}

int TrackerDelIter(CTracker * I, int iter_id)
{
  int result = false;
  if(iter_id >= 0) {
    auto iterIndexIt = I->id2info.find(iter_id);
    auto I_info = I->info.data();
    if (iterIndexIt != I->id2info.end()) {
      auto iter_index = iterIndexIt->second;
      /* remove the iter from the iter chain */
      TrackerInfo *iter_info = I_info + iter_index;

      int prev = iter_info->prev;
      int next = iter_info->next;

      if(prev) {
        I->info[prev].next = next;
      } else {
        I->iter_start = next;
      }

      if(next) {
        I->info[next].prev = prev;
      }

      /* delete the iter id */
      I->id2info.erase(iter_id);

      /* shorten length of the list */
      I->n_iter--;

      result = true;

      ReleaseInfo(I, iter_index);
    }
  }
  return result;
}

int TrackerGetNList(CTracker * I)
{
  return I->n_list;
}

int TrackerGetNCand(CTracker * I)
{
  return I->n_cand;
}

int TrackerGetNLink(CTracker * I)
{
  return I->n_link;
}

int TrackerGetNIter(CTracker * I)
{
  return I->n_iter;
}

int TrackerGetCandRef(CTracker * I, int cand_id, TrackerRef ** ref_ret)
{
  auto candIndexIt = I->id2info.find(cand_id);
  auto I_info = I->info.data();

  if (candIndexIt != I->id2info.end()) {
    TrackerInfo *info = I_info + candIndexIt->second;
    if(info->type == CAND_INFO) {
      *ref_ret = info->ref;
      return true;
    }
  }
  return false;
}

int TrackerGetNListForCand(CTracker * I, int cand_id)
{
  auto candIndexIt = I->id2info.find(cand_id);
  auto I_info = I->info.data();

  if (candIndexIt != I->id2info.end()) {
    TrackerInfo *info = I_info + candIndexIt->second;
    if(info->type == CAND_INFO)
      return info->length;
  }
  return -1;
}

int TrackerGetNCandForList(CTracker * I, int list_id)
{
  auto listIndexIt = I->id2info.find(list_id);
  auto I_info = I->info.data();

  if (listIndexIt != I->id2info.end()) {
    TrackerInfo *info = I_info + listIndexIt->second;
    if(info->type == LIST_INFO)
      return info->length;
  }
  return -1;
}

int TrackerDelList(CTracker * I, int list_id)
{
  int result = false;
  if(list_id >= 0) {
    auto listIndexIt = I->id2info.find(list_id);
    auto I_info = I->info.data();

  if (listIndexIt != I->id2info.end()) {
      auto list_index = listIndexIt->second;
      TrackerInfo *list_info = I_info + list_index;

      if(list_info->type == LIST_INFO) {

        result = true;

        {                       /* first release all the members */
          int iter_start = I->iter_start;
          auto I_member = I->member.data();
          int member_index = list_info->first;

          while(member_index) {
            TrackerMember *member = I_member + member_index;
            TrackerInfo *cand_info = I_info + member->cand_index;
            int cand_id = member->cand_id, list_id = member->list_id;

            if(iter_start) {
              ProtectIterators(I, member_index);
            }

            {
              /* extract from hash chain */
              int prev = member->hash_prev;
              int next = member->hash_next;
              int hash_key = TRACKER_HASH_KEY(cand_id, list_id);

              if(prev) {
                I_member[prev].hash_next = next;
              } else {
                I->hash2member.erase(hash_key);
                if(member->hash_next) {
                  I->hash2member[hash_key] = member->hash_next;
                  /* ASSUMING SUCCESS -- NOT TESTING FOR ERROR */
                }
              }
              if(next) {
                I_member[next].hash_prev = prev;
              }
            }

            {
              /* extract from cand chain */
              int prev = member->cand_prev;
              int next = member->cand_next;

              if(prev) {
                I_member[prev].cand_next = next;
              } else {
                cand_info->first = next;
              }

              if(next) {
                I_member[next].cand_prev = prev;
              } else {
                cand_info->last = prev;
              }
              /* shorten length of the list */
              cand_info->length--;
            }

            /* if any iterators exist, then make sure they don't point at
               this member */

            {                   /* continue along the chain */
              int member_index_copy = member_index;
              member_index = member->list_next;
              ReleaseMember(I, member_index_copy);
            }
          }
        }

        /* delete the list id */
        I->id2info.erase(list_id);

        /* remove the list from the list chain */
        {
          int prev = list_info->prev;
          int next = list_info->next;

          if(prev) {
            I->info[prev].next = next;
          } else {
            I->list_start = next;
          }

          if(next) {
            I->info[next].prev = prev;
          }

          /* shorten length of the list */
          I->n_list--;
        }

        /* and release the info record */

        ReleaseInfo(I, list_index);

      }
    }
  }
  return result;
}

int TrackerLink(CTracker * I, int cand_id, int list_id, int priority)
{
  int result = false;
  int hash_key = TRACKER_HASH_KEY(cand_id, list_id);
  int already_linked = false;
  auto hashStartIt = I->hash2member.find(hash_key);
  ov_word hash_start{};

  {
    if (hashStartIt != I->hash2member.end()) {
      hash_start = hashStartIt->second;
      int member_index = hash_start;
      auto I_member = I->member.data();
      TrackerMember *member;
      while(member_index) {
        member = I_member + member_index;
        if((member->cand_id == cand_id) && (member->list_id == list_id)) {
          already_linked = true;
          break;
        }
        member_index = member->hash_next;
      }
    } else {
      hash_start = 0;      /* will be first entry in the hash chain */
    }
  }

  if(!already_linked) {

    auto candIndexIt = I->id2info.find(cand_id);
    auto listIndexIt = I->id2info.find(list_id);

    if (candIndexIt != I->id2info.end() && listIndexIt != I->id2info.end()) {
      auto cand_index = candIndexIt->second;
      auto list_index = listIndexIt->second;
      auto I_info = I->info.data();
      TrackerInfo *cand_info = I_info + cand_index;
      TrackerInfo *list_info = I_info + list_index;

      if(!already_linked) {
        int member_index = GetNewMember(I);
        if(member_index) {

          if(!hash_start) {        /* not a */
            I->hash2member[hash_key] = member_index;
            hash_start = member_index;
          }

          if(!hash_start) {
            ReleaseMember(I, member_index);
          } else {
            /* cannot fail now */

            auto I_member = I->member.data();
            TrackerMember *member = I_member + member_index;

            result = true;

            cand_info->length++;
            list_info->length++;

            member = &I->member[member_index];

            member->priority = priority;
            member->cand_id = cand_id;
            member->cand_index = cand_index;
            member->list_id = list_id;
            member->list_index = list_index;

            /* insert into the hash chain at start in second spot */

            if(hash_start != member_index) {       /* in the second spot */
              member->hash_prev = hash_start;
              member->hash_next = I_member[hash_start].hash_next;
              I_member[hash_start].hash_next = member_index;
              if(member->hash_next) {
                I_member[member->hash_next].hash_prev = member_index;
              }
            }

            /* else, member is sole member of the hash chain, so no links are required */
            /* insert at end of cand chain */
            member->cand_prev = cand_info->last;
            cand_info->last = member_index;
            if(member->cand_prev) {
              I_member[member->cand_prev].cand_next = member_index;
            } else {
              cand_info->first = member_index;
            }

            /* insert at end of list chain */
            member->list_prev = list_info->last;
            list_info->last = member_index;
            if(member->list_prev) {
              I_member[member->list_prev].list_next = member_index;
            } else {
              list_info->first = member_index;
            }
          }
        }
      }
    }
  }
  return result;
}

int TrackerUnlink(CTracker * I, int cand_id, int list_id)
{
  int result = false;
  int hash_key = TRACKER_HASH_KEY(cand_id, list_id);
  int already_linked = false;
  auto hashStartIt = I->hash2member.find(hash_key);
  auto I_member = I->member.data();
  TrackerMember *member = NULL;
  ov_word member_index{};

  {
    if (hashStartIt != I->hash2member.end()) {
      member_index = hashStartIt->second;
      while(member_index) {
        member = I_member + member_index;
        if((member->cand_id == cand_id) && (member->list_id == list_id)) {
          already_linked = true;
          break;
        }
        member_index = member->hash_next;
      }
    }
  }

  if(already_linked) {          /* member and member_index will be valid */

    auto I_info = I->info.data();
    TrackerInfo *cand_info = I_info + member->cand_index;
    TrackerInfo *list_info = I_info + member->list_index;

    result = true;

    /* if any iterators exist, then make sure they don't point at
       this member */
    if(I->iter_start) {
      ProtectIterators(I, member_index);
    }

    {
      /* extract from hash chain */
      int prev = member->hash_prev;
      int next = member->hash_next;

      if(prev) {
        I_member[prev].hash_next = next;
      } else {
        I->hash2member.erase(hash_key);
        if(member->hash_next) {
          I->hash2member[hash_key] = member->hash_next;
          /* ASSUMING SUCCESS -- NOT TESTING FOR ERROR */
        }
      }
      if(next) {
        I_member[next].hash_prev = prev;
      }
    }

    {
      /* extract from cand chain */
      int prev = member->cand_prev;
      int next = member->cand_next;

      if(prev) {
        I_member[prev].cand_next = next;
      } else {
        cand_info->first = next;
      }

      if(next) {
        I_member[next].cand_prev = prev;
      } else {
        cand_info->last = prev;
      }
      /* shorten length of the list */
      cand_info->length--;
    }

    {
      /* extract from list chain */
      int prev = member->list_prev;
      int next = member->list_next;

      if(prev) {
        I_member[prev].list_next = next;
      } else {
        list_info->first = next;
      }

      if(next) {
        I_member[next].list_prev = prev;
      } else {
        list_info->last = prev;
      }
      /* shorten length of the list */
      list_info->length--;
    }
    /* release the member for reuse */

    ReleaseMember(I, member_index);
  }
  return result;
}

void TrackerFree(CTracker * I)
{
  DeleteP(I);
}

#ifdef TRACKER_UNIT_TEST

#include <random>

#define N_ID 100

int TrackerUnitTest(PyMOLGlobals * G)
{
  int result = true;
  int cand_id[N_ID], list_id[N_ID], iter_id[N_ID], iter_start[N_ID];
  int a;
  int tmp_int;

  std::mt19937 mt{12345678};
  std::uniform_real_distribution<float> dist{};
  CTracker *I = TrackerNew(G);

  for(a = 0; a < N_ID; a++) {
    iter_id[a] = 0;
  }
  /* first test simple new/del */

  for(a = 0; a < N_ID; a++) {
    cand_id[a] = TrackerNewCand(I, NULL);
    if(!cand_id[a])
      fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; %d==0\n", __LINE__,
              cand_id[a]);
  }

  for(a = 0; a < N_ID; a++) {
    list_id[a] = TrackerNewList(I, NULL);
    if(!list_id[a])
      fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; %d==0\n", __LINE__,
              list_id[a]);
  }

  for(a = 0; a < N_ID; a++) {
    if(cand_id[a]) {
      if(TrackerDelCand(I, cand_id[a]))
        cand_id[a] = 0;
      else
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d\n", __LINE__);
    }
    if(list_id[a]) {
      if(TrackerDelList(I, list_id[a]))
        list_id[a] = 0;
      else
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d\n", __LINE__);
    }
  }

  if((tmp_int = TrackerGetNList(I))) {
    fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; %d!=0\n", __LINE__, tmp_int);
  }

  if((tmp_int = TrackerGetNCand(I))) {
    fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; %d!=0\n", __LINE__, tmp_int);
  }

  for(a = 0; a < N_ID; a++) {
    cand_id[a] = TrackerNewCand(I, NULL);
    list_id[a] = TrackerNewList(I, NULL);
  }

  /* test simple serial linking and unlinking */

  {

    for(a = 0; a < N_ID; a++) {
      if(!TrackerLink(I, cand_id[a], list_id[a], 1)) {
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; a=%d \n", __LINE__, a);
      }
      if(!TrackerUnlink(I, cand_id[a], list_id[a])) {
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; a=%d\n", __LINE__, a);
      }
    }

    for(a = 0; a < N_ID; a++) {
      if(!TrackerLink(I, cand_id[a], list_id[a], 1)) {
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; a=%d\n", __LINE__, a);
      }
    }

    if(N_ID != TrackerGetNLink(I)) {
      fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; %d!=%d\n", __LINE__, N_ID,
              TrackerGetNLink(I));
    }

    for(a = 0; a < N_ID; a++) {
      if(!TrackerUnlink(I, cand_id[a], list_id[a])) {
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; a=%d\n", __LINE__, a);
      }
    }
  }

  /* next test random linking and unlinking, followed by systematic unlinking */

  {
    int n_link = 0;
    int list_idx, cand_idx;
    int b;

    for(a = 0; a < (N_ID * N_ID); a++) {
      cand_idx = (int) (N_ID * dist(mt));
      list_idx = (int) (N_ID * dist(mt));

      if(TrackerLink(I, cand_id[cand_idx], list_id[list_idx], 0))
        n_link++;
    }

    if(n_link != TrackerGetNLink(I)) {
      fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; %d!=%d\n", __LINE__, n_link,
              TrackerGetNLink(I));
    }

    for(a = 0; a < (N_ID * N_ID); a++) {
      cand_idx = (int) (N_ID * dist(mt));
      list_idx = (int) (N_ID * dist(mt));

      if(TrackerUnlink(I, cand_id[cand_idx], list_id[list_idx]))
        n_link--;
    }

    if(n_link != TrackerGetNLink(I)) {
      fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; %d!=%d\n", __LINE__, n_link,
              TrackerGetNLink(I));
    }

    for(a = 0; a < N_ID; a++) {
      for(b = 0; b < N_ID; b++) {
        if(TrackerUnlink(I, cand_id[a], list_id[b]))
          n_link--;
      }
    }
    if(n_link != TrackerGetNLink(I)) {
      fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; %d!=%d\n", __LINE__, n_link,
              TrackerGetNLink(I));
    }
    if(n_link != 0) {
      fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; 0!=%d\n", __LINE__, n_link);
    }
  }

  /* next test random linking and list deletion */

  {
    int n_link = 0;
    int list_idx, cand_idx;

    for(a = 0; a < (N_ID * N_ID); a++) {
      cand_idx = (int) (N_ID * dist(mt));
      list_idx = (int) (N_ID * dist(mt));

      if(TrackerLink(I, cand_id[cand_idx], list_id[list_idx], 0))
        n_link++;
    }

    if(n_link != TrackerGetNLink(I)) {
      fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; %d!=%d\n", __LINE__, n_link,
              TrackerGetNLink(I));
    }
  }

  for(a = 0; a < N_ID; a++) {
    if(cand_id[a]) {
      if(TrackerDelCand(I, cand_id[a]))
        cand_id[a] = 0;
      else
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d\n\n", __LINE__);
    }

    if(list_id[a]) {
      if(TrackerDelList(I, list_id[a]))
        list_id[a] = 0;
      else
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d\n\n", __LINE__);
    }
  }

  if(TrackerGetNLink(I)) {
    fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; 0!=%d\n", __LINE__,
            TrackerGetNLink(I));
  }

  /* make sure everyone was deleted */

  for(a = 0; a < N_ID; a++) {
    if(cand_id[a])
      fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d;  %d!=0\n", __LINE__,
              cand_id[a]);
    if(list_id[a])
      fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d;  %d!=0\n", __LINE__,
              list_id[a]);
  }

  /* next test random linking and unlinking and list deletion */

  {
    int n_link = 0;
    int list_idx, cand_idx;
    int b;

    for(a = 0; a < N_ID; a++) {

      for(b = 0; b < N_ID; b++) {
        cand_idx = (int) (N_ID * dist(mt));
        list_idx = (int) (N_ID * dist(mt));

        if(!cand_id[cand_idx]) {
          cand_id[cand_idx] = TrackerNewCand(I, NULL);
        }
        if(!list_id[list_idx]) {
          list_id[list_idx] = TrackerNewList(I, NULL);
        }

        if(cand_id[cand_idx] && list_id[list_idx]) {
          if(TrackerLink(I, cand_id[cand_idx], list_id[list_idx], 0))
            n_link++;
        }
      }

      cand_idx = (int) (N_ID * dist(mt));
      list_idx = (int) (N_ID * dist(mt));

      if(cand_id[cand_idx] && list_id[list_idx]) {
        if(TrackerUnlink(I, cand_id[cand_idx], list_id[list_idx]))
          n_link--;
      }

      cand_idx = (int) (N_ID * dist(mt));
      list_idx = (int) (N_ID * dist(mt));

      if(cand_id[cand_idx]) {
        int len = TrackerGetNListForCand(I, cand_id[cand_idx]);
        if(len >= 0) {
          if(TrackerDelCand(I, cand_id[cand_idx])) {
            cand_id[cand_idx] = 0;
            n_link -= len;
          } else
            fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d\n\n", __LINE__);
        }
      }

      if(list_id[list_idx]) {
        int len = TrackerGetNCandForList(I, list_id[list_idx]);
        if(len >= 0) {
          if(TrackerDelList(I, list_id[list_idx])) {
            list_id[list_idx] = 0;
            n_link -= len;
          } else
            fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d\n\n", __LINE__);
        }
      }

    }

    if(n_link != TrackerGetNLink(I)) {
      fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; %d!=%d\n", __LINE__, n_link,
              TrackerGetNLink(I));
    }
  }

  for(a = 0; a < N_ID; a++) {
    if(cand_id[a]) {
      if(TrackerDelCand(I, cand_id[a]))
        cand_id[a] = 0;
      else
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d\n\n", __LINE__);
    }
    if(list_id[a]) {
      if(TrackerDelList(I, list_id[a]))
        list_id[a] = 0;
      else
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d\n\n", __LINE__);
    }
  }

  if(TrackerGetNLink(I)) {
    fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; 0!=%d\n", __LINE__,
            TrackerGetNLink(I));
  }

  /* make sure everyone was deleted */

  for(a = 0; a < N_ID; a++) {
    if(cand_id[a])
      fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d;  %d!=0\n", __LINE__,
              cand_id[a]);
    if(list_id[a])
      fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d;  %d!=0\n", __LINE__,
              list_id[a]);
  }

  /* okay, now let's mix in some iterators... */

  for(a = 0; a < N_ID; a++) {
    cand_id[a] = TrackerNewCand(I, NULL);
    list_id[a] = TrackerNewList(I, NULL);
  }

  if(TrackerGetNCand(I) != N_ID) {
    fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; %d!=%d\n", __LINE__, N_ID,
            TrackerGetNCand(I));
  }

  if(TrackerGetNList(I) != N_ID) {
    fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; %d!=%d\n", __LINE__, N_ID,
            TrackerGetNList(I));
  }

  {
    int n_link = 0;
    int list_idx, cand_idx;

    for(a = 0; a < (N_ID * N_ID); a++) {
      cand_idx = (int) (N_ID * dist(mt));
      list_idx = (int) (N_ID * dist(mt));

      if(TrackerLink(I, cand_id[cand_idx], list_id[list_idx], 0))
        n_link++;
    }

    if(n_link != TrackerGetNLink(I)) {
      fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; %d!=%d\n", __LINE__, n_link,
              TrackerGetNLink(I));
    }
  }

  /* do iters iterate over the expected number of members? */

  {
    int len;
    int list_idx, cand_idx;

    for(a = 0; a < N_ID; a++) {
      cand_idx = (int) (N_ID * dist(mt));
      if(!(iter_id[a] = TrackerNewIter(I, cand_id[cand_idx], 0))) {
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; 0==%d\n",
                __LINE__, iter_id[a]);
      }

      len = 0;
      while(TrackerIterNextListInCand(I, iter_id[a], NULL))
        len++;

      if(len != TrackerGetNListForCand(I, cand_id[cand_idx]))
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; %d!=%d\n",
                __LINE__, len, TrackerGetNListForCand(I, cand_id[cand_idx]));

      if(TrackerDelIter(I, iter_id[a])) {
        iter_id[a] = 0;
      } else {
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; \n", __LINE__);
      }
    }

    for(a = 0; a < N_ID; a++) {
      list_idx = (int) (N_ID * dist(mt));
      if(!(iter_id[a] = TrackerNewIter(I, list_id[list_idx], 0))) {
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; 0==%d\n",
                __LINE__, iter_id[a]);
      }

      len = 0;
      while(TrackerIterNextCandInList(I, iter_id[a], NULL))
        len++;

      if(len != TrackerGetNCandForList(I, list_id[list_idx]))
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; %d!=%d\n",
                __LINE__, len, TrackerGetNCandForList(I, list_id[list_idx]));

      if(TrackerDelIter(I, iter_id[a])) {
        iter_id[a] = 0;
      } else {
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; \n", __LINE__);
      }
    }
  }
  /* are iterators robust to member deletion? */

  {
    int list_idx, cand_idx;
    int b;

    for(a = 0; a < N_ID; a++) {
      cand_idx = (int) (N_ID * dist(mt));
      if(!(iter_id[a] = TrackerNewIter(I, cand_id[cand_idx], 0))) {
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; 0==%d\n",
                __LINE__, iter_id[a]);
      }
    }

    {
      int cnt = 0;
      const int expected_cnt = 6119;    /* THIS TEST IS FRAGILE -- result depends on
                                           N_ID, random seem, tests that have been run above and
                                           of course, the iterator recovery behavior */

      for(a = 0; a < N_ID; a++) {
        cand_idx = (int) (N_ID * dist(mt));
        list_idx = (int) (N_ID * dist(mt));

        if(cand_id[cand_idx]) {
          if(TrackerDelCand(I, cand_id[cand_idx]))
            cand_id[cand_idx] = 0;
          else
            fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d\n\n", __LINE__);
        }

        if(list_id[list_idx]) {
          if(TrackerDelList(I, list_id[list_idx]))
            list_id[list_idx] = 0;
          else
            fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d\n\n", __LINE__);
        }

        for(b = 0; b < N_ID; b++) {
          if(TrackerIterNextListInCand(I, iter_id[b], NULL))
            cnt++;
        }
      }

      if(cnt != expected_cnt) {
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; %d!=%d\n",
                __LINE__, expected_cnt, cnt);
      }
    }
  }

  if(TrackerGetNIter(I) != N_ID) {
    fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; %d!=%d\n", __LINE__, N_ID,
            TrackerGetNIter(I));
  }

  /* delete iters */

  for(a = 0; a < N_ID; a++) {
    if(iter_id[a]) {
      if(TrackerDelIter(I, iter_id[a])) {
        iter_id[a] = 0;
      } else {
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; \n", __LINE__);
      }
    }
  }

  /* recreate cand, lists, links */

  for(a = 0; a < N_ID; a++) {
    if(!cand_id[a])
      cand_id[a] = TrackerNewCand(I, NULL);
    if(!list_id[a])
      list_id[a] = TrackerNewList(I, NULL);
  }

  {
    int n_link = 0;
    int list_idx, cand_idx;

    for(a = 0; a < (N_ID * N_ID); a++) {
      cand_idx = (int) (N_ID * dist(mt));
      list_idx = (int) (N_ID * dist(mt));

      if(TrackerLink(I, cand_id[cand_idx], list_id[list_idx], 0))
        n_link++;
    }

  }

  /* do iters iterate over the expected number of members? */

  {
    int len;
    int list_idx, cand_idx;
    int diff;

    for(a = 0; a < N_ID; a++) {
      cand_idx = (int) (N_ID * dist(mt));
      iter_start[a] = cand_id[cand_idx];
      if(!(iter_id[a] = TrackerNewIter(I, cand_id[cand_idx], 0))) {
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; 0==%d\n",
                __LINE__, iter_id[a]);
      }
      TrackerIterNextListInCand(I, iter_id[a], NULL);
    }

    for(a = 0; a < N_ID; a++) {
      list_idx = (int) (N_ID * dist(mt));
      if(list_id[list_idx]) {
        if(TrackerDelList(I, list_id[list_idx]))
          list_id[list_idx] = 0;
      }
    }

    for(a = 0; a < N_ID; a++) {
      len = 1;
      while(TrackerIterNextListInCand(I, iter_id[a], NULL))
        len++;

      diff = abs(len - TrackerGetNListForCand(I, iter_start[a]));
      /* shouldn't vary by more than 1 */

      if(diff > 1) {
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; %d>1\n", __LINE__, diff);
      }

      if(TrackerDelIter(I, iter_id[a])) {
        iter_id[a] = 0;
      } else {
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; \n", __LINE__);
      }
    }
  }

  /* delete everyone */

  for(a = 0; a < N_ID; a++) {
    if(cand_id[a]) {
      if(TrackerDelCand(I, cand_id[a]))
        cand_id[a] = 0;
      else
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d\n\n", __LINE__);
    }
    if(list_id[a]) {
      if(TrackerDelList(I, list_id[a]))
        list_id[a] = 0;
      else
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d\n\n", __LINE__);
    }
    if(iter_id[a]) {
      if(TrackerDelIter(I, iter_id[a])) {
        iter_id[a] = 0;
      } else {
        fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; \n", __LINE__);
      }
    }
  }

  if(TrackerGetNLink(I)) {
    fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; 0!=%d\n", __LINE__,
            TrackerGetNLink(I));
  }

  if(TrackerGetNCand(I)) {
    fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; 0!=%d\n", __LINE__,
            TrackerGetNLink(I));
  }

  if(TrackerGetNList(I)) {
    fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; 0!=%d\n", __LINE__,
            TrackerGetNList(I));
  }

  if(TrackerGetNIter(I)) {
    fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d; 0!=%d\n", __LINE__,
            TrackerGetNIter(I));
  }

  /* make sure everyone was deleted */

  for(a = 0; a < N_ID; a++) {
    if(cand_id[a])
      fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d;  %d!=0\n", __LINE__,
              cand_id[a]);
    if(list_id[a])
      fprintf(stderr, "TRACKER_UNIT_TEST FAILED AT LINE %d;  %d!=0\n", __LINE__,
              list_id[a]);
  }

  TrackerFree(I);

  if(!result) {
    fprintf(stderr, "TRACKER_UNIT_TEST FAILED -- EXITING\n");
    exit(0);
  }
  printf("Tracker unit tests SUCCESSFUL!\n");
  return result;
}

#endif
