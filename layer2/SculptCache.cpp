
/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
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

#include <unordered_map>

#include"os_predef.h"
#include"SculptCache.h"

struct SculptCacheKey {
  int rest_type;
  int id0, id1, id2, id3;

  bool operator==(const SculptCacheKey& other) const noexcept
  {
    return rest_type == other.rest_type && id0 == other.id0 &&
           id1 == other.id1 && id2 == other.id2 && id3 == other.id3;
  }

  struct Hash {
    size_t operator()(const SculptCacheKey& key) const noexcept
    {
      // Notes: id0 - id3 are unique atom indices which for typical use cases
      // (sculpting of small systems) don't exceed 16bit. Also, id2 is zero
      // for distance restraints, and id3 is zero for distance and angle
      // restraints. rest_type is currently <= 10.
      constexpr auto N = sizeof(size_t);
      return size_t(key.id0) << (4 * N) ^ size_t(key.id1) ^
             size_t(key.id2) << (6 * N) ^ size_t(key.id3) << (2 * N) ^
             size_t(key.id2) >> (2 * N) ^ size_t(key.rest_type) << (3 * N);
    }
  };
};

struct _CSculptCache {
  std::unordered_map<SculptCacheKey, float, SculptCacheKey::Hash> m_data;
};

int SculptCacheInit(PyMOLGlobals * G)
{
  G->SculptCache = new CSculptCache();
  return 1;
}

void SculptCachePurge(PyMOLGlobals * G)
{
  CSculptCache *I = G->SculptCache;
  I->m_data.clear();
}

void SculptCacheFree(PyMOLGlobals * G)
{
  delete G->SculptCache;
  G->SculptCache = nullptr;
}

int SculptCacheQuery(PyMOLGlobals * G, int rest_type, int id0, int id1, int id2, int id3,
                     float *value)
{
  CSculptCache *I = G->SculptCache;
  SculptCacheKey key = {rest_type, id0, id1, id2, id3};
  auto it = I->m_data.find(key);
  if (it != I->m_data.end()) {
    *value = it->second;
    return true;
  }
  return false;
}

void SculptCacheStore(PyMOLGlobals * G, int rest_type, int id0, int id1, int id2, int id3,
                      float value)
{
  CSculptCache *I = G->SculptCache;
  SculptCacheKey key = {rest_type, id0, id1, id2, id3};
  I->m_data[key] = value;
}
