
/*
 * (c) Schrodinger, Inc.
 */

#include"BondTypeHistory.h"
#include"MemoryDebug.h"

#define COPY_ATTR(attr_name) dest->attr_name = src->attr_name
#define COPY_ATTR_ARR_2(attr_name) dest->attr_name[0] = src->attr_name[0]; dest->attr_name[1] = src->attr_name[1]
#define COPY_ATTR_N(attr_name, N) memcpy( dest->attr_name, src->attr_name, N)

template <typename fromVersion, typename toVersion>
void Copy_BondType(const fromVersion *src,  toVersion *dest){
  COPY_ATTR(index[0]);
  COPY_ATTR(index[1]);
  COPY_ATTR(order);
  COPY_ATTR(id);
  COPY_ATTR(unique_id);
  COPY_ATTR(stereo);
  COPY_ATTR(has_setting);
}

template <typename fromVersion, typename toVersion>
void CopyN_BondType(const fromVersion *src,  toVersion *dest, int NBond){
  for (int a = 0; a < NBond; ++a){
    Copy_BondType(src++, dest++);
  }
}

template <typename toVersion>
toVersion *CreateAndCopyN_BondType(const BondType *src, int NBond){
  toVersion *dest = VLACalloc(toVersion, NBond);
  toVersion *origdest = dest;
  for (int a = 0; a < NBond; ++a){
    Copy_BondType(src++, dest++);
  }
  return origdest;
}

void Copy_Into_BondType_From_Version(const void *src, int bondInfo_version, BondType *Bond, int NBond){
  switch(bondInfo_version){
  case 176:
    CopyN_BondType((BondType_1_7_6*)src, Bond, NBond);
    break;
  case 177:
    CopyN_BondType((BondType_1_7_7*)src, Bond, NBond);
    break;
  case 181:
    CopyN_BondType((BondType_1_8_1*)src, Bond, NBond);
    break;
  default:
      printf("ERROR: unknown bondInfo_version=%d from BondInfoVERSION=%d\n", bondInfo_version, BondInfoVERSION);
  }
}

void *Copy_To_BondType_Version(int bondInfo_version, BondType *Bond, int NBond){
  switch (bondInfo_version){
  case 176:
    return CreateAndCopyN_BondType<BondType_1_7_6>(Bond, NBond);
    break;
  case 177:
    return CreateAndCopyN_BondType<BondType_1_7_7>(Bond, NBond);
    break;
  case 181:
    return CreateAndCopyN_BondType<BondType_1_8_1>(Bond, NBond);
    break;
  default:
      printf("ERROR: Copy_To_BondType_Version: unknown bondInfo_version=%d from BondInfoVERSION=%d\n", bondInfo_version, BondInfoVERSION);
  }
  return NULL;

}
