
/*
 * (c) Schrodinger, Inc.
 */

#include"BondTypeHistory.h"
#include"MemoryDebug.h"

#define COPY_ATTR(attr_name) dest->attr_name = src->attr_name
#define COPY_ATTR_ARR_2(attr_name) dest->attr_name[0] = src->attr_name[0]; dest->attr_name[1] = src->attr_name[1]
#define COPY_ATTR_N(attr_name, N) memcpy( dest->attr_name, src->attr_name, N)

/**
 * Use the legacy `stereo` attribute to distinguish between PyMOL 1.8.1 and
 * PyMOL 2.5 PSE binary dump formats. This allows us to introduce the `symop_2`
 * attribute without a backwards incompatible format change.
 */
namespace symop_stereo_magic_helper
{
constexpr signed char magic = 1 << 6;

static void set(void*) {}
template <typename toVersion, typename = decltype(&toVersion::stereo)>
static void set(toVersion* dest)
{
  dest->stereo = magic;
}

static bool get(const void*)
{
  return true;
}
template <typename fromVersion, typename = decltype(&fromVersion::stereo)>
static bool get(const fromVersion* src)
{
  return src->stereo == magic;
}
} // namespace symop_stereo_magic_helper

static void Copy_Attr_SymOp2(const void*, void*) {}
template <typename fromVersion, typename toVersion,
    typename = decltype(&fromVersion::symop_2),
    typename = decltype(&toVersion::symop_2)>
static void Copy_Attr_SymOp2(const fromVersion* src, toVersion* dest)
{
  if (src->symop_2 && symop_stereo_magic_helper::get(src)) {
    COPY_ATTR(symop_2);
    symop_stereo_magic_helper::set(dest);
  }
}

template <typename fromVersion, typename toVersion>
void Copy_BondType(const fromVersion *src,  toVersion *dest){
  COPY_ATTR(index[0]);
  COPY_ATTR(index[1]);
  COPY_ATTR(order);
  COPY_ATTR(unique_id);
  COPY_ATTR(has_setting);
  Copy_Attr_SymOp2(src, dest);
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
  return nullptr;

}
