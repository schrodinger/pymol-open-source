
/*
 * (c) Schrodinger, Inc.
 */

#include"AtomInfoHistory.h"
#include"MemoryDebug.h"

#define COPY_ATTR(attr_name) dest->attr_name = src->attr_name
#define COPY_ATTR_ARR_2(attr_name) dest->attr_name[0] = src->attr_name[0]; dest->attr_name[1] = src->attr_name[1]
#define COPY_ATTR_N(attr_name, N) memcpy( dest->attr_name, src->attr_name, N)
#define COPY_ATTR_S(attr_name) copy_attr_s(dest->attr_name, src->attr_name);

template <typename A>   inline float get_anisou_factor() { return 1.f; }
template <>             inline float get_anisou_factor<AtomInfoType_1_8_1>() { return 10000.f; }

template <typename D, typename S>
void AtomInfoTypeConverter::copy1(D *dest, const S *src) {
  COPY_ATTR(resv);
  COPY_ATTR(customType);
  COPY_ATTR(priority);
  COPY_ATTR(b);
  COPY_ATTR(q);
  COPY_ATTR(vdw);
  COPY_ATTR(partialCharge);
  COPY_ATTR(color);
  COPY_ATTR(id);
  COPY_ATTR(flags);
  COPY_ATTR(unique_id);
  COPY_ATTR(discrete_state);
  COPY_ATTR(elec_radius);
  COPY_ATTR(rank);
  COPY_ATTR(visRep);
  COPY_ATTR(formalCharge);
  COPY_ATTR(stereo);
  COPY_ATTR(cartoon);
  COPY_ATTR(hetatm);
  COPY_ATTR(bonded);
  COPY_ATTR(chemFlag);
  COPY_ATTR(geom);
  COPY_ATTR(valence);
  COPY_ATTR(masked);
  COPY_ATTR(protekted);
  COPY_ATTR(protons);
  COPY_ATTR(hb_donor);
  COPY_ATTR(hb_acceptor);
  COPY_ATTR(has_setting);

  COPY_ATTR_ARR_2(alt);

  COPY_ATTR_S(textType);
  COPY_ATTR_S(custom);
  COPY_ATTR_S(label);
  COPY_ATTR_S(segi);
  COPY_ATTR_S(chain);
  COPY_ATTR_S(resn);
  COPY_ATTR_S(name);

  dest->setResi(src->resv, src->getInscode());

  COPY_ATTR_N(elem, sizeof(ElemName));
  COPY_ATTR_ARR_2(ssType);

  if (src->has_anisou()){
    auto d_U = dest->get_anisou();
    auto s_U = src->get_anisou();
    if (d_U) {
      for (int i = 0; i < 6; ++i) {
        d_U[i] = s_U[i] * (get_anisou_factor<D>() / get_anisou_factor<S>());
      }
    }
  }
}

template <typename D, typename S>
void AtomInfoTypeConverter::copyN(D *dest, const S *src) {
  for (int a = 0; a < NAtom; ++a){
    copy1(dest++, src++);
  }
}

void AtomInfoTypeConverter::copy(AtomInfoType *dest, const void *src, int srcversion) {
  switch (srcversion){
  case 176:
    copyN(dest, (AtomInfoType_1_7_6*)src);
    break;
  case 177:
    copyN(dest, (AtomInfoType_1_7_7*)src);
    break;
  case 181:
    copyN(dest, (AtomInfoType_1_8_1*)src);
    break;
  default:
    printf("ERROR: Copy_Into_AtomInfoType_From_Version: unknown srcversion=%d from AtomInfoVERSION=%d\n", srcversion, AtomInfoVERSION);
  }
}

template <typename D>
D * AtomInfoTypeConverter::allocCopy(const AtomInfoType *src) {
  D * dest = VLACalloc(D, NAtom);
  copyN(dest, src);
  return dest;
}

void * AtomInfoTypeConverter::allocCopy(int destversion, const AtomInfoType *src) {
  switch (destversion){
  case 176:
    return allocCopy<AtomInfoType_1_7_6>(src);
  case 177:
    return allocCopy<AtomInfoType_1_7_7>(src);
  case 181:
    return allocCopy<AtomInfoType_1_8_1>(src);
  }
  printf("ERROR: AtomInfoTypeConverter: unknown destversion=%d from AtomInfoVERSION=%d\n", destversion, AtomInfoVERSION);
  return nullptr;
}
