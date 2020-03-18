#include "Test.h"

#include "pymol/type_traits.h"

TEST_CASE("reshape", "[type_traits]")
{
  float flat[] = {1.f, 2.f, 3.f, 4.f, 5.f, 6.f};
  auto shaped = pymol::reshape<3>(flat);
  static_assert(std::is_same<decltype(shaped), float(*)[3]>::value, "");

  typedef float f3array[3];
  f3array* foo = shaped;
  static_assert(std::is_same<decltype(shaped), decltype(foo)>::value, "");

  REQUIRE(pymol::equal(shaped[0], shaped[0] + 3, flat + 0));
  REQUIRE(pymol::equal(shaped[1], shaped[1] + 3, flat + 3));
}

TEST_CASE("flatten", "[type_traits]")
{
  float shaped[2][3] = {{1.f, 2.f, 3.f}, {4.f, 5.f, 6.f}};
  const float* flat = pymol::flatten(shaped);
  REQUIRE(pymol::equal(shaped[0], shaped[0] + 3, flat + 0));
  REQUIRE(pymol::equal(shaped[1], shaped[1] + 3, flat + 3));
  auto* shaped_ptr = shaped;
  float* flat_from_ptr = pymol::flatten(shaped_ptr);
  REQUIRE(pymol::equal(flat, flat + 6, flat_from_ptr));
}
