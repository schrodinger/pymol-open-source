#include <string>

#include "Test.h"

#include "pymol/algorithm.h"

TEST_CASE("Clamp", "[Algorithm]")
{
  REQUIRE(pymol::clamp(-1, 0, 10) == 0);
  REQUIRE(pymol::clamp(3, 0, 10) == 3);
  REQUIRE(pymol::clamp(13, 0, 10) == 10);
}

TEST_CASE("Equal", "[Algorithm]")
{
  std::vector<int> a{1, 2, 3, 4};
  std::vector<int> b{1, 2, 3, 4, 5};
  REQUIRE(pymol::equal(a.begin(), a.end(), b.begin()));
  b[3] = 42;
  REQUIRE(!pymol::equal(a.begin(), a.end(), b.begin()));
}

