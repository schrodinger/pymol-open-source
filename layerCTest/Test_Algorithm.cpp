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

TEST_CASE("Almost Equal", "[Algorithm]")
{
  std::vector<int> a{1, 2, 3, 4};
  std::vector<int> b{1, 2, 3, 4, 5};
  REQUIRE(pymol::equal(a.begin(), a.end(), b.begin()));
  b[3] = 42;
  REQUIRE(!pymol::equal(a.begin(), a.end(), b.begin()));
}

TEST_CASE("Range Equal", "[Algorithm]")
{
  std::array<int, 4> a{1, 2, 3, 4};
  int b[4] = {1, 2, 3, 4};
  static_assert(a.size() == (sizeof(b)/sizeof(*b)), "Arrays not equal size!");
  REQUIRE(pymol::ranges::equal(a, b));
  a[0] = 2;
  REQUIRE(!pymol::ranges::equal(a, b));
}

TEST_CASE("Contains", "[Algorithm]")
{
  std::array<int, 4> a{1, 2, 3, 4};
  REQUIRE(pymol::ranges::contains(a, 1));
  REQUIRE(!pymol::ranges::contains(a, 5));
}

TEST_CASE("Contains If", "[Algorithm]")
{
  std::array<int, 4> a{1, 2, 3, 4};
  std::array<int, 5> b{1, 2, 3, 4, 5};
  constexpr int N = 5;
  auto is_N_divisible = [N](const int i) { return i % N == 0; };
  REQUIRE(!pymol::ranges::contains_if(a, is_N_divisible));
  REQUIRE(pymol::ranges::contains_if(b, is_N_divisible));
}

TEST_CASE("Left Fold", "[Algorithm]")
{
  std::array<int, 5> a{1, 2, 3, 4, 5};
  REQUIRE(pymol::ranges::left_fold(a) == 15u);
  REQUIRE(pymol::ranges::left_fold(a, 10) == 25u);
  REQUIRE(pymol::ranges::left_fold(a, 1, std::multiplies<int>()) == 120u);

  struct S { std::size_t value; S(std::size_t v) : value{v}{} };
  std::array<S, 5> b{ S{ 1 }, S{ 2 }, S{ 3 }, S{ 4 }, S{ 5 }};
  REQUIRE(pymol::ranges::left_fold(b, 1u, [](std::size_t acc, const S& s){ return acc + s.value; }) == 16u);
}
