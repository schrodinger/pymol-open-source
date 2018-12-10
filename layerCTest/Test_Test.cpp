#include <algorithm>
#include <numeric>

#include "Test.h"

using namespace pymol::test;

namespace {
struct MyTestClass
{
  int i;
  float f;
  double d;
  bool b;
  int* p;
  int a[16];
};
}
// Creates a non-zero buffer
template <const int N> std::array<int, N> nonZeroBuffer()
{
  std::array<int, N> buffer;
  std::iota(buffer.begin(), buffer.end(), 1);
  return buffer;
}

TEST_CASE("ZeroStruct", "[TEST]")
{
  MyTestClass obj{}; // Uniform Initialization required
  REQUIRE(isStructZero(obj));
  obj.i = 1;
  REQUIRE(!isStructZero(obj));
}

TEST_CASE("isArrayZero", "[TEST]"){
  std::array<int, 10> a{}; //Uniform Initialization required
  REQUIRE(isArrayZero(a.data(), a.size()));
  a[0] = 1;
  REQUIRE(!isArrayZero(a.data(), a.size()));
}

TEST_CASE("isArrayEqual", "[TEST]")
{
  std::array<int, 4> a{1, 2, 3, 4};
  int b[4] = {1, 2, 3, 4};
  static_assert(a.size() == (sizeof(b)/sizeof(*b)), "Arrays not equal size!");
  REQUIRE(isArrayEqual(a.data(), b, a.size()));
  a[0] = 2;
  REQUIRE(!isArrayEqual(a.data(), b, a.size()));
}

TEST_CASE("isNullptr", "[TEST]")
{
  int a = 1;
  int* p = nullptr;
  REQUIRE(isNullptr(p));
  p = &a;
  REQUIRE(!isNullptr(p));
}
