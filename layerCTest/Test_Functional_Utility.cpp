#include <tuple>

#include "pymol/tuple.h"
#include "pymol/functional.h"
#include "Test.h"

TEST_CASE("Integer Sequence", "[FuncUtil]")
{
  pymol::integer_sequence<std::size_t, 1, 2, 3, 4> intSeq;
  static_assert(intSeq.size() == 4, "Incorrect integer sequence size");
  REQUIRE(intSeq.size() == 4);
}

TEST_CASE("Make Integer Sequence", "[FuncUtil]")
{
  auto intSeq = pymol::make_integer_sequence<std::size_t, 10>();
  static_assert(intSeq.size() == 10, "Incorrect integer sequence size");
  REQUIRE(intSeq.size() == 10);
}

TEST_CASE("Index Sequence", "[FuncUtil]")
{
  pymol::index_sequence<0, 1, 2, 3, 4, 5> idxSeq;
  static_assert(idxSeq.size() == 6, "Incorrect integer sequence size");
  REQUIRE(idxSeq.size() == 6);
}

TEST_CASE("Make Index Sequence", "[FuncUtil]")
{
  auto idxSeq = pymol::make_index_sequence<10>();
  static_assert(idxSeq.size() == 10, "Incorrect integer sequence size");
  REQUIRE(idxSeq.size() == 10);
}

TEST_CASE("Index Sequence For", "[FuncUtil]")
{
  auto idxSeq = pymol::index_sequence_for<int, const char*, float&>();
  static_assert(idxSeq.size() == 3, "Incorrect integer sequence size");
  REQUIRE(idxSeq.size() == 3);
}

constexpr int constexprFunc()
{
  return 42;
};

TEST_CASE("Invoke", "[FuncUtil]")
{
  auto returnValue = constexprFunc();
  auto returnValue2 = pymol::invoke(constexprFunc);
  static_assert(std::is_same<decltype(returnValue), int>::value, "");
  static_assert(std::is_same<decltype(returnValue2), int>::value, "");
  static_assert(constexprFunc() == 42, "Incorrect result");
  REQUIRE(std::is_same<decltype(returnValue), int>::value);
  REQUIRE(std::is_same<decltype(returnValue2), int>::value);
  REQUIRE(constexprFunc() == 42);
  REQUIRE(pymol::invoke(constexprFunc) == 42);
}

int applyme(std::string& s, int, float)
{
  return 42;
}

TEST_CASE("Apply", "[FuncUtil]")
{
  std::string s{"Hello World"};
  auto value = pymol::apply(applyme, std::forward_as_tuple(s, 100, 3.14f));
  REQUIRE(value == 42);
}
