#include <array>
#include "Test.h"
#include "vla.h"

using namespace pymol::test;
using pymol::vla;

TEST_CASE("VLA Default", "[VLA]")
{
  vla<int> myVLA;
  REQUIRE(isNullptr(myVLA.data()));
}

TEST_CASE("VLA Alloc And Size", "[VLA]")
{
  vla<int> myVLA(5);
  REQUIRE(myVLA.size() == 5);
}

TEST_CASE("VLA Default Val", "[VLA]")
{
  vla<int> myVLA(5, 0);
  REQUIRE(myVLA.size() == 5);
  REQUIRE(isArrayZero(myVLA.data(), myVLA.size()));
}

TEST_CASE("VLA Initializer List", "[VLA]")
{
  vla<int> myVLA{1, 2, 3, 4, 5};
  const std::array<int, 5> myArr{1, 2, 3, 4, 5};
  REQUIRE(isArrayEqual(myVLA.data(), myArr.data(), myVLA.size()));
}

TEST_CASE("VLA Special Member Functions", "[VLA]")
{
  REQUIRE(isRegular<vla<int>>());
}

TEST_CASE("VLA Copy Construct", "[VLA]")
{
  vla<int> myVLA{1, 2, 3, 4, 5};
  auto copyVLA = myVLA;
  REQUIRE(myVLA.size() == copyVLA.size());
  REQUIRE(isArrayEqual(myVLA.data(), copyVLA.data(), myVLA.size()));
}

TEST_CASE("VLA Copy Assign", "[VLA]")
{
  vla<int> myVLA{1, 2, 3, 4, 5};
  vla<int> copyVLA(1);
  REQUIRE(myVLA.size() != copyVLA.size());
  copyVLA = myVLA;
  REQUIRE(myVLA.size() == copyVLA.size());
  REQUIRE(isArrayEqual(myVLA.data(), copyVLA.data(), myVLA.size()));
}

TEST_CASE("VLA Move Construct", "[VLA]")
{
  vla<int> myVLA{1, 2, 3, 4, 5};
  auto copyVLA = std::move(myVLA);
  REQUIRE(copyVLA.size() == 5);
  const std::array<int, 5> myArr{1, 2, 3, 4, 5};
  REQUIRE(isArrayEqual(copyVLA.data(), myArr.data(), copyVLA.size()));
}

TEST_CASE("VLA Move Assign", "[VLA]")
{
  vla<int> myVLA{1, 2, 3, 4, 5};
  vla<int> copyVLA(3);
  copyVLA = std::move(myVLA);
  REQUIRE(copyVLA.size() == 5);
  const std::array<int, 5> myArr{1, 2, 3, 4, 5};
  REQUIRE(isArrayEqual(copyVLA.data(), myArr.data(), copyVLA.size()));
}

TEST_CASE("Vector_Resize", "[VLA]")
{
  vla<int> myVLA(5);
  REQUIRE(myVLA.size() == 5);
  myVLA.resize(3);
  REQUIRE(myVLA.size() == 3);
  vla<int> myVLA2;
  myVLA2.resize(3);
  REQUIRE(myVLA2.size() == 3);
}

TEST_CASE("FreeP", "[VLA]")
{
  vla<int> myVLA(10);
  REQUIRE(myVLA.size() == 10);
  myVLA.freeP();
  REQUIRE(myVLA.size() == 0);
  REQUIRE(isNullptr(myVLA.data()));
}

TEST_CASE("Bool Cast", "[VLA]")
{
  vla<int> myVLA;
  REQUIRE(myVLA == nullptr);
  vla<int> myVLA2(2);
  REQUIRE(myVLA2 != nullptr);
  myVLA2.freeP();
  REQUIRE(myVLA2 == nullptr);
}

TEST_CASE("Index", "[VLA]")
{
  vla<int> myVLA{1, 2, 3, 4, 5};
  REQUIRE(myVLA[2] == 3);
}

TEST_CASE("Range Based For", "[VLA]")
{
  vla<int> myVLA{0, 1, 2, 3, 4};
  std::size_t i{0u};
  for (auto& m : myVLA) {
    REQUIRE(myVLA[i] == i);
    i++;
  }
  i = 0;
  for (const auto& m : myVLA) {
    REQUIRE(myVLA[i] == i);
    i++;
  }
}

TEST_CASE("To_StdVector", "[VLA]")
{
  vla<int> myVLA{1, 2, 3, 4, 5};
  auto myStdVec = myVLA.toStdVector();
  REQUIRE(myStdVec.size() == myVLA.size());
  REQUIRE(isArrayEqual(myStdVec.data(), myVLA.data(), myStdVec.size()));
}

TEST_CASE("From_StdVector", "[VLA]")
{
  std::vector<int> myStdVec{1, 2, 3, 4, 5};
  vla<int> myVLA(myStdVec);
  REQUIRE(myStdVec.size() == myVLA.size());
  REQUIRE(isArrayEqual(myStdVec.data(), myVLA.data(), myStdVec.size()));
}

TEST_CASE("From_VLACalloc", "[VLA]")
{
  vla<int> myVLA(VLACalloc(int, 5));
  REQUIRE(myVLA.size() == 5);
  REQUIRE(isArrayZero(myVLA.data(), 5));
}

TEST_CASE("Classic_Copy", "[VLA]")
{
  vla<int> myVLA(5, 10);
  auto myVLACopy = VLACopy2(myVLA);
  REQUIRE(isArrayEqual(myVLA.data(), myVLACopy.data(), myVLA.size()));
  myVLACopy[1] = 100;
  REQUIRE(myVLA[1] != myVLACopy[1]);
}

TEST_CASE("Classic_Check2", "[VLA]")
{
  vla<int> myVLA(5);
  VLACheck2(myVLA, 10);
  REQUIRE(myVLA.size() >= 11);
}

TEST_CASE("Classic_Size2", "[VLA]")
{
  vla<int> myVLA(5);
  VLASize2(myVLA, 10);
  REQUIRE(myVLA.size() == 10);
}

TEST_CASE("Classic_Free", "[VLA]")
{
  vla<int> myVLA(5);
  VLAFreeP(myVLA);
  REQUIRE(isNullptr(myVLA.data()));
}

// vi:sw=2:expandtab
