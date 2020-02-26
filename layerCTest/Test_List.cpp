#include "List.h"
#include "ListMacros.h"
#include "Test.h"

#include <numeric>

struct Foo {
  int someInt;
  Foo* next;
};

/**
 * Mimics ListElemCalloc without PyMOLGlobals
 * @return a zero-initialized instance on the free store
 */
template <typename T> static T* ListElemCallocNoG()
{
  static_assert(std::is_trivial<T>::value, "");
  return pymol::calloc<T>(1);
}

TEST_CASE("Legacy List", "[List]")
{
  Foo* myList;
  ListInit(myList);
  for (int i = 0; i < 4; ++i) {
    auto bar = ListElemCallocNoG<Foo>();
    bar->someInt = i;
    ListAppend(myList, bar, next, Foo);
  }

  Foo* fooPtr = nullptr;
  int count = 0;
  while (ListIterate(myList, fooPtr, next)) {
    count += fooPtr->someInt;
  }
  REQUIRE(count == 6);
  ListFree(myList, next, Foo);
}


TEST_CASE("pymol::ListAdapter type checks", "[List]")
{
  Foo* myList;
  ListInit(myList);
  auto myAdapter = pymol::make_list_adapter(myList);
  static_assert(
      std::is_same<decltype(myAdapter.begin())::value_type, Foo>::value,
      "auto begin wrong type");
  static_assert(
      std::is_same<decltype(myAdapter.cbegin())::value_type, const Foo>::value,
      "auto cbegin wrong type");
  static_assert(std::is_same<decltype(myAdapter.end())::value_type, Foo>::value,
      "auto end wrong type");
  static_assert(
      std::is_same<decltype(myAdapter.cend())::value_type, const Foo>::value,
      "auto cend wrong type");
  ListFree(myList, next, Foo);
}

TEST_CASE("pymol::ListAdapter For pre", "[List]")
{
  Foo* myList;
  ListInit(myList);
  for (int i = 0; i < 4; ++i) {
    auto bar = ListElemCallocNoG<Foo>();
    bar->someInt = i;
    ListAppend(myList, bar, next, Foo);
  }

  auto myAdapter = pymol::make_list_adapter(myList);

  int count = 0;
  for (auto it = myAdapter.begin(); it != myAdapter.end(); ++it) {
    count += it->someInt;
  }
  REQUIRE(count == 6);
  ListFree(myList, next, Foo);
}

TEST_CASE("pymol::ListAdapter For post", "[List]")
{
  Foo* myList;
  ListInit(myList);
  for (int i = 0; i < 4; ++i) {
    auto bar = ListElemCallocNoG<Foo>();
    bar->someInt = i;
    ListAppend(myList, bar, next, Foo);
  }

  auto myAdapter = pymol::make_list_adapter(myList);

  int count = 0;
  for (auto it = myAdapter.begin(); it != myAdapter.end(); it++) {
    count += it->someInt;
  }
  REQUIRE(count == 6);
  ListFree(myList, next, Foo);
}

TEST_CASE("pymol::ListAdapter For cbegin", "[List]")
{
  Foo* myList;
  ListInit(myList);
  for (int i = 0; i < 4; ++i) {
    auto bar = ListElemCallocNoG<Foo>();
    bar->someInt = i;
    ListAppend(myList, bar, next, Foo);
  }

  auto myAdapter = pymol::make_list_adapter(myList);

  int count = 0;
  for (auto it = myAdapter.cbegin(); it != myAdapter.cend(); it++) {
    count += it->someInt;
  }
  REQUIRE(count == 6);

  ListFree(myList, next, Foo);
}

TEST_CASE("pymol::ListAdapter For const_iterator", "[List]")
{
  Foo* myList;
  ListInit(myList);
  for (int i = 0; i < 4; ++i) {
    auto bar = ListElemCallocNoG<Foo>();
    bar->someInt = i;
    ListAppend(myList, bar, next, Foo);
  }

  auto myAdapter = pymol::make_list_adapter(myList);

  int count = 0;
  for (pymol::ListAdapter<Foo>::const_iterator it = myAdapter.begin();
       it != myAdapter.end(); it++) {
    count += it->someInt;
  }
  REQUIRE(count == 6);
  ListFree(myList, next, Foo);
}

TEST_CASE("pymol::ListAdapter Ranged For", "[List]")
{
  Foo* myList;
  ListInit(myList);
  for (int i = 0; i < 4; ++i) {
    auto bar = ListElemCallocNoG<Foo>();
    bar->someInt = i;
    ListAppend(myList, bar, next, Foo);
  }

  int count = 0;
  for (const auto& foo : pymol::ListAdapter<Foo>(myList)) {
    count += foo.someInt;
  }
  REQUIRE(count == 6);
  ListFree(myList, next, Foo);
}

TEST_CASE("pymol::ListAdapter Std Algorithm InputIterator", "[List]")
{
  Foo* myList;
  ListInit(myList);
  for (int i = 0; i < 4; ++i) {
    auto bar = ListElemCallocNoG<Foo>();
    bar->someInt = i;
    ListAppend(myList, bar, next, Foo);
  }

  auto myAdapter = pymol::make_list_adapter(myList);

  int count = std::accumulate(myAdapter.begin(), myAdapter.end(), 0,
      [](int sum, const Foo& b) { return sum + b.someInt; });
  REQUIRE(count == 6);
  ListFree(myList, next, Foo);
}

TEST_CASE("pymol::ListAdapter Std Algorithm OutputIterator", "[List]")
{
  Foo* myList;
  ListInit(myList);
  for (int i = 0; i < 4; ++i) {
    auto bar = ListElemCallocNoG<Foo>();
    bar->someInt = i;
    ListAppend(myList, bar, next, Foo);
  }

  auto myAdapter = pymol::make_list_adapter(myList);

  std::transform(
      myAdapter.begin(), myAdapter.end(), myAdapter.begin(), [](Foo& foo) {
        foo.someInt += 1;
        return foo;
      });

  int count = std::accumulate(myAdapter.begin(), myAdapter.end(), 0,
      [](int sum, const Foo& b) { return sum + b.someInt; });
  REQUIRE(count == 10);
  ListFree(myList, next, Foo);
}

// vi:sw=2:expandtab
