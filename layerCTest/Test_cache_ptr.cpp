#include "Test.h"

#include "pymol/memory.h"

using namespace pymol::test;

TEST_CASE("cache_ptr move construct", "[cache_ptr]")
{
  pymol::cache_ptr<int> p1;

  {
    pymol::cache_ptr<int> p2(std::move(p1));
    REQUIRE(p1.get() == nullptr);
    REQUIRE(p2.get() == nullptr);
  }

  auto ptr = new int(123);
  p1.reset(ptr);

  {
    pymol::cache_ptr<int> p2(std::move(p1));
    REQUIRE(p1.get() == nullptr);
    REQUIRE(p2.get() == ptr);
  }
}

TEST_CASE("cache_ptr copy construct", "[cache_ptr]")
{
  pymol::cache_ptr<int> p1;

  {
    pymol::cache_ptr<int> p2(p1);
    REQUIRE(p1.get() == nullptr);
    REQUIRE(p2.get() == nullptr);
  }

  auto ptr = new int(123);
  p1.reset(ptr);

  {
    pymol::cache_ptr<int> p2(p1);
    REQUIRE(p1.get() == ptr);
    REQUIRE(p2.get() == nullptr);
  }
}

TEST_CASE("cache_ptr move assign", "[cache_ptr]")
{
  pymol::cache_ptr<int> p1;

  {
    pymol::cache_ptr<int> p2;
    p2 = std::move(p1);
    REQUIRE(p1.get() == nullptr);
    REQUIRE(p2.get() == nullptr);
  }

  auto ptr = new int(123);
  p1.reset(ptr);

  {
    pymol::cache_ptr<int> p2;
    p2 = std::move(p1);
    REQUIRE(p1.get() == nullptr);
    REQUIRE(p2.get() == ptr);
  }
}

TEST_CASE("cache_ptr copy assign", "[cache_ptr]")
{
  pymol::cache_ptr<int> p1;

  {
    pymol::cache_ptr<int> p2;
    p2 = p1;
    REQUIRE(p1.get() == nullptr);
    REQUIRE(p2.get() == nullptr);
  }

  auto ptr = new int(123);
  p1.reset(ptr);

  {
    pymol::cache_ptr<int> p2;
    p2 = p1;
    REQUIRE(p1.get() == ptr);
    REQUIRE(p2.get() == nullptr);
  }
}
