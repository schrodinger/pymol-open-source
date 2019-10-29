#include "Test.h"

#include "pymol/memory.h"

using namespace pymol::test;

TEST_CASE("copyable_ptr move construct", "[copyable_ptr]")
{
  pymol::copyable_ptr<int> p1;

  {
    pymol::copyable_ptr<int> p2(std::move(p1));
    REQUIRE(p1.get() == nullptr);
    REQUIRE(p2.get() == nullptr);
  }

  auto ptr = new int(123);
  p1.reset(ptr);

  {
    pymol::copyable_ptr<int> p2(std::move(p1));
    REQUIRE(p1.get() == nullptr);
    REQUIRE(p2.get() == ptr);
  }
}

TEST_CASE("copyable_ptr copy construct", "[copyable_ptr]")
{
  pymol::copyable_ptr<int> p1;

  {
    pymol::copyable_ptr<int> p2(p1);
    REQUIRE(p1.get() == nullptr);
    REQUIRE(p2.get() == nullptr);
  }

  auto ptr = new int(123);
  p1.reset(ptr);

  {
    pymol::copyable_ptr<int> p2(p1);
    REQUIRE(p1.get() == ptr);
    REQUIRE(p2.get() != ptr);
    REQUIRE(p2.get() != nullptr);
    REQUIRE(*p2 == 123);
  }
}

TEST_CASE("copyable_ptr move assign", "[copyable_ptr]")
{
  pymol::copyable_ptr<int> p1;

  {
    pymol::copyable_ptr<int> p2;
    p2 = std::move(p1);
    REQUIRE(p1.get() == nullptr);
    REQUIRE(p2.get() == nullptr);
  }

  auto ptr = new int(123);
  p1.reset(ptr);

  {
    pymol::copyable_ptr<int> p2;
    p2 = std::move(p1);
    REQUIRE(p1.get() == nullptr);
    REQUIRE(p2.get() == ptr);
  }
}

TEST_CASE("copyable_ptr copy assign", "[copyable_ptr]")
{
  pymol::copyable_ptr<int> p1;

  {
    pymol::copyable_ptr<int> p2;
    p2 = p1;
    REQUIRE(p1.get() == nullptr);
    REQUIRE(p2.get() == nullptr);
  }

  auto ptr = new int(123);
  p1.reset(ptr);

  {
    pymol::copyable_ptr<int> p2;
    p2 = p1;
    REQUIRE(p1.get() == ptr);
    REQUIRE(p2.get() != ptr);
    REQUIRE(p2.get() != nullptr);
    REQUIRE(*p2 == 123);
  }
}
