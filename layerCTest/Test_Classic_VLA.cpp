#include "MemoryDebug.h"
#include "Test.h"

TEST_CASE("VLA Classic Alloc", "[VLA_Classic]")
{
  const std::size_t size = 10;
  auto ptr = VLACalloc(int, size);
  REQUIRE(pymol::test::isArrayZero(ptr, size));
  VLAFree(ptr);
}
