#include "Picking.h"
#include "Test.h"

using namespace pymol::test;

TEST_CASE("getTotalBits", "[Picking]")
{
  PickColorConverter pickconv;
  REQUIRE(pickconv.getTotalBits() == 12);

  const int bits1[] = {1, 2, 3, 4};
  pickconv.setRgbaBits(bits1);
  REQUIRE(pickconv.getTotalBits() == 10);

  const int bits2[] = {2, 3, 4, 5};
  pickconv.setRgbaBits(bits2);
  REQUIRE(pickconv.getTotalBits() == 14);

  // check pick-though bit with 32bit
  const int bits3[] = {8, 8, 8, 8};
  pickconv.setRgbaBits(bits3);
  REQUIRE(pickconv.getTotalBits() == 31);
}

TEST_CASE("validateCheckBits", "[Picking]")
{
  PickColorConverter pickconv;

  const int bits1[] = {2, 3, 4, 5};
  pickconv.setRgbaBits(bits1, 2);
  REQUIRE(pickconv.getTotalBits() == 8);

  const int bits2[] = {8, 8, 8, 4};
  pickconv.setRgbaBits(bits2, 3);
  REQUIRE(pickconv.getTotalBits() == 17);

  unsigned char color2[4] = {0xfc, 0x04, 0x0c, 0x2f};
  REQUIRE(pickconv.validateCheckBits(color2));
  color2[1] |= 0x02;
  REQUIRE(!pickconv.validateCheckBits(color2));
}

TEST_CASE("indexFromColor", "[Picking]")
{
  PickColorConverter pickconv;

  unsigned index1 = 0x816u; // 12-bit: hightest bit, lowest bit, middle bits
  unsigned char color1[4];
  pickconv.colorFromIndex(color1, index1);
  unsigned index2 = pickconv.indexFromColor(color1);
  REQUIRE(index1 == index2);

  index1 = 12345; // 14-bit
  const int bits[] = {2, 3, 4, 5};
  pickconv.setRgbaBits(bits);

  unsigned char color2[4];
  pickconv.colorFromIndex(color2, index1);
  index2 = pickconv.indexFromColor(color2);
  REQUIRE(index1 == index2);

  // with different bit depth, colors should be different
  REQUIRE((color1[0] != color2[0] || color1[1] != color2[1] ||
           color1[2] != color2[2] || color1[3] != color2[3]));
}

TEST_CASE("rounding", "[Picking]")
{
  PickColorConverter pickconv;

  unsigned index1 = 1234; // 12-bit
  unsigned char color1[4];
  pickconv.colorFromIndex(color1, index1);

  // withing rounding tolerance
  color1[0] += 1;
  color1[1] -= 1;
  color1[2] += 2;
  unsigned index2 = pickconv.indexFromColor(color1);
  REQUIRE(index1 == index2);

  // beyond rounding tolernace
  color1[0] += 1 << 4;
  unsigned index3 = pickconv.indexFromColor(color1);
  REQUIRE(index1 != index3);
}

TEST_CASE("pickmgr", "[Picking]")
{
  PickColorManager pickmgr;
  PickContext ctx{};
  unsigned char color[4];

  REQUIRE(pickmgr.count() == 0);
  pickmgr.colorNext(color, &ctx, 123, 0);
  REQUIRE(pickmgr.count() == 1);
  pickmgr.colorNext(color, &ctx, 123, 0); // repeated
  REQUIRE(pickmgr.count() == 1);
  pickmgr.colorNext(color, &ctx, 124, 0);
  pickmgr.colorNext(color, &ctx, 124, 5);
  REQUIRE(pickmgr.count() == 3);
  pickmgr.colorNext(color, &ctx, 0, cPickableNoPick); // no count increment
  REQUIRE(pickmgr.count() == 3);

  pickmgr.m_valid = true;

  REQUIRE(pickmgr.getIdentifier(0) == nullptr);
  REQUIRE(pickmgr.getIdentifier(1)->src.index == 123);
  REQUIRE(pickmgr.getIdentifier(2)->src.index == 124);
  REQUIRE(pickmgr.getIdentifier(3)->src.bond == 5);
  REQUIRE(pickmgr.getIdentifier(4) == nullptr);

  pickmgr.resetCount();

  REQUIRE(!pickmgr.m_valid);
  REQUIRE(pickmgr.count() == 0);
  REQUIRE(pickmgr.getIdentifier(1) != nullptr); // still available

  pickmgr.m_valid = true;

  pickmgr.invalidate();

  REQUIRE(!pickmgr.m_valid);
  REQUIRE(pickmgr.getIdentifier(1) == nullptr);
}
