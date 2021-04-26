#include <algorithm>
#include <iterator>

#include "CGO.h"
#include "ScrollBar.h"
#include "Test.h"

using namespace pymol::test;

namespace
{
ScrollBar getValuedScrollBar()
{
  ScrollBar sb{nullptr, true};
  sb.setLimits(11, 1); // m_MaxValue == 10
  return sb;
}
} // namespace

static
bool operator==(const BlockRect& lhs, const BlockRect& rhs)
{
  return lhs.top == rhs.top && lhs.left == rhs.left &&
         lhs.bottom == rhs.bottom && lhs.right == rhs.right;
}

TEST_CASE("ScrollBar Default", "[ScrollBar]")
{
  ScrollBar sb(nullptr, true);
  REQUIRE(sizeof(sb) == sizeof(ScrollBar));
}

TEST_CASE("ScrollBar getValue", "[ScrollBar]")
{
  ScrollBar sb(nullptr, true);
  REQUIRE(pymol::almost_equal(sb.getValue(), 0.0f));
}

TEST_CASE("ScrollBar setLimits & isMaxed", "[ScrollBar]")
{
  auto sb = getValuedScrollBar();
  REQUIRE(!sb.isMaxed());
  REQUIRE(pymol::almost_equal(sb.getValue(), 0.0f));
}

TEST_CASE("ScrollBar maxOut", "[ScrollBar]")
{
  auto sb = getValuedScrollBar();
  sb.maxOut();
  REQUIRE(sb.isMaxed());
  REQUIRE(pymol::almost_equal(sb.getValue(), 10.0f));
}

TEST_CASE("ScrollBar setValue (Clamped)", "[ScrollBar]")
{
  auto sb = getValuedScrollBar();
  sb.setValue(-2.0f);
  REQUIRE(pymol::almost_equal(sb.getValue(), 0.0f));
  sb.setValue(8.0f);
  REQUIRE(pymol::almost_equal(sb.getValue(), 8.0f));
  sb.setValue(12.0f);
  REQUIRE(pymol::almost_equal(sb.getValue(), 10.0f));
}

TEST_CASE("ScrollBar setValue (Unclamped)", "[ScrollBar]")
{
  auto sb = getValuedScrollBar();
  sb.setValueNoCheck(-2.0f);
  REQUIRE(pymol::almost_equal(sb.getValue(), -2.0f));
  sb.setValueNoCheck(8.0f);
  REQUIRE(pymol::almost_equal(sb.getValue(), 8.0f));
  sb.setValueNoCheck(12.0f);
  REQUIRE(pymol::almost_equal(sb.getValue(), 12.0f));
}

TEST_CASE("ScrollBar moveBy", "[ScrollBar]")
{
  auto sb = getValuedScrollBar();
  sb.moveBy(123.0f);
  REQUIRE(pymol::almost_equal(sb.getValue(), 10.0f));
}

TEST_CASE("ScrollBar setBox", "[ScrollBar]")
{
  auto sb = getValuedScrollBar();
  BlockRect br{5, 10, 15, 20};
  sb.setBox(5, 10, 15, 20);
  REQUIRE(br == sb.rect);
}

/**
 * TODO: Member functions to test
 *
 * int release(int, int, int, int) override;
 * int click(int, int, int, int) override;
 * int drag(int, int, int) override;
 * int draw(CGO*) override
 * int resphape(int, int) override
 *
 * void fill(CGO*);
 * void drawNoFill(CGO*);
 * int grabbed();
 * void drawHandle(float, CGO*);
 */

