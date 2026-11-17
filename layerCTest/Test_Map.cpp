#include <algorithm>
#include <set>

#include "Test.h"

#include "Map.h"
#include "PyMOLGlobals.h"

using namespace pymol;

namespace
{
// two close points and one far away
constexpr float s_xyz[][3] = {
    {0.f, 0.f, 0.f},
    {1.f, 0.f, 0.f},
    {20.f, 20.f, 20.f},
};
const float* const s_coords = &s_xyz[0][0];
constexpr int s_ncoords = 3;
constexpr float s_range = 2.f;

std::set<int> neighbors(MapType& map, const float* v)
{
  std::set<int> found;
  for (const auto j : MapEIter(map, v)) {
    found.insert(j);
  }
  return found;
}

struct InterruptGuard {
  PyMOLGlobals* G;
  explicit InterruptGuard(PyMOLGlobals* G_)
      : G(G_)
  {
    G->Interrupt = true;
  }
  ~InterruptGuard() { G->Interrupt = false; }
};
} // namespace

TEST_CASE("MapEIter finds neighbors", "[Map]")
{
  PyMOLInstance pymol;
  auto G = pymol.G();
  MapType map(G, s_range, s_coords, s_ncoords);
  REQUIRE(neighbors(map, s_coords) == std::set<int>{0, 1});
  REQUIRE(neighbors(map, s_coords + 6) == std::set<int>{2});
}

// https://github.com/schrodinger/pymol-open-source/issues/533
TEST_CASE("MapEIter ignores interrupt flag", "[Map]")
{
  PyMOLInstance pymol;
  auto G = pymol.G();
  MapType map(G, s_range, s_coords, s_ncoords);
  InterruptGuard interrupt(G);
  REQUIRE(neighbors(map, s_coords) == std::set<int>{0, 1});
  REQUIRE(neighbors(map, s_coords + 6) == std::set<int>{2});
}

// https://github.com/schrodinger/pymol-open-source/issues/533
TEST_CASE("MapSetupExpress interrupted leaves no neighbors", "[Map]")
{
  PyMOLInstance pymol;
  auto G = pymol.G();
  MapType map(G, s_range, s_coords, s_ncoords);

  {
    InterruptGuard interrupt(G);
    REQUIRE(!MapSetupExpress(&map));
  }

  // no starts into a list which was never stored
  REQUIRE(map.EList.empty());
  REQUIRE(
      map.EHead.size() == std::size_t(map.Dim[0]) * map.Dim[1] * map.Dim[2]);
  REQUIRE(std::all_of(
      map.EHead.begin(), map.EHead.end(), [](int i) { return i == 0; }));

  // lazy setup in MapEIter recovers
  REQUIRE(neighbors(map, s_coords) == std::set<int>{0, 1});
}

TEST_CASE("MapSetupExpress non-interruptible", "[Map]")
{
  PyMOLInstance pymol;
  auto G = pymol.G();
  MapType map(G, s_range, s_coords, s_ncoords);

  {
    InterruptGuard interrupt(G);
    REQUIRE(MapSetupExpress(&map, false));
  }

  REQUIRE(!map.EList.empty());
  REQUIRE(neighbors(map, s_coords) == std::set<int>{0, 1});
}
