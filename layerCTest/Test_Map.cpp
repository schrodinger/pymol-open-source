#include <cmath>
#include <random>
#include <set>

#include "Test.h"

#include "Map.h"
#include "Vector.h"

using namespace pymol;

/* The list MapSetupExpress files under a point holds every point within the
   cutoff of it, once, and ends with -1. */
static void check_neighbor_lists(PyMOLGlobals* G, int n)
{
  const float cutoff = 1.5F;
  std::mt19937 gen(n);
  std::uniform_real_distribution<float> coord(0.F, std::cbrt(float(n)));
  std::vector<float> v(3 * n);
  for (auto& x : v)
    x = coord(gen);

  MapType map(G, cutoff, v.data(), n, nullptr);
  REQUIRE(MapSetupExpress(&map));

  int missing = 0, repeated = 0;
  for (int a = 0; a < n; a += n / 100) {
    const float* va = v.data() + 3 * a;
    std::set<int> listed;
    int i = *MapLocusEStart(&map, va);
    REQUIRE(i > 0);
    for (int j = map.EList[i]; j >= 0; j = map.EList[++i])
      repeated += !listed.insert(j).second;
    for (int j = 0; j < n; j++)
      if (within3f(v.data() + 3 * j, va, cutoff))
        missing += !listed.count(j);
  }
  REQUIRE(missing == 0);
  REQUIRE(repeated == 0);
}

TEST_CASE("MapSetupExpress lists every neighbor", "[Map]")
{
  PyMOLInstance pymol;
  check_neighbor_lists(pymol.G(), 500);   // one thread
  check_neighbor_lists(pymol.G(), 20000); // above the threading threshold
}

/* MapSetupExpressPerp files under a square the points of the 27 squares
   around it; points flagged as non-spanners only from the square's own Z
   layer. */
static void check_perp_lists(PyMOLGlobals* G, int n)
{
  std::mt19937 gen(n);
  std::uniform_real_distribution<float> coord(0.F, std::cbrt(float(n)));
  std::vector<float> v(3 * n);
  for (int i = 0; i < 3 * n; i++)
    v[i] = i % 3 == 2 ? -10.F - coord(gen) : coord(gen); // in front of camera
  std::vector<int> spanner(n);
  for (int i = 0; i < n; i++)
    spanner[i] = i % 3 == 0;

  MapType map(G, 1.5F, v.data(), n, nullptr);
  REQUIRE(MapSetupExpressPerp(&map, v.data(), 1.F, true, spanner.data()));

  int wrong = 0;
  for (int a = map.iMin[0]; a <= map.iMax[0]; a++)
    for (int b = map.iMin[1]; b <= map.iMax[1]; b++)
      for (int c = map.iMin[2]; c <= map.iMax[2]; c++) {
        std::multiset<int> expected, listed;
        bool any = false;
        for (int d = a - 1; d <= a + 1; d++)
          for (int e = b - 1; e <= b + 1; e++)
            for (int f = c - 1; f <= c + 1; f++)
              for (int i = *MapFirst((&map), d, e, f); i >= 0;
                   i = MapNext((&map), i)) {
                any = true;
                if (f == c || spanner[i])
                  expected.insert(i);
              }
        int st = -*MapEStart((&map), a, b, c); // negative_start
        wrong += any != (st > 0);
        if (st > 0)
          for (int j = map.EList[st]; j >= 0; j = map.EList[++st])
            listed.insert(j);
        wrong += listed != expected;
      }
  REQUIRE(wrong == 0);
}

TEST_CASE("MapSetupExpressPerp lists the 27 squares", "[Map]")
{
  PyMOLInstance pymol;
  check_perp_lists(pymol.G(), 500);   // one thread
  check_perp_lists(pymol.G(), 20000); // above the threading threshold
}
