#include "Test.h"

#include "Crystal.h"
#include "Vector.h"

TEST_CASE("cryst-misc", "[CCrystal]")
{
  CCrystal cryst(nullptr);
  REQUIRE(cryst.dims()[0] == Approx(1.0));
  REQUIRE(cryst.dims()[1] == Approx(1.0));
  REQUIRE(cryst.dims()[2] == Approx(1.0));
  REQUIRE(cryst.angles()[0] == Approx(90.0));
  REQUIRE(cryst.angles()[1] == Approx(90.0));
  REQUIRE(cryst.angles()[2] == Approx(90.0));
  REQUIRE(is_identityf(3, cryst.fracToReal()));
  REQUIRE(is_identityf(3, cryst.realToFrac()));

  // set dims
  cryst.setDims(4, 5, 8);
  float const f2r_458[] = {4, 0, 0, 0, 5, 0, 0, 0, 8};
  float const r2f_458[] = {0.25, 0, 0, 0, 0.2, 0, 0, 0, 0.125};

  REQUIRE(cryst.dims()[0] == Approx(4.0));
  REQUIRE(cryst.dims()[1] == Approx(5.0));
  REQUIRE(cryst.dims()[2] == Approx(8.0));
  REQUIRE(cryst.angles()[0] == Approx(90.0));
  REQUIRE(cryst.angles()[1] == Approx(90.0));
  REQUIRE(cryst.angles()[2] == Approx(90.0));
  REQUIRE(is_allclosef(3, cryst.fracToReal(), 3, f2r_458, 3));
  REQUIRE(is_allclosef(3, cryst.realToFrac(), 3, r2f_458, 3));

  // set angles
  cryst.setAngles(60, 70, 80);
  float const f2r_458_607080[] = {
      4., 0.8682, 2.7362, 0., 4.9240, 3.5792, 0., 0., 6.6108};
  float const r2f_458_607080[] = {
      0.25, -0.0441, -0.0796, 0., 0.2031, -0.11, 0., 0., 0.1513};

  REQUIRE(cryst.dims()[0] == Approx(4.0));
  REQUIRE(cryst.dims()[1] == Approx(5.0));
  REQUIRE(cryst.dims()[2] == Approx(8.0));
  REQUIRE(cryst.angles()[0] == Approx(60.0));
  REQUIRE(cryst.angles()[1] == Approx(70.0));
  REQUIRE(cryst.angles()[2] == Approx(80.0));
  REQUIRE(is_allclosef(3, cryst.fracToReal(), 3, f2r_458_607080, 3, 1e-3));
  REQUIRE(is_allclosef(3, cryst.realToFrac(), 3, r2f_458_607080, 3, 1e-3));

  // volume
  REQUIRE(cryst.unitCellVolume() == Approx(130.20694));

  // cross-check
  REQUIRE(!is_identityf(3, cryst.fracToReal()));
  REQUIRE(!is_allclosef(3, cryst.fracToReal(), 3, r2f_458_607080, 3, 1e-3));

  // copy
  auto crystcopy = cryst;
  REQUIRE(crystcopy.dims()[0] == Approx(4.0));
  REQUIRE(crystcopy.dims()[1] == Approx(5.0));
  REQUIRE(crystcopy.dims()[2] == Approx(8.0));
  REQUIRE(crystcopy.angles()[0] == Approx(60.0));
  REQUIRE(crystcopy.angles()[1] == Approx(70.0));
  REQUIRE(crystcopy.angles()[2] == Approx(80.0));
  REQUIRE(is_allclosef(3, crystcopy.fracToReal(), 3, f2r_458_607080, 3, 1e-3));
  REQUIRE(is_allclosef(3, crystcopy.realToFrac(), 3, r2f_458_607080, 3, 1e-3));

  // reset
  cryst = CCrystal(nullptr);
  REQUIRE(cryst.dims()[0] == Approx(1.0));
  REQUIRE(cryst.dims()[1] == Approx(1.0));
  REQUIRE(cryst.dims()[2] == Approx(1.0));
  REQUIRE(cryst.angles()[0] == Approx(90.0));
  REQUIRE(cryst.angles()[1] == Approx(90.0));
  REQUIRE(cryst.angles()[2] == Approx(90.0));
  REQUIRE(is_identityf(3, cryst.fracToReal()));
  REQUIRE(is_identityf(3, cryst.realToFrac()));

  // set frac-to-real
  cryst.setFracToReal(f2r_458_607080);
  REQUIRE(cryst.dims()[0] == Approx(4.0));
  REQUIRE(cryst.dims()[1] == Approx(5.0));
  REQUIRE(cryst.dims()[2] == Approx(8.0));
  REQUIRE(cryst.angles()[0] == Approx(60.0));
  REQUIRE(cryst.angles()[1] == Approx(70.0));
  REQUIRE(cryst.angles()[2] == Approx(80.0));
  REQUIRE(is_allclosef(3, cryst.fracToReal(), 3, f2r_458_607080, 3, 1e-3));
  REQUIRE(is_allclosef(3, cryst.realToFrac(), 3, r2f_458_607080, 3, 1e-3));
}

// vi:sw=2:expandtab
