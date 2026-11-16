#include "Test.h"

#include "Color.h"
#include "Setting.h"

TEST_CASE("ColorSetAutoColors replaces the object color order", "[Color]")
{
  pymol::PyMOLInstance pymol;
  auto G = pymol.G();
  SettingSetGlobal_i(G, cSetting_auto_color_next, 0);
  ColorSetAutoColors({7, 8, 9});
  for (int expected : {7, 8, 9, 7}) { // wraps
    REQUIRE(ColorGetNext(G) == expected);
  }
  ColorSetAutoColors({}); // the PyMOL order, also after a longer table
  SettingSetGlobal_i(G, cSetting_auto_color_next, 0);
  REQUIRE(ColorGetNext(G) == 26);
}
