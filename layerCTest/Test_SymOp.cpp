#include "SymOp.h"
#include "Test.h"
#include "vla.h"

using pymol::SymOp;

TEST_CASE("symop-zero-initialization", "[SymOp]")
{
  SymOp* symop;
  char buffer[] = "ABCD";
  static_assert(sizeof(buffer) >= sizeof(SymOp), "");

  symop = new (buffer) SymOp();

  REQUIRE(symop->to_string() == "1_555");
  REQUIRE(symop->index == 0);
  REQUIRE(symop->x == 0);
  REQUIRE(symop->y == 0);
  REQUIRE(symop->z == 0);
}

TEST_CASE("symop-reset", "[SymOp]")
{
  SymOp symop;
  symop.reset("3_458");
  REQUIRE(symop.to_string() == "3_458");
  REQUIRE(symop.index == 2);
  REQUIRE(symop.x == -1);
  REQUIRE(symop.y == 0);
  REQUIRE(symop.z == 3);
}

TEST_CASE("symop-vla", "[SymOp]")
{
  // does zero-initialization
  pymol::vla<SymOp> symops;
  symops.resize(1);
  REQUIRE(symops[0].to_string() == "1_555");
}
