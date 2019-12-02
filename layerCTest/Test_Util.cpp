#include <string>

#include "Test.h"

#include "Util2.h"

TEST_CASE("join_to_string", "[Util]")
{
  auto ch = '!';
  auto str = pymol::join_to_string("Hello ", 42, ch);
  REQUIRE(str == "Hello 42!");
}

TEST_CASE("string_format", "[Util]")
{
  auto str = pymol::string_format("%s %d!", "Hello", 42);
  REQUIRE(str == "Hello 42!");

  auto str2 = pymol::string_format("%s %d!", std::string("Hello"), 42);
  REQUIRE(str2 == "Hello 42!");

  auto str3 = pymol::string_format("");
  REQUIRE(str3 == "");

  auto str4 = pymol::string_format("%s %d!", "ThisStringWillNotBeSmallStringOptimized", 42);
  REQUIRE(str4 == "ThisStringWillNotBeSmallStringOptimized 42!");
}
