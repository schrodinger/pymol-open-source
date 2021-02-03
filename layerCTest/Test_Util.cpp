#include <string>

#include "Test.h"

#include "Util.h"
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

#if 0
  // We have a static assert which won't allow empty string. Allowing empty
  // string would produce a compiler warning with clang.
  auto str3 = pymol::string_format("");
  REQUIRE(str3 == "");
#endif

  auto str4 = pymol::string_format("%s %d!", "ThisStringWillNotBeSmallStringOptimized", 42);
  REQUIRE(str4 == "ThisStringWillNotBeSmallStringOptimized 42!");
}

TEST_CASE("CleanStr", "[Util]")
{
  char str[256];
  strcpy(str, "Hello");
  UtilCleanStr(str);
  REQUIRE(strcmp(str, "Hello") == 0);

  char str2[256];
  strcpy(str2, "  Hello  ");
  UtilCleanStr(str2);
  REQUIRE(strcmp(str2, "Hello") == 0);

  char str3[256];
  strcpy(str3, "  Hello  42  ");
  UtilCleanStr(str3);
  REQUIRE(strcmp(str3, "Hello  42") == 0);
}

TEST_CASE("CleanStdStr", "[Util]")
{
  std::string str = "Hello";
  str = UtilCleanStdStr(str);
  REQUIRE(str == "Hello");

  std::string str2 = "  Hello  ";
  str2 = UtilCleanStdStr(str2);
  REQUIRE(str2 == "Hello");

  std::string str3 = "  Hello  42  ";
  str3 = UtilCleanStdStr(str3);
  REQUIRE(str3 == "Hello  42");
}

TEST_CASE("StartsWith", "[Util]")
{
  char str[256];
  strcpy(str, "_Hello");
  REQUIRE(pymol::starts_with(str, "_H"));

  std::string str2 = "_Jello";
  REQUIRE(pymol::starts_with(str2, "_J"));

  std::string str3 = "_Fello";
  REQUIRE(!pymol::starts_with(str2, "F"));
}
