#include "Test.h"

#include "pymol/zstring_view.h"

#include <map>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

using namespace pymol::test;

TEST_CASE("special member functions", "[zstring_view]")
{
  std::string s1 = "foobar";
  std::string s2 = s1; // same value, different data() pointer
  std::string s3 = "other value";

  // view from std::string
  pymol::zstring_view v1(s1);
  REQUIRE(v1 == s1);
  REQUIRE(s1 == v1); // swap lhs rhs
  REQUIRE(v1 == s2);
  REQUIRE(v1 != s3);
  REQUIRE(v1 == "foobar");
  REQUIRE(v1 != "abc");

  // view from view
  pymol::zstring_view v2(v1);
  REQUIRE(v2 == v1);
  REQUIRE(v2 == s1);
  REQUIRE(v2 == s2);
  REQUIRE(v2 != s3);

  // view assigment
  v1 = v2;
  REQUIRE(v1 == v2);
  REQUIRE(v1 == "foobar");

  // string assignment
  v1 = s3;
  REQUIRE(v1 != v2);
  REQUIRE(v1 != s1);
  REQUIRE(v1 != s2);
  REQUIRE(v1 == s3);

  // static string constant assingment
  v1 = "constant";
  REQUIRE(v1 == "constant");
  REQUIRE(v1 == std::string("constant"));

  // swap
  std::swap(v1, v2);
  REQUIRE(v1 == "foobar");
  REQUIRE(v2 == "constant");

  // view to string
  auto s4 = std::string(v1.c_str());
  REQUIRE(s4 == "foobar");
}

TEST_CASE("nullsafe", "[zstring_view]")
{
  std::string s3 = "other value";
  const char* c1 = nullptr;

  auto nullsafe1 = pymol::null_safe_zstring_view(c1);
  auto nullsafe2 = pymol::null_safe_zstring_view(nullsafe1);
  auto nullsafe3 = pymol::null_safe_zstring_view(s3);
  REQUIRE(nullsafe1.c_str() != nullptr);
  REQUIRE(nullsafe1.c_str() == std::string(""));
  REQUIRE(nullsafe1 == "");
  REQUIRE("" == nullsafe1); // swap lhs rhs
  REQUIRE(nullsafe1 == std::string(""));
  REQUIRE(nullsafe1 == nullsafe2);
  REQUIRE(nullsafe3 == s3);
  REQUIRE(nullsafe3 == std::string(s3).c_str());
  REQUIRE(nullsafe3 != nullsafe1);

  // assign view from null-safe
  pymol::zstring_view v1 = nullsafe1;
  v1 = nullsafe1;

  // construct view from null-safe
  auto v4 = pymol::zstring_view(nullsafe3);
  REQUIRE(v4 == s3);
}

TEST_CASE("substr", "[zstring_view]")
{
  pymol::zstring_view v1 = "foobar";
  pymol::zstring_view v2 = v1;
  v1.remove_prefix(3);
  REQUIRE(v1 == "bar");
  REQUIRE(v2.substr(3) == "bar");
  REQUIRE(v2 == "foobar");
}

TEST_CASE("size", "[zstring_view]")
{
  auto v1 = pymol::zstring_view("foo");
  REQUIRE(v1.size() == 3);
  REQUIRE(!v1.empty());
  v1 = "";
  REQUIRE(v1.size() == 0);
  REQUIRE(v1.empty());
}

TEST_CASE("sets", "[zstring_view]")
{
  std::set<pymol::zstring_view> set1;
  std::unordered_set<pymol::zstring_view> set2;

  pymol::zstring_view v1 =
      "some very long string which can't be small string optimized";
  std::string s1 = v1.c_str(); // copy
  pymol::null_safe_zstring_view nullsafe1 = s1;

  REQUIRE(v1 == s1);
  REQUIRE(v1 == nullsafe1);
  REQUIRE(s1 == nullsafe1);

  set1.insert(v1);
  set1.insert(s1);
  set1.insert(nullsafe1);

  set2.insert(v1);
  set2.insert(s1);
  set2.insert(nullsafe1);

  set1.insert("foo");
  set2.insert("foo");

  std::string s2("foo");
  set1.insert(s2);
  set2.insert(s2);

  REQUIRE(set1.size() == 2);
  REQUIRE(set2.size() == 2);

  REQUIRE(set1.count("foo") == 1);
  REQUIRE(set2.count("foo") == 1);

  REQUIRE(set1.count("bar") == 0);
  REQUIRE(set2.count("bar") == 0);

  REQUIRE(set1.count(s1) == 1);
  REQUIRE(set2.count(s1) == 1);
}

TEST_CASE("maps", "[zstring_view]")
{
  std::map<pymol::zstring_view, int> map1;
  std::unordered_map<pymol::zstring_view, int> map2;

  map1["foo"] = 12;
  map2["foo"] = 12;
  REQUIRE(map1["foo"] == 12);
  REQUIRE(map2["foo"] == 12);
  REQUIRE(map1.size() == 1);
  REQUIRE(map2.size() == 1);

  auto s1 = std::string("foo");

  REQUIRE(map1[s1] == 12);
  REQUIRE(map2[s1] == 12);
  map1[s1] = 34;
  map2[s1] = 34;
  REQUIRE(map1.size() == 1);
  REQUIRE(map2.size() == 1);
  REQUIRE(map1["foo"] == 34);
  REQUIRE(map2["foo"] == 34);

  auto s2 = std::string("bar");

  REQUIRE(map1[s2] == 0);
  REQUIRE(map2[s2] == 0);
  map1[s2] = 56;
  map2[s2] = 56;
  REQUIRE(map1.size() == 2);
  REQUIRE(map2.size() == 2);
  REQUIRE(map1[s2] == 56);
  REQUIRE(map2[s2] == 56);
}

TEST_CASE("ostream", "[zstring_view]")
{
  pymol::zstring_view v1("foobar");
  std::ostringstream out;
  out << v1;
  REQUIRE(out.str() == "foobar");
}

TEST_CASE("copy", "[zstring_view]")
{
  pymol::zstring_view v1 = "foobar";
  char buffer[16];
  buffer[v1.copy(buffer, sizeof(buffer) - 1)] = 0;
  REQUIRE(v1 == buffer);
  buffer[v1.copy(buffer, 3)] = 0;
  REQUIRE(std::string(v1.c_str(), 3) == buffer);
  buffer[v1.copy(buffer, sizeof(buffer) - 1, 3)] = 0;
  REQUIRE(std::string(v1.c_str() + 3) == buffer);
}

TEST_CASE("starts or ends with", "[zstring_view]")
{
  pymol::zstring_view v1 = "foobar";
  REQUIRE(v1.starts_with("foo"));
  REQUIRE(!v1.starts_with("bar"));
  REQUIRE(!v1.ends_with("foo"));
  REQUIRE(v1.ends_with("bar"));
}

TEST_CASE("find", "[zstring_view]")
{
  pymol::zstring_view v1 = "foobar";
  REQUIRE(v1.find('a') == 4);
  REQUIRE(v1.find("ob") == 2);
  REQUIRE(v1.find("obo") == pymol::zstring_view::npos);
  REQUIRE(v1.find_first_of('b') == 3);
  REQUIRE(v1.find_first_of("abc") == 3);
  REQUIRE(v1.find_first_not_of("abc") == 0);
  REQUIRE(v1.find_first_not_of("abfor") == pymol::zstring_view::npos);
}

TEST_CASE("explicit bool cast", "[zstring_view]")
{
  // default constructed is false
  pymol::zstring_view v1;
  REQUIRE(!v1);
  // nullsafe always true
  pymol::null_safe_zstring_view v2;
  REQUIRE(v2);

  // empty string is true
  v1 = "";
  REQUIRE(v1);
  v2 = "";
  REQUIRE(v2);

  // null assignment is legal
  v1 = nullptr;
  REQUIRE(!v1);
  v2 = nullptr;
  REQUIRE(v2);

  v1 = "foo";
  REQUIRE(v1);
  v2 = "foo";
  REQUIRE(v2);
}

TEST_CASE("implicit cast", "[zstring_view]")
{
  static_assert(std::is_convertible<const char*, pymol::zstring_view>::value,
      "from const char*");
  static_assert(!std::is_convertible<pymol::zstring_view, const char*>::value,
      "to const char*");
  static_assert(std::is_convertible<pymol::null_safe_zstring_view,
                    pymol::zstring_view>::value,
      "from null_safe_zstring_view");
  static_assert(std::is_convertible<pymol::zstring_view,
                    pymol::null_safe_zstring_view>::value,
      "to null_safe_zstring_view");
  static_assert(std::is_convertible<std::string, pymol::zstring_view>::value,
      "from std::string");
  static_assert(!std::is_convertible<pymol::zstring_view, std::string>::value,
      "to std::string");
  static_assert(
      !std::is_convertible<bool, pymol::zstring_view>::value, "from bool");
  static_assert(
      !std::is_convertible<pymol::zstring_view, bool>::value, "to bool");
}
