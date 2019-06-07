#include "Test.h"

#include "CifFile.h"

using namespace pymol::test;

const char* SAMPLE_CIF_STR = R"""(
data_FOO
_cat1.key1 noquotes
_cat1.key2 "two words"
_cat1.key3 ? # unknown
_cat1.key4 . # inapplicable
_cat1.KEY5 "UPPER CASE key"
loop_
_cat2.key1
_cat2.key2
_cat2.key3
_cat2.key4
10 0.1 . foo
11 0.2 ? "TWO WORDS"
12  ?  ?
;multi
line
value
; . 0.4 . .
data_bar
data_baz
_undotted_key "why not"
_typed_float1 1.23(45)e3
_typed_float2 1.234(5)e1
_typed_float3 1.23456789
)""";

TEST_CASE("misc", "[CifFile]")
{
  // syntax 1
  pymol::cif_file cf1(nullptr, SAMPLE_CIF_STR);
  // syntax 2 (requires move constructor)
  auto cf2 = pymol::cif_file(nullptr, SAMPLE_CIF_STR);
  // move assign
  pymol::cif_file cf3;
  cf3 = pymol::cif_file(nullptr, SAMPLE_CIF_STR);

  // check all three instances have same data
  REQUIRE(cf1.datablocks().size() == 3);
  REQUIRE(cf2.datablocks().size() == 3);
  REQUIRE(cf3.datablocks().size() == 3);
  REQUIRE(cf1.datablocks()[2].get_opt("_undotted_key")->as_s() == std::string("why not"));
  REQUIRE(cf2.datablocks()[2].get_opt("_undotted_key")->as_s() == std::string("why not"));
  REQUIRE(cf3.datablocks()[2].get_opt("_undotted_key")->as_s() == std::string("why not"));

  auto& blocks = cf1.datablocks();

  REQUIRE(blocks[0].code() == std::string("FOO"));
  REQUIRE(blocks[1].code() == std::string("bar"));
  REQUIRE(blocks[2].code() == std::string("baz"));

  auto* data = &blocks.front();

  REQUIRE(data->get_arr("_cat1.key3") != nullptr);
  REQUIRE(data->get_arr("_cat1.key3") == data->get_opt("_cat1.key3"));
  REQUIRE(data->get_arr("_cat1.key6") == nullptr);

  REQUIRE(data->get_opt("_cat1.key1")->is_missing() == false);
  REQUIRE(data->get_opt("_cat1.key2")->is_missing() == false);
  REQUIRE(data->get_opt("_cat1.key3")->is_missing());
  REQUIRE(data->get_opt("_cat1.key4")->is_missing());
  REQUIRE(data->get_opt("_cat1.key5")->is_missing() == false);

  REQUIRE(data->get_opt("_cat1.key4")->is_missing_all());
  REQUIRE(data->get_opt("_cat1.key5")->is_missing_all() == false);

  // looped data

  REQUIRE(data->get_opt("_cat2.key1")->is_missing_all() == false);
  REQUIRE(data->get_opt("_cat2.key3")->is_missing_all());

  // template getters

  std::vector<int> vec1{10, 11, 12, 0};
  std::vector<float> vec2{0.1f, 0.2f, 99.f, 0.4f};

  REQUIRE(data->get_opt("_cat2.key1")->to_vector<int>() == vec1);
  REQUIRE(data->get_opt("_cat2.key2")->to_vector<float>(99.f) == vec2);

  REQUIRE(data->get_opt("_cat2.key4")->as<std::string>(0) == "foo");
  REQUIRE(data->get_opt("_cat2.key4")->as<std::string>(1) == "TWO WORDS");
  REQUIRE(data->get_opt("_cat2.key4")->as<std::string>(2) == "multi\nline\nvalue");
  REQUIRE(data->get_opt("_cat2.key4")->as<std::string>(3) == "");

  REQUIRE(data->get_opt("_cat2.key4")->as<const char*>(0) == std::string("foo"));
  REQUIRE(data->get_opt("_cat2.key4")->as<const char*>(3) == nullptr);

  REQUIRE(data->get_opt("_cat2.key4")->to_vector<const char*>()[0] == std::string("foo"));
  REQUIRE(data->get_opt("_cat2.key4")->to_vector<const char*>()[3] == nullptr);
  REQUIRE(data->get_opt("_cat2.key4")->to_vector<const char*>("ABC")[0] == std::string("foo"));
  REQUIRE(data->get_opt("_cat2.key4")->to_vector<const char*>("ABC")[3] == std::string("ABC"));

  // type deducted from default value

  REQUIRE(data->get_opt("_cat2.key1")->as(0, 99) / 3 == 3); // int
  REQUIRE(data->get_opt("_cat2.key1")->as(0, 99) / 3 != Approx(10. / 3.)); // int
  REQUIRE(data->get_opt("_cat2.key1")->as(0, 99.) / 3 == Approx(10. / 3.)); // double
  REQUIRE(data->get_opt("_cat2.key2")->as(0, 99.) == 0.1);
  REQUIRE(data->get_opt("_cat2.key3")->as(0, 99.f) == 99.f);
  REQUIRE(data->get_opt("_cat2.key4")->as(0, std::string("type deducted")) == "foo");
  REQUIRE(data->get_opt("_cat2.key4")->as(3, std::string("type deducted")) == "type deducted");

  // as_X getters

  REQUIRE(data->get_opt("_cat2.key4")->as_s(0, "ABC") == std::string("foo"));
  REQUIRE(data->get_opt("_cat2.key4")->as_s(3, "ABC") == std::string("ABC")); // missing

  REQUIRE(data->get_opt("_cat2.key1")->as_i(0, 99) == 10);
  REQUIRE(data->get_opt("_cat2.key1")->as_i(1, 99) == 11);
  REQUIRE(data->get_opt("_cat2.key1")->as_i(3, 99) == 99); // missing

  REQUIRE(data->get_opt("_cat2.key1")->as_d(0, 99.) == 10.);
  REQUIRE(data->get_opt("_cat2.key1")->as_d(1, 99.) == 11.);
  REQUIRE(data->get_opt("_cat2.key1")->as_d(3, 99.) == 99.);  // missing

  REQUIRE(data->get_opt("_cat2.key2")->as_d(0, 99.) == 0.1);
  REQUIRE(data->get_opt("_cat2.key2")->as_d(2, 99.) == 99.f); // missing
  REQUIRE(data->get_opt("_cat2.key2")->as_d(3, 99.) == 0.4);

  // out of bounds is default

  REQUIRE(data->get_opt("_cat2.key1")->as_i(50, 99) == 99);

  // alternate names

  REQUIRE(data->get_opt("_cat2.key1", "_other_name")->as_i(0, 99) == 10);
  REQUIRE(data->get_opt("_other_name", "_cat2.key1")->as_i(0, 99) == 10);
  REQUIRE(data->get_opt("_other_name", "_cat2_key1")->as_i(0, 99) == 99);

  // wildcard lookup

  REQUIRE(data->get_arr("_cat2_key1") == nullptr);
  REQUIRE(data->get_opt("_cat2?key1")->as_i(0, 99) == 10);
  REQUIRE(blocks[2].get_arr("_undotted.key") == nullptr);
  REQUIRE(blocks[2].get_opt("_undotted?key")->as_s() == std::string("why not"));

  // float parsing

  REQUIRE(blocks[2].get_opt("_typed_float1")->as<float>() == Approx(1230.f));
  REQUIRE(blocks[2].get_opt("_typed_float1")->as<double>() == Approx(1230.00000));
  REQUIRE(blocks[2].get_opt("_typed_float2")->as<double>() == Approx(12.3400000));
  REQUIRE(blocks[2].get_opt("_typed_float3")->as<double>() == Approx(1.23456789));
}

// vi:sw=2:expandtab
